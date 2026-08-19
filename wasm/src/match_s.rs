//! Radiative resistance, simultaneous port match, and power-wave \(S\).
//!
//! §6–8 of `docs/approximate_matched_basis.md`, plus a phenomenological
//! mutual reactance \(X(\Delta x,\Delta y)\), Kurokawa \(S\) at complex \(z_0\),
//! and an optional common complex reference \(z_c\) (§7.1). Operates on an
//! already-formed radiated-power Gram \(P_H\). Internals are f64; the Gram
//! itself stays f32.

pub const Z_REF: f64 = 50.0;
pub const EPS_Z: f64 = 1e-9;
pub const K_MAX: u32 = 200;
pub const TAU: f64 = 1e-3;

const PIVOT_FLOOR: f64 = 1e-18;
const MATCH_BETA: f64 = 0.5;

/// Result of \(P_H \to Z=R+jX \to z_0 \to S\).
pub struct MatchedS {
	#[allow(dead_code)]
	pub n: usize,
	#[allow(dead_code)]
	pub z_ref: f64,
	pub z0: Vec<f64>,
	pub z0_im: Vec<f64>,
	pub r_re: Vec<f64>,
	pub r_im: Vec<f64>,
	pub s_re: Vec<f64>,
	pub s_im: Vec<f64>,
	pub t_re: Vec<f64>,
	pub t_im: Vec<f64>,
	pub iterations: u32,
	pub residual: f64,
}

impl MatchedS {
	pub fn empty() -> Self {
		Self {
			n: 0,
			z_ref: Z_REF,
			z0: Vec::new(),
			z0_im: Vec::new(),
			r_re: Vec::new(),
			r_im: Vec::new(),
			s_re: Vec::new(),
			s_im: Vec::new(),
			t_re: Vec::new(),
			t_im: Vec::new(),
			iterations: 0,
			residual: 0.0,
		}
	}

	/// \(R = 2 Z_\mathrm{ref} P_H\), match on \(\Re(R)\), \(S\) from Hermitian \(R\).
	#[allow(dead_code)]
	pub fn from_gram(p_re: &[f32], p_im: &[f32], n: usize, z_ref: f64) -> Self {
		Self::from_gram_coupled(p_re, p_im, n, z_ref, &[], &[], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	}

	/// Same as [`from_gram`], then \(Z = R + jX(\Delta x,\Delta y)\). Per-port
	/// conjugate match unless \(\Re(z_c)>0\), in which case every port uses
	/// that common \(z_c\).
	pub fn from_gram_coupled(
		p_re: &[f32],
		p_im: &[f32],
		n: usize,
		z_ref: f64,
		x: &[f32],
		y: &[f32],
		x_nn: f64,
		alpha: f64,
		beta: f64,
		aniso: f64,
		z_common_re: f64,
		z_common_im: f64,
	) -> Self {
		let z_ref = if z_ref.is_finite() && z_ref > 0.0 {
			z_ref
		} else {
			Z_REF
		};
		if n == 0 || p_re.len() != n * n {
			return Self::empty();
		}
		let nn = n * n;
		let scale = 2.0 * z_ref;
		let mut r_re = vec![0.0f64; nn];
		let mut r_im = vec![0.0f64; nn];
		for i in 0..nn {
			r_re[i] = scale * (p_re[i] as f64);
			if p_im.len() == nn {
				r_im[i] = scale * (p_im[i] as f64);
			}
		}
		for p in 0..n {
			r_im[p * n + p] = 0.0;
		}

		let coupled = add_mutual_reactance(&mut r_im, x, y, n, x_nn, alpha, beta, aniso);
		if z_common_re.is_finite() && z_common_re > 0.0 {
			let zc_im = if z_common_im.is_finite() {
				z_common_im
			} else {
				0.0
			};
			Self::from_z_common(r_re, r_im, n, z_ref, z_common_re, zc_im, coupled)
		} else if coupled {
			Self::from_z_complex(r_re, r_im, n, z_ref)
		} else {
			Self::from_z_real(r_re, r_im, n, z_ref)
		}
	}

	/// Real \(z_0\) on \(\Re(Z)\); Hermitian Cholesky for \(S\) and \(T\).
	fn from_z_real(r_re: Vec<f64>, r_im: Vec<f64>, n: usize, z_ref: f64) -> Self {
		let nn = n * n;
		let mut z0 = vec![0.0f64; n];
		for p in 0..n {
			z0[p] = r_re[p * n + p].max(EPS_Z);
		}

		let mut iterations = 0u32;
		let mut work = vec![0.0f64; nn];
		for _ in 0..K_MAX {
			iterations += 1;
			let mut z_next = vec![0.0f64; n];
			fill_rre_plus_d(&r_re, &z0, n, &mut work);
			let ydiag = inverse_diag_spd(&mut work, n);
			for p in 0..n {
				let ypp = ydiag[p];
				let zin = if ypp.abs() < 1e-30 {
					f64::INFINITY
				} else {
					1.0 / ypp - z0[p]
				};
				z_next[p] = zin.max(EPS_Z);
			}
			let mut delta = 0.0f64;
			for p in 0..n {
				delta = delta.max((z_next[p] - z0[p]).abs());
			}
			z0 = z_next;
			if delta < TAU {
				break;
			}
		}

		fill_rre_plus_d(&r_re, &z0, n, &mut work);
		let ydiag = inverse_diag_spd(&mut work, n);
		let mut residual = 0.0f64;
		for p in 0..n {
			let ypp = ydiag[p];
			let zin = if ypp.abs() < 1e-30 {
				f64::INFINITY
			} else {
				1.0 / ypp - z0[p]
			};
			residual = residual.max((zin.max(EPS_Z) - z0[p]).abs());
		}

		let (s_re, s_im, t_re, t_im) = form_s_and_t(&r_re, &r_im, &z0, n, z_ref);
		Self {
			n,
			z_ref,
			z0,
			z0_im: vec![0.0f64; n],
			r_re,
			r_im,
			s_re,
			s_im,
			t_re,
			t_im,
			iterations,
			residual,
		}
	}

	/// Conjugate match \(z_0 = Z_\mathrm{in}^*\) and Kurokawa \(S\), \(T\).
	fn from_z_complex(r_re: Vec<f64>, r_im: Vec<f64>, n: usize, z_ref: f64) -> Self {
		let nn = n * n;
		let mut z0_re = vec![0.0f64; n];
		let mut z0_im = vec![0.0f64; n];
		for p in 0..n {
			z0_re[p] = r_re[p * n + p].max(EPS_Z);
		}

		let mut iterations = 0u32;
		let mut work_re = vec![0.0f64; nn];
		let mut work_im = vec![0.0f64; nn];
		let mut beta = 1.0f64;
		for k in 0..K_MAX {
			iterations += 1;
			if k == K_MAX / 4 {
				beta = MATCH_BETA;
			}
			fill_z_plus_d(&r_re, &r_im, &z0_re, &z0_im, n, &mut work_re, &mut work_im);
			let (ydiag_re, ydiag_im) = inverse_diag_cplx(&mut work_re, &mut work_im, n);
			let mut z_next_re = vec![0.0f64; n];
			let mut z_next_im = vec![0.0f64; n];
			let mut delta = 0.0f64;
			for p in 0..n {
				let (zin_re, zin_im) = zin_from_ypp(
					ydiag_re[p],
					ydiag_im[p],
					z0_re[p],
					z0_im[p],
				);
				// \(z^\star = Z_\mathrm{in}^*\) with \(\Re(z^\star)\ge\varepsilon_z\).
				let star_re = zin_re.max(EPS_Z);
				let star_im = -zin_im;
				z_next_re[p] = (1.0 - beta) * z0_re[p] + beta * star_re;
				z_next_im[p] = (1.0 - beta) * z0_im[p] + beta * star_im;
				if z_next_re[p] < EPS_Z {
					z_next_re[p] = EPS_Z;
				}
				delta = delta.max((star_re - z0_re[p]).hypot(star_im - z0_im[p]));
			}
			z0_re = z_next_re;
			z0_im = z_next_im;
			if delta < TAU {
				break;
			}
		}

		fill_z_plus_d(&r_re, &r_im, &z0_re, &z0_im, n, &mut work_re, &mut work_im);
		let (ydiag_re, ydiag_im) = inverse_diag_cplx(&mut work_re, &mut work_im, n);
		let mut residual = 0.0f64;
		for p in 0..n {
			let (zin_re, zin_im) = zin_from_ypp(
				ydiag_re[p],
				ydiag_im[p],
				z0_re[p],
				z0_im[p],
			);
			let star_re = zin_re.max(EPS_Z);
			let star_im = -zin_im;
			residual = residual.max((star_re - z0_re[p]).hypot(star_im - z0_im[p]));
		}

		let (s_re, s_im, t_re, t_im) =
			form_s_and_t_kurokawa(&r_re, &r_im, &z0_re, &z0_im, n, z_ref);
		Self {
			n,
			z_ref,
			z0: z0_re,
			z0_im,
			r_re,
			r_im,
			s_re,
			s_im,
			t_re,
			t_im,
			iterations,
			residual,
		}
	}

	/// Skip the solver: every port uses the same complex \(z_c\) (\(\Re(z_c)>0\)).
	fn from_z_common(
		r_re: Vec<f64>,
		r_im: Vec<f64>,
		n: usize,
		z_ref: f64,
		zc_re: f64,
		zc_im: f64,
		coupled: bool,
	) -> Self {
		let nn = n * n;
		let zc_re = zc_re.max(EPS_Z);
		let z0 = vec![zc_re; n];
		let z0_im = vec![zc_im; n];
		let use_kurokawa = coupled || zc_im != 0.0;

		let mut residual = 0.0f64;
		let (s_re, s_im, t_re, t_im) = if use_kurokawa {
			let mut work_re = vec![0.0f64; nn];
			let mut work_im = vec![0.0f64; nn];
			fill_z_plus_d(&r_re, &r_im, &z0, &z0_im, n, &mut work_re, &mut work_im);
			let (ydiag_re, ydiag_im) = inverse_diag_cplx(&mut work_re, &mut work_im, n);
			for p in 0..n {
				let (zin_re, zin_im) =
					zin_from_ypp(ydiag_re[p], ydiag_im[p], z0[p], z0_im[p]);
				residual = residual.max((zin_re - zc_re).hypot(zin_im - zc_im));
			}
			form_s_and_t_kurokawa(&r_re, &r_im, &z0, &z0_im, n, z_ref)
		} else {
			let mut work = vec![0.0f64; nn];
			fill_rre_plus_d(&r_re, &z0, n, &mut work);
			let ydiag = inverse_diag_spd(&mut work, n);
			for p in 0..n {
				let ypp = ydiag[p];
				let zin = if ypp.abs() < 1e-30 {
					f64::INFINITY
				} else {
					1.0 / ypp - z0[p]
				};
				residual = residual.max((zin - zc_re).abs());
			}
			form_s_and_t(&r_re, &r_im, &z0, n, z_ref)
		};

		Self {
			n,
			z_ref,
			z0,
			z0_im,
			r_re,
			r_im,
			s_re,
			s_im,
			t_re,
			t_im,
			iterations: 0,
			residual,
		}
	}
}

/// \(X_{pq}=X_{nn}(d_\min/\rho)^\alpha\cos(\beta(\rho/d_\min-1))(1+A\cos 2\varphi)\).
/// Returns true if any off-diagonal was written.
fn add_mutual_reactance(
	r_im: &mut [f64],
	x: &[f32],
	y: &[f32],
	n: usize,
	x_nn: f64,
	alpha: f64,
	beta: f64,
	aniso: f64,
) -> bool {
	if n < 2 || x.len() != n || y.len() != n {
		return false;
	}
	if !x_nn.is_finite() || x_nn == 0.0 {
		return false;
	}
	let alpha = if alpha.is_finite() && alpha >= 0.0 {
		alpha
	} else {
		return false;
	};
	let beta = if beta.is_finite() { beta } else { 0.0 };
	let aniso = if aniso.is_finite() { aniso } else { 0.0 };

	let mut d_min = f64::INFINITY;
	for p in 0..n {
		for q in 0..p {
			let dx = x[p] as f64 - x[q] as f64;
			let dy = y[p] as f64 - y[q] as f64;
			let rho = dx.hypot(dy);
			if rho > 0.0 && rho < d_min {
				d_min = rho;
			}
		}
	}
	if !d_min.is_finite() || d_min <= 0.0 {
		return false;
	}

	let mut wrote = false;
	for p in 0..n {
		for q in 0..p {
			let dx = x[p] as f64 - x[q] as f64;
			let dy = y[p] as f64 - y[q] as f64;
			let rho = dx.hypot(dy);
			let xpq = if rho > 0.0 {
				let envelope = x_nn * (d_min / rho).powf(alpha);
				let osc = (beta * (rho / d_min - 1.0)).cos();
				let phi = dy.atan2(dx);
				let loc = 1.0 + aniso * (2.0 * phi).cos();
				envelope * osc * loc
			} else {
				0.0
			};
			if xpq == 0.0 {
				continue;
			}
			r_im[p * n + q] += xpq;
			r_im[q * n + p] += xpq;
			wrote = true;
		}
	}
	wrote
}

fn zin_from_ypp(ypp_re: f64, ypp_im: f64, z0_re: f64, z0_im: f64) -> (f64, f64) {
	if ypp_re.hypot(ypp_im) < 1e-30 {
		return (f64::INFINITY, 0.0);
	}
	let (inv_re, inv_im) = cdiv(1.0, 0.0, ypp_re, ypp_im);
	(inv_re - z0_re, inv_im - z0_im)
}

fn fill_rre_plus_d(r_re: &[f64], z0: &[f64], n: usize, out: &mut [f64]) {
	out.copy_from_slice(r_re);
	for p in 0..n {
		out[p * n + p] += z0[p];
	}
}

fn fill_z_plus_d(
	z_re: &[f64],
	z_im: &[f64],
	z0_re: &[f64],
	z0_im: &[f64],
	n: usize,
	out_re: &mut [f64],
	out_im: &mut [f64],
) {
	out_re.copy_from_slice(z_re);
	out_im.copy_from_slice(z_im);
	for p in 0..n {
		out_re[p * n + p] += z0_re[p];
		out_im[p * n + p] += z0_im[p];
	}
}

#[inline]
fn cmul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
	(ar * br - ai * bi, ar * bi + ai * br)
}

#[inline]
fn cdiv(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
	let d = br * br + bi * bi;
	if d < 1e-30 {
		return (0.0, 0.0);
	}
	((ar * br + ai * bi) / d, (ai * br - ar * bi) / d)
}

/// In-place real Cholesky: lower triangle of `a` becomes \(L\) with \(A = L L^T\).
fn chol_factor_real(a: &mut [f64], n: usize) {
	for i in 0..n {
		for j in 0..=i {
			let mut sum = a[i * n + j];
			for k in 0..j {
				sum -= a[i * n + k] * a[j * n + k];
			}
			if i == j {
				a[i * n + i] = sum.max(PIVOT_FLOOR).sqrt();
			} else {
				a[i * n + j] = sum / a[j * n + j];
			}
		}
	}
}

fn chol_solve_real(l: &[f64], n: usize, b: &mut [f64]) {
	for i in 0..n {
		let mut s = b[i];
		for k in 0..i {
			s -= l[i * n + k] * b[k];
		}
		b[i] = s / l[i * n + i];
	}
	for i in (0..n).rev() {
		let mut s = b[i];
		for k in i + 1..n {
			s -= l[k * n + i] * b[k];
		}
		b[i] = s / l[i * n + i];
	}
}

/// Diagonal of \(A^{-1}\) for real SPD `a` (destroyed).
fn inverse_diag_spd(a: &mut [f64], n: usize) -> Vec<f64> {
	chol_factor_real(a, n);
	let mut ydiag = vec![0.0f64; n];
	let mut rhs = vec![0.0f64; n];
	for j in 0..n {
		rhs.fill(0.0);
		rhs[j] = 1.0;
		chol_solve_real(a, n, &mut rhs);
		ydiag[j] = rhs[j];
	}
	ydiag
}

fn chol_factor_herm(re: &mut [f64], im: &mut [f64], n: usize) {
	for i in 0..n {
		for j in 0..=i {
			let mut sum_re = re[i * n + j];
			let mut sum_im = im[i * n + j];
			for k in 0..j {
				let lir = re[i * n + k];
				let lii = im[i * n + k];
				let ljr = re[j * n + k];
				let lji = im[j * n + k];
				sum_re -= lir * ljr + lii * lji;
				sum_im -= lii * ljr - lir * lji;
			}
			if i == j {
				re[i * n + i] = sum_re.max(PIVOT_FLOOR).sqrt();
				im[i * n + i] = 0.0;
			} else {
				let d = re[j * n + j];
				re[i * n + j] = sum_re / d;
				im[i * n + j] = sum_im / d;
			}
		}
	}
}

fn chol_solve_herm(l_re: &[f64], l_im: &[f64], n: usize, b_re: &mut [f64], b_im: &mut [f64]) {
	for i in 0..n {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for k in 0..i {
			let lr = l_re[i * n + k];
			let li = l_im[i * n + k];
			sr -= lr * b_re[k] - li * b_im[k];
			si -= lr * b_im[k] + li * b_re[k];
		}
		let d = l_re[i * n + i];
		b_re[i] = sr / d;
		b_im[i] = si / d;
	}
	for i in (0..n).rev() {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for k in i + 1..n {
			let lr = l_re[k * n + i];
			let li = l_im[k * n + i];
			sr -= lr * b_re[k] + li * b_im[k];
			si -= lr * b_im[k] - li * b_re[k];
		}
		let d = l_re[i * n + i];
		b_re[i] = sr / d;
		b_im[i] = si / d;
	}
}

/// Complex LU with partial pivoting. `re`/`im` become \(L+U\) (unit \(L\) diagonal).
/// `piv[k]` is the row swapped with \(k\) at step \(k\).
fn lu_factor_cplx(re: &mut [f64], im: &mut [f64], n: usize, piv: &mut [usize]) {
	for i in 0..n {
		piv[i] = i;
	}
	for k in 0..n {
		let mut p = k;
		let mut max_abs = re[k * n + k].hypot(im[k * n + k]);
		for i in (k + 1)..n {
			let a = re[i * n + k].hypot(im[i * n + k]);
			if a > max_abs {
				max_abs = a;
				p = i;
			}
		}
		piv[k] = p;
		if p != k {
			for j in 0..n {
				re.swap(k * n + j, p * n + j);
				im.swap(k * n + j, p * n + j);
			}
		}
		let d_re = re[k * n + k];
		let d_im = im[k * n + k];
		if d_re.hypot(d_im) < PIVOT_FLOOR {
			re[k * n + k] = PIVOT_FLOOR;
			im[k * n + k] = 0.0;
		}
		let pk_re = re[k * n + k];
		let pk_im = im[k * n + k];
		for i in (k + 1)..n {
			let (mr, mi) = cdiv(re[i * n + k], im[i * n + k], pk_re, pk_im);
			re[i * n + k] = mr;
			im[i * n + k] = mi;
			for j in (k + 1)..n {
				let (pr, pi) = cmul(mr, mi, re[k * n + j], im[k * n + j]);
				re[i * n + j] -= pr;
				im[i * n + j] -= pi;
			}
		}
	}
}

fn lu_solve_cplx(
	lu_re: &[f64],
	lu_im: &[f64],
	n: usize,
	piv: &[usize],
	b_re: &mut [f64],
	b_im: &mut [f64],
) {
	for k in 0..n {
		let p = piv[k];
		if p != k {
			b_re.swap(k, p);
			b_im.swap(k, p);
		}
	}
	for i in 0..n {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for j in 0..i {
			let (pr, pi) = cmul(lu_re[i * n + j], lu_im[i * n + j], b_re[j], b_im[j]);
			sr -= pr;
			si -= pi;
		}
		b_re[i] = sr;
		b_im[i] = si;
	}
	for i in (0..n).rev() {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for j in (i + 1)..n {
			let (pr, pi) = cmul(lu_re[i * n + j], lu_im[i * n + j], b_re[j], b_im[j]);
			sr -= pr;
			si -= pi;
		}
		let (xr, xi) = cdiv(sr, si, lu_re[i * n + i], lu_im[i * n + i]);
		b_re[i] = xr;
		b_im[i] = xi;
	}
}

fn inverse_diag_cplx(re: &mut [f64], im: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
	let mut piv = vec![0usize; n];
	lu_factor_cplx(re, im, n, &mut piv);
	let mut y_re = vec![0.0f64; n];
	let mut y_im = vec![0.0f64; n];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		br.fill(0.0);
		bi.fill(0.0);
		br[j] = 1.0;
		lu_solve_cplx(re, im, n, &piv, &mut br, &mut bi);
		y_re[j] = br[j];
		y_im[j] = bi[j];
	}
	(y_re, y_im)
}

/// \(S = D^{-1/2}(R-D)\,\mathrm{solve}(R+D, D^{1/2})\),
/// \(T = 2\sqrt{Z_\mathrm{ref}}\,\mathrm{solve}(R+D, D^{1/2})\).
fn form_s_and_t(
	r_re: &[f64],
	r_im: &[f64],
	z0: &[f64],
	n: usize,
	z_ref: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut sqrtz = vec![0.0f64; n];
	for p in 0..n {
		sqrtz[p] = z0[p].max(EPS_Z).sqrt();
	}

	let mut l_re = r_re.to_vec();
	let mut l_im = r_im.to_vec();
	for p in 0..n {
		l_re[p * n + p] += z0[p];
		l_im[p * n + p] = 0.0;
	}
	chol_factor_herm(&mut l_re, &mut l_im, n);

	let mut x_re = vec![0.0f64; nn];
	let mut x_im = vec![0.0f64; nn];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		br.fill(0.0);
		bi.fill(0.0);
		br[j] = sqrtz[j];
		chol_solve_herm(&l_re, &l_im, n, &mut br, &mut bi);
		for i in 0..n {
			x_re[i * n + j] = br[i];
			x_im[i * n + j] = bi[i];
		}
	}

	let t_scale = 2.0 * z_ref.max(EPS_Z).sqrt();
	let mut t_re = vec![0.0f64; nn];
	let mut t_im = vec![0.0f64; nn];
	for i in 0..nn {
		t_re[i] = t_scale * x_re[i];
		t_im[i] = t_scale * x_im[i];
	}

	let mut rmd_re = r_re.to_vec();
	let rmd_im = r_im.to_vec();
	for p in 0..n {
		rmd_re[p * n + p] -= z0[p];
	}

	let mut s_re = vec![0.0f64; nn];
	let mut s_im = vec![0.0f64; nn];
	for i in 0..n {
		let inv_sqrt = 1.0 / sqrtz[i];
		for j in 0..n {
			let mut re = 0.0;
			let mut im = 0.0;
			for k in 0..n {
				let ar = rmd_re[i * n + k];
				let ai = rmd_im[i * n + k];
				let xr = x_re[k * n + j];
				let xi = x_im[k * n + j];
				re += ar * xr - ai * xi;
				im += ar * xi + ai * xr;
			}
			s_re[i * n + j] = re * inv_sqrt;
			s_im[i * n + j] = im * inv_sqrt;
		}
	}
	(s_re, s_im, t_re, t_im)
}

/// Kurokawa: \(S = G^{-1}(Z-D^*)(Z+D)^{-1}G\), \(T = 2\sqrt{Z_\mathrm{ref}}(Z+D)^{-1}G\),
/// \(G=\mathrm{diag}\sqrt{\Re(z_0)}\).
fn form_s_and_t_kurokawa(
	z_re: &[f64],
	z_im: &[f64],
	z0_re: &[f64],
	z0_im: &[f64],
	n: usize,
	z_ref: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut g = vec![0.0f64; n];
	for p in 0..n {
		g[p] = z0_re[p].max(EPS_Z).sqrt();
	}

	let mut lu_re = z_re.to_vec();
	let mut lu_im = z_im.to_vec();
	for p in 0..n {
		lu_re[p * n + p] += z0_re[p];
		lu_im[p * n + p] += z0_im[p];
	}
	let mut piv = vec![0usize; n];
	lu_factor_cplx(&mut lu_re, &mut lu_im, n, &mut piv);

	let mut w_re = vec![0.0f64; nn];
	let mut w_im = vec![0.0f64; nn];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		br.fill(0.0);
		bi.fill(0.0);
		br[j] = g[j];
		lu_solve_cplx(&lu_re, &lu_im, n, &piv, &mut br, &mut bi);
		for i in 0..n {
			w_re[i * n + j] = br[i];
			w_im[i * n + j] = bi[i];
		}
	}

	let t_scale = 2.0 * z_ref.max(EPS_Z).sqrt();
	let mut t_re = vec![0.0f64; nn];
	let mut t_im = vec![0.0f64; nn];
	for i in 0..nn {
		t_re[i] = t_scale * w_re[i];
		t_im[i] = t_scale * w_im[i];
	}

	// \(Z - D^*\): diagonal \(\,Z_{pp}-\overline{z_{0,p}}\).
	let mut zmd_re = z_re.to_vec();
	let mut zmd_im = z_im.to_vec();
	for p in 0..n {
		zmd_re[p * n + p] -= z0_re[p];
		zmd_im[p * n + p] += z0_im[p];
	}

	let mut s_re = vec![0.0f64; nn];
	let mut s_im = vec![0.0f64; nn];
	for i in 0..n {
		let inv_g = 1.0 / g[i];
		for j in 0..n {
			let mut re = 0.0;
			let mut im = 0.0;
			for k in 0..n {
				let ar = zmd_re[i * n + k];
				let ai = zmd_im[i * n + k];
				let wr = w_re[k * n + j];
				let wi = w_im[k * n + j];
				re += ar * wr - ai * wi;
				im += ar * wi + ai * wr;
			}
			s_re[i * n + j] = re * inv_g;
			s_im[i * n + j] = im * inv_g;
		}
	}
	(s_re, s_im, t_re, t_im)
}

#[cfg(test)]
mod tests {
	use super::*;

	fn close(a: f64, b: f64, tol: f64, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	fn expected_xpq(
		dx: f64,
		dy: f64,
		d_min: f64,
		x_nn: f64,
		alpha: f64,
		beta: f64,
		aniso: f64,
	) -> f64 {
		let rho = dx.hypot(dy);
		if rho <= 0.0 {
			return 0.0;
		}
		let envelope = x_nn * (d_min / rho).powf(alpha);
		let osc = (beta * (rho / d_min - 1.0)).cos();
		let phi = dy.atan2(dx);
		envelope * osc * (1.0 + aniso * (2.0 * phi).cos())
	}

	fn pair_d_min(x: &[f32], y: &[f32]) -> f64 {
		let n = x.len();
		let mut d_min = f64::INFINITY;
		for p in 0..n {
			for q in 0..p {
				let dx = x[p] as f64 - x[q] as f64;
				let dy = y[p] as f64 - y[q] as f64;
				let rho = dx.hypot(dy);
				if rho > 0.0 && rho < d_min {
					d_min = rho;
				}
			}
		}
		d_min
	}

	fn assert_x_kernel(
		m: &MatchedS,
		x: &[f32],
		y: &[f32],
		x_nn: f64,
		alpha: f64,
		beta: f64,
		aniso: f64,
	) {
		let n = x.len();
		let d_min = pair_d_min(x, y);
		for p in 0..n {
			close(m.r_im[p * n + p], 0.0, 0.0, "X_pp");
			for q in 0..p {
				let dx = x[p] as f64 - x[q] as f64;
				let dy = y[p] as f64 - y[q] as f64;
				let want = expected_xpq(dx, dy, d_min, x_nn, alpha, beta, aniso);
				close(m.r_im[p * n + q], want, 1e-9, "X_pq");
				close(m.r_im[q * n + p], want, 1e-9, "X_qp");
			}
		}
	}

	#[test]
	fn chol_solves_spd_identity() {
		let n = 3;
		#[rustfmt::skip]
		let a0 = [
			4.0, 1.0, 0.5,
			1.0, 3.0, 0.2,
			0.5, 0.2, 2.0,
		];
		let mut l = a0.to_vec();
		chol_factor_real(&mut l, n);
		for j in 0..n {
			let mut b = vec![0.0; n];
			b[j] = 1.0;
			chol_solve_real(&l, n, &mut b);
			let mut ax = vec![0.0; n];
			for i in 0..n {
				for k in 0..n {
					ax[i] += a0[i * n + k] * b[k];
				}
			}
			for i in 0..n {
				let expect = if i == j { 1.0 } else { 0.0 };
				close(ax[i], expect, 1e-12, "A x = e_j");
			}
		}
	}

	#[test]
	fn lu_solves_complex() {
		let n = 2;
		#[rustfmt::skip]
		let a_re = [
			1.0, 2.0,
			3.0, 4.0,
		];
		#[rustfmt::skip]
		let a_im = [
			1.0, 0.0,
			0.0, -1.0,
		];
		let mut lu_re = a_re.to_vec();
		let mut lu_im = a_im.to_vec();
		let mut piv = vec![0usize; n];
		lu_factor_cplx(&mut lu_re, &mut lu_im, n, &mut piv);
		let mut br = vec![1.0, 0.0];
		let mut bi = vec![0.0, 0.0];
		lu_solve_cplx(&lu_re, &lu_im, n, &piv, &mut br, &mut bi);
		let mut ax_re = [0.0; 2];
		let mut ax_im = [0.0; 2];
		for i in 0..n {
			for k in 0..n {
				let (pr, pi) = cmul(a_re[i * n + k], a_im[i * n + k], br[k], bi[k]);
				ax_re[i] += pr;
				ax_im[i] += pi;
			}
		}
		close(ax_re[0], 1.0, 1e-12, "lu Ax re0");
		close(ax_re[1], 0.0, 1e-12, "lu Ax re1");
		close(ax_im[0], 0.0, 1e-12, "lu Ax im0");
		close(ax_im[1], 0.0, 1e-12, "lu Ax im1");
	}

	#[test]
	fn n1_gram_gives_zref_and_s0() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let m = MatchedS::from_gram(&p_re, &p_im, 1, Z_REF);
		assert_eq!(m.n, 1);
		close(m.z_ref, Z_REF, 0.0, "z_ref");
		close(m.r_re[0], Z_REF, 1e-12, "R11");
		close(m.r_im[0], 0.0, 0.0, "R11 im");
		close(m.z0[0], Z_REF, 1e-9, "z0");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		close(m.s_re[0], 0.0, 1e-12, "S11 re");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		close(m.t_im[0], 0.0, 1e-12, "T11 im");
		assert!(m.residual < TAU, "residual {}", m.residual);
	}

	#[test]
	fn two_port_real_r_matches_and_s_symmetric() {
		// Isolated P_H so R = [[50, 12.5], [12.5, 50]].
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let m = MatchedS::from_gram(&p_re, &p_im, 2, Z_REF);
		close(m.r_re[0], 50.0, 1e-12, "R11");
		close(m.r_re[1], 12.5, 1e-12, "R12");
		assert!(m.residual < TAU, "residual {}", m.residual);
		close(m.z0[0], m.z0[1], 1e-9, "equal ports");
		close(m.s_re[0], 0.0, 2e-4, "S11 ~ 0");
		close(m.s_re[3], 0.0, 2e-4, "S22 ~ 0");
		close(m.s_re[1], m.s_re[2], 1e-12, "S12 = S21");
		close(m.s_im[1], 0.0, 1e-12, "S12 im");
		assert!(m.s_re[1].abs() > 0.05, "coupled |S12|={}", m.s_re[1]);
	}

	#[test]
	fn uncoupled_two_port_is_matched_open() {
		let p_re = [0.5f32, 0.0, 0.0, 0.5];
		let p_im = [0.0f32; 4];
		let m = MatchedS::from_gram(&p_re, &p_im, 2, Z_REF);
		close(m.z0[0], Z_REF, 1e-9, "z0_1");
		close(m.z0[1], Z_REF, 1e-9, "z0_2");
		for v in m.s_re.iter().chain(m.s_im.iter()) {
			close(*v, 0.0, 1e-12, "S");
		}
		close(m.t_re[0], 1.0, 1e-12, "T11");
		close(m.t_re[3], 1.0, 1e-12, "T22");
		close(m.t_re[1], 0.0, 1e-12, "T12");
		close(m.t_re[2], 0.0, 1e-12, "T21");
	}

	fn gemv(m_re: &[f64], m_im: &[f64], a_re: &[f64], a_im: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
		let mut w_re = vec![0.0; n];
		let mut w_im = vec![0.0; n];
		for i in 0..n {
			let mut re = 0.0;
			let mut im = 0.0;
			for j in 0..n {
				let mr = m_re[i * n + j];
				let mi = m_im[i * n + j];
				re += mr * a_re[j] - mi * a_im[j];
				im += mr * a_im[j] + mi * a_re[j];
			}
			w_re[i] = re;
			w_im[i] = im;
		}
		(w_re, w_im)
	}

	#[test]
	fn two_port_t_maps_incident_to_weights() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let m = MatchedS::from_gram(&p_re, &p_im, 2, Z_REF);
		close(m.t_im.iter().fold(0.0f64, |a, v| a.max(v.abs())), 0.0, 1e-12, "T im");
		assert!((m.t_re[0] - 1.0).abs() > 1e-4, "coupled T11 != 1");
		let a_re = [1.0, 0.0];
		let a_im = [0.0, 0.0];
		let (w_re, w_im) = gemv(&m.t_re, &m.t_im, &a_re, &a_im, 2);
		close(w_re[0], m.t_re[0], 1e-12, "T col0");
		close(w_re[1], m.t_re[2], 1e-12, "T col0 row1");
		close(w_im[0], 0.0, 1e-12, "w im");
		close(w_im[1], 0.0, 1e-12, "w im");
	}

	#[test]
	fn herm_chol_agrees_with_real_on_real_matrix() {
		let n = 2;
		#[rustfmt::skip]
		let a = [
			4.0, 1.0,
			1.0, 3.0,
		];
		let mut l_re = a.to_vec();
		let mut l_im = vec![0.0; 4];
		chol_factor_herm(&mut l_re, &mut l_im, n);
		let mut br = vec![1.0, 0.0];
		let mut bi = vec![0.0, 0.0];
		chol_solve_herm(&l_re, &l_im, n, &mut br, &mut bi);
		let mut ax = [0.0; 2];
		for i in 0..2 {
			for k in 0..2 {
				ax[i] += a[i * 2 + k] * br[k];
			}
		}
		close(ax[0], 1.0, 1e-12, "herm Ax0");
		close(ax[1], 0.0, 1e-12, "herm Ax1");
		close(bi[0], 0.0, 1e-12, "herm im0");
		close(bi[1], 0.0, 1e-12, "herm im1");
	}

	#[test]
	fn xnn_zero_matches_from_gram() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let a = MatchedS::from_gram(&p_re, &p_im, 2, Z_REF);
		let b = MatchedS::from_gram_coupled(
			&p_re,
			&p_im,
			2,
			Z_REF,
			&[0.0, 0.5],
			&[0.0, 0.0],
			0.0,
			2.0,
			0.0,
			0.0,
			0.0,
			0.0,
		);
		close(a.z0[0], b.z0[0], 0.0, "z0");
		close(a.s_re[1], b.s_re[1], 0.0, "S12");
		close(b.z0_im[0], 0.0, 0.0, "z0 im");
		close(b.r_im[1], 0.0, 0.0, "X12");
	}

	#[test]
	fn two_port_reactance_is_symmetric_and_matched() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let x_nn = 10.0;
		let m = MatchedS::from_gram_coupled(
			&p_re,
			&p_im,
			2,
			Z_REF,
			&[0.0, 0.5],
			&[0.0, 0.0],
			x_nn,
			2.0,
			0.0,
			0.0,
			0.0,
			0.0,
		);
		close(m.r_re[0], 50.0, 1e-12, "R11");
		close(m.r_im[1], x_nn, 1e-12, "X12");
		close(m.r_im[2], x_nn, 1e-12, "X21");
		close(m.r_im[0], 0.0, 0.0, "X11");
		assert!(m.residual < TAU, "residual {}", m.residual);
		let s11 = m.s_re[0].hypot(m.s_im[0]);
		let s22 = m.s_re[3].hypot(m.s_im[3]);
		assert!(s11 < 2e-3, "|S11|={s11}");
		assert!(s22 < 2e-3, "|S22|={s22}");
		assert!(m.z0_im[0].abs() > 1e-6, "z0 imag {}", m.z0_im[0]);
		close(m.z0_im[0], m.z0_im[1], 1e-9, "equal Im z0");
	}

	#[test]
	fn three_irregular_decay_and_alpha_zero() {
		let n = 3;
		let mut p_re = vec![0.5f32; n * n];
		for p in 0..n {
			p_re[p * n + p] = 0.5;
		}
		// Off-diag 0.05 → R_pq = 5 Ω.
		for p in 0..n {
			for q in 0..n {
				if p != q {
					p_re[p * n + q] = 0.05;
				}
			}
		}
		let p_im = vec![0.0f32; n * n];
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.25, -0.25];
		let x_nn = 8.0;
		let alpha = 2.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, alpha, 0.0, 0.0, 0.0, 0.0,
		);
		assert!(m.residual < TAU, "residual {}", m.residual);

		fn rho(i: usize, j: usize, x: &[f32], y: &[f32]) -> f64 {
			let dx = x[i] as f64 - x[j] as f64;
			let dy = y[i] as f64 - y[j] as f64;
			dx.hypot(dy)
		}
		let r01 = rho(0, 1, &x, &y);
		let r02 = rho(0, 2, &x, &y);
		let r12 = rho(1, 2, &x, &y);
		let d_min = r01.min(r02).min(r12);
		let expect = |rij: f64| x_nn * (d_min / rij).powf(alpha);
		close(m.r_im[0 * n + 1], expect(r01), 1e-9, "X01");
		close(m.r_im[1 * n + 0], expect(r01), 1e-9, "X10");
		close(m.r_im[0 * n + 2], expect(r02), 1e-9, "X02");
		close(m.r_im[1 * n + 2], expect(r12), 1e-9, "X12");
		close(m.r_im[0], 0.0, 0.0, "X00");
		let mag00 = m.s_re[0].hypot(m.s_im[0]);
		assert!(mag00 < 2e-3, "|S00|={mag00}");
		assert!(
			expect(r01).abs() >= expect(r02).abs() - 1e-12,
			"closest pair is strongest"
		);

		let m0 = MatchedS::from_gram_coupled(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, 0.0, 0.0, 0.0, 0.0, 0.0,
		);
		close(m0.r_im[1], x_nn, 1e-12, "alpha0 X01");
		close(m0.r_im[2], x_nn, 1e-12, "alpha0 X02");
		close(m0.r_im[5], x_nn, 1e-12, "alpha0 X12");
	}

	#[test]
	fn coincident_pair_is_finite() {
		let p_re = [0.5f32, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
		let p_im = [0.0f32; 9];
		let m = MatchedS::from_gram_coupled(
			&p_re,
			&p_im,
			3,
			Z_REF,
			&[0.0, 0.0, 1.0],
			&[0.0, 0.0, 0.0],
			12.0,
			3.0,
			0.0,
			0.0,
			0.0,
			0.0,
		);
		assert!(m.r_im.iter().all(|v| v.is_finite()), "X finite");
		close(m.r_im[1], 0.0, 0.0, "coincident X01");
		close(m.r_im[2], 12.0, 1e-12, "X02 = Xnn");
		close(m.r_im[5], 12.0, 1e-12, "X12 = Xnn");
		assert!(m.residual.is_finite());
	}

	#[test]
	fn n1_common_zref_is_open() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 1, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, Z_REF, 0.0,
		);
		close(m.z0[0], Z_REF, 0.0, "z0");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		close(m.s_re[0], 0.0, 1e-12, "S11 re");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		assert_eq!(m.iterations, 0);
	}

	#[test]
	fn n1_common_real_mismatch() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let zc = 40.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 1, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, zc, 0.0,
		);
		let r = 50.0;
		let s = (r - zc) / (r + zc);
		close(m.s_re[0], s, 1e-12, "S11");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.z0[0], zc, 0.0, "z0");
	}

	#[test]
	fn n1_common_complex_kurokawa() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let zc_re = 50.0;
		let zc_im = 10.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 1, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, zc_re, zc_im,
		);
		// S = (Z - zc*)/(Z + zc) with Z = 50.
		let num_re = 50.0 - zc_re;
		let num_im = zc_im;
		let den_re = 50.0 + zc_re;
		let den_im = zc_im;
		let d2 = den_re * den_re + den_im * den_im;
		let s_re = (num_re * den_re + num_im * den_im) / d2;
		let s_im = (num_im * den_re - num_re * den_im) / d2;
		close(m.s_re[0], s_re, 1e-12, "S11 re");
		close(m.s_im[0], s_im, 1e-12, "S11 im");
		close(m.z0[0], zc_re, 0.0, "z0 re");
		close(m.z0_im[0], zc_im, 0.0, "z0 im");
	}

	#[test]
	fn two_port_common_z0_equal_and_sii_not_zero() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let per = MatchedS::from_gram(&p_re, &p_im, 2, Z_REF);
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 2, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, Z_REF, 0.0,
		);
		close(m.z0[0], Z_REF, 0.0, "z0_0");
		close(m.z0[1], Z_REF, 0.0, "z0_1");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		let sii = m.s_re[0].hypot(m.s_im[0]).max(m.s_re[3].hypot(m.s_im[3]));
		let per_sii = per.s_re[0].hypot(per.s_im[0]);
		assert!(sii > 0.01, "common |Sii|={sii}");
		assert!(sii > per_sii * 10.0, "common |Sii| vs per-port {sii} vs {per_sii}");

		// Closed form S = (Z - zc I)(Z + zc I)^{-1} at zc = 40.
		let zc = 40.0;
		let m40 = MatchedS::from_gram_coupled(
			&p_re, &p_im, 2, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, zc, 0.0,
		);
		let det = 90.0 * 90.0 - 12.5 * 12.5;
		let s11 = (10.0 * 90.0 + 12.5 * (-12.5)) / det;
		let s12 = (10.0 * (-12.5) + 12.5 * 90.0) / det;
		close(m40.s_re[0], s11, 1e-12, "S11");
		close(m40.s_re[1], s12, 1e-12, "S12");
		close(m40.s_re[2], s12, 1e-12, "S21");
		close(m40.s_re[3], s11, 1e-12, "S22");
		close(m40.s_im.iter().fold(0.0f64, |a, v| a.max(v.abs())), 0.0, 1e-12, "S im");
	}

	#[test]
	fn collinear_beta_pi_flips_next_nearest() {
		let n = 3;
		let p_re = [0.5f32, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
		let p_im = [0.0f32; 9];
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.0, 0.0];
		let x_nn = 10.0;
		let alpha = 2.0;
		let beta = std::f64::consts::PI;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, alpha, beta, 0.0, 0.0, 0.0,
		);
		assert_x_kernel(&m, &x, &y, x_nn, alpha, beta, 0.0);
		close(m.r_im[1], x_nn, 1e-12, "X01 nn");
		close(m.r_im[5], x_nn, 1e-12, "X12 nn");
		close(m.r_im[2], -x_nn * 0.25, 1e-12, "X02 opposite sign");
		assert!(m.r_im[2] * m.r_im[1] < 0.0, "next-nearest opposite sign");
	}

	#[test]
	fn right_angle_aniso_splits_equal_distance() {
		let n = 3;
		let p_re = [0.5f32, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
		let p_im = [0.0f32; 9];
		let x = [0.0f32, 1.0, 0.0];
		let y = [0.0f32, 0.0, 1.0];
		let x_nn = 10.0;
		let aniso = 0.5;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, 2.0, 0.0, aniso, 0.0, 0.0,
		);
		assert_x_kernel(&m, &x, &y, x_nn, 2.0, 0.0, aniso);
		close(m.r_im[1], x_nn * (1.0 + aniso), 1e-12, "X along +x");
		close(m.r_im[2], x_nn * (1.0 - aniso), 1e-12, "X along +y");
		assert!((m.r_im[1].abs() - m.r_im[2].abs()).abs() > 1.0);
		close(m.r_im[1], m.r_im[3], 1e-12, "X01 = X10");
		close(m.r_im[2], m.r_im[6], 1e-12, "X02 = X20");
	}

	#[test]
	fn sunflower_like_cloud_matches_pair_formula() {
		let n_vogel = 7;
		let n = n_vogel + 3;
		let gd = (5.0f64.sqrt() - 1.0) / 2.0;
		let mut x = vec![0.0f32; n];
		let mut y = vec![0.0f32; n];
		for i in 0..n_vogel {
			let t = 2.0 * std::f64::consts::PI * gd * (i + 1) as f64;
			let r = t.sqrt() * 0.22;
			x[i] = (r * t.cos()) as f32;
			y[i] = (r * t.sin()) as f32;
		}
		x[n_vogel] = 0.0;
		y[n_vogel] = 0.0;
		x[n_vogel + 1] = 2.0;
		y[n_vogel + 1] = 0.0;
		x[n_vogel + 2] = 0.0;
		y[n_vogel + 2] = 2.0;
		let mut p_re = vec![0.0f32; n * n];
		for p in 0..n {
			p_re[p * n + p] = 0.5;
		}
		let p_im = vec![0.0f32; n * n];
		let x_nn = 8.0;
		let alpha = 2.0;
		let beta = std::f64::consts::FRAC_PI_2;
		let aniso = 0.4;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, alpha, beta, aniso, 0.0, 0.0,
		);
		assert!(m.r_im.iter().all(|v| v.is_finite()), "X finite");
		assert_x_kernel(&m, &x, &y, x_nn, alpha, beta, aniso);
		let o = n_vogel;
		assert!(
			(m.r_im[o * n + o + 1] - m.r_im[o * n + o + 2]).abs() > 1e-9,
			"equal-ρ x/y arms split"
		);
	}
}
