//! Radiative resistance, simultaneous port match, and power-wave \(S\).
//!
//! §6–8 of `docs/approximate_matched_basis.md`: \(Z=R+jX_\mathrm{mutual}+jX_\mathrm{self}I\),
//! a real per-port or common \(z_0\), and real-reference \(S\) and \(T\). Operates
//! Internals are f64; the Gram itself stays f32. Dense solves go through
//! [`crate::linalg`] (faer). Textbook LU/Cholesky remains in
//! [`crate::legacy_linalg`] for parity tests only.

pub const Z_REF: f64 = 50.0;
pub const EPS_Z: f64 = 1e-9;
pub const K_MAX: u32 = 200;
pub const TAU: f64 = 1e-3;

const MATCH_BETA: f64 = 0.5;

use crate::linalg as la;

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

	/// Same as [`from_gram`], then \(Z = R + jX(\Delta x,\Delta y) + j X_\mathrm{self} I\).
	/// Per-port real match on \(\Re(Z_\mathrm{in})\) unless \(\Re(z_c)>0\), in which
	/// case every port uses that common real \(z_c\).
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
		x_self: f64,
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

		let mutual = add_mutual_reactance(&mut r_im, x, y, n, x_nn, alpha, beta, aniso);
		let x_self = if x_self.is_finite() { x_self } else { 0.0 };
		if x_self != 0.0 {
			for p in 0..n {
				r_im[p * n + p] = x_self;
			}
		}
		let has_x = mutual || x_self != 0.0;
		if z_common_re.is_finite() && z_common_re > 0.0 {
			Self::from_z_common(r_re, r_im, n, z_ref, z_common_re, has_x)
		} else if has_x {
			Self::from_z_real_reactive(r_re, r_im, n, z_ref)
		} else {
			Self::from_z_real(r_re, r_im, n, z_ref)
		}
	}

	/// Same Gram \(R\) as [`from_gram_coupled`], then the propagation overlay
	/// \(X(\lvert\Delta x\rvert,\lvert\Delta y\rvert;\varepsilon_x,\varepsilon_y,\alpha_\lambda,f)\).
	/// Always a common real \(z_c\) (non-positive or non-finite \(z_c\) becomes \(z_\mathrm{ref}\)).
	pub fn from_gram_propagation(
		p_re: &[f32],
		p_im: &[f32],
		n: usize,
		z_ref: f64,
		x: &[f32],
		y: &[f32],
		x_nn: f64,
		att: f64,
		eps_x: f64,
		eps_y: f64,
		freq: f64,
		z_common_re: f64,
		x_self: f64,
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

		let mutual = add_propagation_reactance(&mut r_im, x, y, n, x_nn, att, eps_x, eps_y, freq);
		let x_self = if x_self.is_finite() { x_self } else { 0.0 };
		if x_self != 0.0 {
			for p in 0..n {
				r_im[p * n + p] = x_self;
			}
		}
		let has_x = mutual || x_self != 0.0;
		let zc = if z_common_re.is_finite() && z_common_re > 0.0 {
			z_common_re
		} else {
			z_ref
		};
		Self::from_z_common(r_re, r_im, n, z_ref, zc, has_x)
	}

	/// Complex \(Z\) already in ohms. Does **not** scale a Gram. Always a
	/// common real \(z_c\) (non-positive or non-finite \(z_c\) becomes
	/// \(z_\mathrm{ref}\)). Finite \(x_\mathrm{self}\) is **added** to every
	/// diagonal of \(\Im(Z)\) (series reactance). Green-mode match entry.
	pub fn from_z(
		z_re: &[f64],
		z_im: &[f64],
		n: usize,
		z_ref: f64,
		z_common_re: f64,
		x_self: f64,
	) -> Self {
		let z_ref = if z_ref.is_finite() && z_ref > 0.0 {
			z_ref
		} else {
			Z_REF
		};
		if n == 0 || z_re.len() != n * n {
			return Self::empty();
		}
		let nn = n * n;
		let r_re = z_re.to_vec();
		let mut r_im = vec![0.0f64; nn];
		if z_im.len() == nn {
			r_im.copy_from_slice(z_im);
		}
		let x_self = if x_self.is_finite() { x_self } else { 0.0 };
		if x_self != 0.0 {
			for p in 0..n {
				r_im[p * n + p] += x_self;
			}
		}
		let has_x = r_im.iter().any(|&x| x != 0.0);
		let zc = if z_common_re.is_finite() && z_common_re > 0.0 {
			z_common_re
		} else {
			z_ref
		};
		Self::from_z_common(r_re, r_im, n, z_ref, zc, has_x)
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
			let ydiag = la::inverse_diag_spd(&work, n);
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

		let st = form_s_and_t(&r_re, &r_im, &z0, n, z_ref);
		let residual = residual_real_port(&st.ydiag_re, &z0);
		Self {
			n,
			z_ref,
			z0,
			z0_im: vec![0.0f64; n],
			r_re,
			r_im,
			s_re: st.s_re,
			s_im: st.s_im,
			t_re: st.t_re,
			t_im: st.t_im,
			iterations,
			residual,
		}
	}

	/// Real \(z_0\) on \(\Re(Z_\mathrm{in})\) of the full complex \(Z\); LU for \(S\), \(T\).
	fn from_z_real_reactive(r_re: Vec<f64>, r_im: Vec<f64>, n: usize, z_ref: f64) -> Self {
		let nn = n * n;
		let mut z0 = vec![0.0f64; n];
		for p in 0..n {
			z0[p] = r_re[p * n + p].max(EPS_Z);
		}
		let z0_im = vec![0.0f64; n];

		let mut iterations = 0u32;
		let mut work_re = vec![0.0f64; nn];
		let mut work_im = vec![0.0f64; nn];
		let mut beta = 1.0f64;
		for k in 0..K_MAX {
			iterations += 1;
			if k == K_MAX / 4 {
				beta = MATCH_BETA;
			}
			fill_z_plus_d(&r_re, &r_im, &z0, &z0_im, n, &mut work_re, &mut work_im);
			let (ydiag_re, ydiag_im) = la::inverse_diag_cplx(&work_re, &work_im, n);
			let mut z_next = vec![0.0f64; n];
			let mut delta = 0.0f64;
			for p in 0..n {
				let (zin_re, _zin_im) = zin_from_ypp(ydiag_re[p], ydiag_im[p], z0[p], 0.0);
				let star = zin_re.max(EPS_Z);
				z_next[p] = (1.0 - beta) * z0[p] + beta * star;
				if z_next[p] < EPS_Z {
					z_next[p] = EPS_Z;
				}
				delta = delta.max((star - z0[p]).abs());
			}
			z0 = z_next;
			if delta < TAU {
				break;
			}
		}

		let st = form_s_and_t_kurokawa(&r_re, &r_im, &z0, &z0_im, n, z_ref);
		let residual = residual_real_port_cplx(&st.ydiag_re, &st.ydiag_im, &z0);
		Self {
			n,
			z_ref,
			z0,
			z0_im,
			r_re,
			r_im,
			s_re: st.s_re,
			s_im: st.s_im,
			t_re: st.t_re,
			t_im: st.t_im,
			iterations,
			residual,
		}
	}

	/// Skip the solver: every port uses the same real \(z_c\) (\(\Re(z_c)>0\)).
	fn from_z_common(
		r_re: Vec<f64>,
		r_im: Vec<f64>,
		n: usize,
		z_ref: f64,
		zc_re: f64,
		has_x: bool,
	) -> Self {
		let zc_re = zc_re.max(EPS_Z);
		let z0 = vec![zc_re; n];
		let z0_im = vec![0.0f64; n];

		let st = if has_x {
			form_s_and_t_kurokawa(&r_re, &r_im, &z0, &z0_im, n, z_ref)
		} else {
			form_s_and_t(&r_re, &r_im, &z0, n, z_ref)
		};
		let residual = if has_x {
			residual_common_cplx(&st.ydiag_re, &st.ydiag_im, &z0, zc_re)
		} else {
			residual_common_real(&st.ydiag_re, zc_re)
		};

		Self {
			n,
			z_ref,
			z0,
			z0_im,
			r_re,
			r_im,
			s_re: st.s_re,
			s_im: st.s_im,
			t_re: st.t_re,
			t_im: st.t_im,
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

/// \(X_{pq}=X_{nn}(\mathrm{env}/\mathrm{env}_\mathrm{ref})\cos(\phi-\phi_\mathrm{ref})\),
/// \(\phi=2\pi f(\sqrt{\varepsilon_x}\lvert\Delta x\rvert+\sqrt{\varepsilon_y}\lvert\Delta y\rvert)\),
/// \(\mathrm{env}=e^{-\alpha_\lambda\rho}\). Ref is the closest pair. Returns true if any
/// off-diagonal was written.
fn add_propagation_reactance(
	r_im: &mut [f64],
	x: &[f32],
	y: &[f32],
	n: usize,
	x_nn: f64,
	att: f64,
	eps_x: f64,
	eps_y: f64,
	freq: f64,
) -> bool {
	if n < 2 || x.len() != n || y.len() != n {
		return false;
	}
	if !x_nn.is_finite() || x_nn == 0.0 {
		return false;
	}
	let att = if att.is_finite() && att >= 0.0 { att } else { 0.0 };
	let eps_x = if eps_x.is_finite() && eps_x >= 0.0 {
		eps_x
	} else {
		1.0
	};
	let eps_y = if eps_y.is_finite() && eps_y >= 0.0 {
		eps_y
	} else {
		1.0
	};
	let freq = if freq.is_finite() && freq > 0.0 { freq } else { 1.0 };
	let kx = 2.0 * std::f64::consts::PI * freq * eps_x.sqrt();
	let ky = 2.0 * std::f64::consts::PI * freq * eps_y.sqrt();

	let mut d_min = f64::INFINITY;
	let mut dx_ref = 0.0;
	let mut dy_ref = 0.0;
	for p in 0..n {
		for q in 0..p {
			let dx = (x[p] as f64 - x[q] as f64).abs();
			let dy = (y[p] as f64 - y[q] as f64).abs();
			let rho = dx.hypot(dy);
			if rho > 0.0 && rho < d_min {
				d_min = rho;
				dx_ref = dx;
				dy_ref = dy;
			}
		}
	}
	if !d_min.is_finite() || d_min <= 0.0 {
		return false;
	}
	let phi_ref = kx * dx_ref + ky * dy_ref;
	let env_ref = (-att * d_min).exp();
	if env_ref == 0.0 || !env_ref.is_finite() {
		return false;
	}

	let mut wrote = false;
	for p in 0..n {
		for q in 0..p {
			let dx = (x[p] as f64 - x[q] as f64).abs();
			let dy = (y[p] as f64 - y[q] as f64).abs();
			let rho = dx.hypot(dy);
			if !(rho > 0.0) {
				continue;
			}
			let env = (-att * rho).exp();
			let phi = kx * dx + ky * dy;
			let xpq = x_nn * (env / env_ref) * (phi - phi_ref).cos();
			if !xpq.is_finite() || xpq == 0.0 {
				continue;
			}
			r_im[p * n + q] += xpq;
			r_im[q * n + p] += xpq;
			wrote = true;
		}
	}
	wrote
}

struct StResult {
	s_re: Vec<f64>,
	s_im: Vec<f64>,
	t_re: Vec<f64>,
	t_im: Vec<f64>,
	ydiag_re: Vec<f64>,
	ydiag_im: Vec<f64>,
}

fn zin_from_ypp_real(ypp: f64, z0: f64) -> f64 {
	if ypp.abs() < 1e-30 {
		f64::INFINITY
	} else {
		1.0 / ypp - z0
	}
}

/// Per-port real match residual: \(\max_p\lvert \operatorname{Re}(Z_\mathrm{in})-z_{0,p}\rvert\).
fn residual_real_port(ydiag: &[f64], z0: &[f64]) -> f64 {
	let n = z0.len();
	let mut residual = 0.0f64;
	for p in 0..n {
		let zin = zin_from_ypp_real(ydiag[p], z0[p]);
		residual = residual.max((zin.max(EPS_Z) - z0[p]).abs());
	}
	residual
}

/// Per-port real match residual from complex \(Y_{pp}\).
fn residual_real_port_cplx(ydiag_re: &[f64], ydiag_im: &[f64], z0: &[f64]) -> f64 {
	let n = z0.len();
	let mut residual = 0.0f64;
	for p in 0..n {
		let (zin_re, _zin_im) = zin_from_ypp(ydiag_re[p], ydiag_im[p], z0[p], 0.0);
		residual = residual.max((zin_re.max(EPS_Z) - z0[p]).abs());
	}
	residual
}

/// Common \(z_c\) residual for real \(Z\): \(\max_p\lvert Z_\mathrm{in}-z_c\rvert\).
fn residual_common_real(ydiag: &[f64], zc: f64) -> f64 {
	let mut residual = 0.0f64;
	for ypp in ydiag {
		residual = residual.max((zin_from_ypp_real(*ypp, zc) - zc).abs());
	}
	residual
}

/// Common \(z_c\) residual for complex \(Z\): \(\max_p\lvert Z_\mathrm{in}-z_c\rvert\).
fn residual_common_cplx(ydiag_re: &[f64], ydiag_im: &[f64], z0: &[f64], zc: f64) -> f64 {
	let n = z0.len();
	let mut residual = 0.0f64;
	for p in 0..n {
		let (zin_re, zin_im) = zin_from_ypp(ydiag_re[p], ydiag_im[p], z0[p], 0.0);
		residual = residual.max((zin_re - zc).hypot(zin_im));
	}
	residual
}

fn zin_from_ypp(ypp_re: f64, ypp_im: f64, z0_re: f64, z0_im: f64) -> (f64, f64) {
	if ypp_re.hypot(ypp_im) < 1e-30 {
		return (f64::INFINITY, 0.0);
	}
	let d = ypp_re * ypp_re + ypp_im * ypp_im;
	let inv_re = ypp_re / d;
	let inv_im = -ypp_im / d;
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

/// \(S = D^{-1/2}(R-D)\,\mathrm{solve}(R+D, D^{1/2})\),
/// \(T = 2\sqrt{Z_\mathrm{ref}}\,\mathrm{solve}(R+D, D^{1/2})\).
fn form_s_and_t(
	r_re: &[f64],
	r_im: &[f64],
	z0: &[f64],
	n: usize,
	z_ref: f64,
) -> StResult {
	form_s_and_t_with(r_re, r_im, z0, n, z_ref, la::solve_herm_multi)
}

fn form_s_and_t_with(
	r_re: &[f64],
	r_im: &[f64],
	z0: &[f64],
	n: usize,
	z_ref: f64,
	solve: fn(&[f64], &[f64], usize, &[f64], &[f64]) -> (Vec<f64>, Vec<f64>),
) -> StResult {
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
	let mut b_re = vec![0.0f64; nn];
	let b_im = vec![0.0f64; nn];
	for j in 0..n {
		b_re[j * n + j] = sqrtz[j];
	}
	let (x_re, x_im) = solve(&l_re, &l_im, n, &b_re, &b_im);

	let mut ydiag_re = vec![0.0f64; n];
	for j in 0..n {
		ydiag_re[j] = x_re[j * n + j] / sqrtz[j];
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

	let (prod_re, prod_im) = la::matmul_cplx(&rmd_re, &rmd_im, &x_re, &x_im, n);
	let mut s_re = vec![0.0f64; nn];
	let mut s_im = vec![0.0f64; nn];
	for i in 0..n {
		let inv_sqrt = 1.0 / sqrtz[i];
		for j in 0..n {
			s_re[i * n + j] = prod_re[i * n + j] * inv_sqrt;
			s_im[i * n + j] = prod_im[i * n + j] * inv_sqrt;
		}
	}
	StResult {
		s_re,
		s_im,
		t_re,
		t_im,
		ydiag_re,
		ydiag_im: vec![0.0f64; n],
	}
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
) -> StResult {
	form_s_and_t_kurokawa_with(z_re, z_im, z0_re, z0_im, n, z_ref, la::solve_cplx_multi)
}

fn form_s_and_t_kurokawa_with(
	z_re: &[f64],
	z_im: &[f64],
	z0_re: &[f64],
	z0_im: &[f64],
	n: usize,
	z_ref: f64,
	solve: fn(&[f64], &[f64], usize, &[f64], &[f64]) -> (Vec<f64>, Vec<f64>),
) -> StResult {
	let nn = n * n;
	let mut g = vec![0.0f64; n];
	for p in 0..n {
		g[p] = z0_re[p].max(EPS_Z).sqrt();
	}

	let mut a_re = z_re.to_vec();
	let mut a_im = z_im.to_vec();
	for p in 0..n {
		a_re[p * n + p] += z0_re[p];
		a_im[p * n + p] += z0_im[p];
	}
	let mut b_re = vec![0.0f64; nn];
	let b_im = vec![0.0f64; nn];
	for j in 0..n {
		b_re[j * n + j] = g[j];
	}
	let (w_re, w_im) = solve(&a_re, &a_im, n, &b_re, &b_im);

	let mut ydiag_re = vec![0.0f64; n];
	let mut ydiag_im = vec![0.0f64; n];
	for j in 0..n {
		ydiag_re[j] = w_re[j * n + j] / g[j];
		ydiag_im[j] = w_im[j * n + j] / g[j];
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

	let (prod_re, prod_im) = la::matmul_cplx(&zmd_re, &zmd_im, &w_re, &w_im, n);
	let mut s_re = vec![0.0f64; nn];
	let mut s_im = vec![0.0f64; nn];
	for i in 0..n {
		let inv_g = 1.0 / g[i];
		for j in 0..n {
			s_re[i * n + j] = prod_re[i * n + j] * inv_g;
			s_im[i * n + j] = prod_im[i * n + j] * inv_g;
		}
	}
	StResult {
		s_re,
		s_im,
		t_re,
		t_im,
		ydiag_re,
		ydiag_im,
	}
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
		x_self: f64,
	) {
		let n = x.len();
		let d_min = pair_d_min(x, y);
		for p in 0..n {
			close(m.r_im[p * n + p], x_self, 0.0, "X_pp");
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
	fn from_z_n1_real_is_open() {
		let z_re = [Z_REF];
		let z_im = [0.0];
		let m = MatchedS::from_z(&z_re, &z_im, 1, Z_REF, Z_REF, 0.0);
		assert_eq!(m.n, 1);
		close(m.z0[0], Z_REF, 0.0, "z0");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		close(m.r_re[0], Z_REF, 0.0, "Z11");
		close(m.s_re[0], 0.0, 1e-12, "S11 re");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		close(m.t_im[0], 0.0, 1e-12, "T11 im");
		assert_eq!(m.iterations, 0);
	}

	#[test]
	fn from_z_adds_xself_on_diag_only() {
		let z_re = [10.0, 1.0, 1.0, 10.0];
		let z_im = [3.0, 0.5, 0.5, 3.0];
		let x_self = 7.0;
		let m = MatchedS::from_z(&z_re, &z_im, 2, Z_REF, 10.0, x_self);
		close(m.r_im[0], 10.0, 0.0, "X11");
		close(m.r_im[3], 10.0, 0.0, "X22");
		close(m.r_im[1], 0.5, 0.0, "X12");
		close(m.r_im[2], 0.5, 0.0, "X21");
		close(m.r_re[0], 10.0, 0.0, "R11");
		close(m.r_re[1], 1.0, 0.0, "R12");
	}

	#[test]
	fn from_z_n1_cancelled_x_is_open() {
		let z_re = [Z_REF];
		let z_im = [-12.0];
		let m = MatchedS::from_z(&z_re, &z_im, 1, Z_REF, Z_REF, 12.0);
		close(m.r_im[0], 0.0, 0.0, "X11");
		close(m.z0[0], Z_REF, 0.0, "z0");
		close(m.s_re[0], 0.0, 1e-12, "S11 re");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		close(m.t_im[0], 0.0, 1e-12, "T11 im");
	}

	#[test]
	fn from_z_invalid_zc_clamps_to_zref() {
		let z_re = [Z_REF];
		let z_im = [0.0];
		let m = MatchedS::from_z(&z_re, &z_im, 1, Z_REF, f64::NAN, 0.0);
		close(m.z0[0], Z_REF, 0.0, "z0");
		close(m.s_re[0], 0.0, 1e-12, "S11");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		let m0 = MatchedS::from_z(&z_re, &z_im, 1, Z_REF, 0.0, 0.0);
		close(m0.z0[0], Z_REF, 0.0, "z0 from 0");
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
	fn two_port_self_and_mutual_x_on_z() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let p_im = [0.0f32; 4];
		let x_nn = 10.0;
		let x_self = 5.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 2, Z_REF, &[0.0, 0.5], &[0.0, 0.0], x_nn, 2.0, 0.0, 0.0, 0.0, x_self,
		);
		close(m.r_im[0], x_self, 0.0, "X11");
		close(m.r_im[3], x_self, 0.0, "X22");
		close(m.r_im[1], x_nn, 1e-12, "X12");
		close(m.r_im[2], x_nn, 1e-12, "X21");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		assert!(m.residual < TAU, "residual {}", m.residual);
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
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		close(m.z0_im[1], 0.0, 0.0, "z0 im 1");
		assert!(m.residual < TAU, "residual {}", m.residual);
		let s11 = m.s_re[0].hypot(m.s_im[0]);
		let s22 = m.s_re[3].hypot(m.s_im[3]);
		assert!(s11 > 1e-3, "|S11| leftover reactance {s11}");
		assert!(s22 > 1e-3, "|S22| leftover reactance {s22}");
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
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
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
	fn n1_common_self_x_on_z() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let zc = 50.0;
		let x_self = 10.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 1, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, zc, x_self,
		);
		let z_re = 50.0;
		let z_im = x_self;
		let num_re = z_re - zc;
		let num_im = z_im;
		let den_re = z_re + zc;
		let den_im = z_im;
		let d2 = den_re * den_re + den_im * den_im;
		let s_re = (num_re * den_re + num_im * den_im) / d2;
		let s_im = (num_im * den_re - num_re * den_im) / d2;
		close(m.r_im[0], x_self, 0.0, "X11");
		close(m.s_re[0], s_re, 1e-12, "S11 re");
		close(m.s_im[0], s_im, 1e-12, "S11 im");
		close(m.z0[0], zc, 0.0, "z0 re");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
	}

	#[test]
	fn n1_per_port_self_x_leaves_sii() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let x_self = 10.0;
		let m = MatchedS::from_gram_coupled(
			&p_re, &p_im, 1, Z_REF, &[], &[], 0.0, 0.0, 0.0, 0.0, 0.0, x_self,
		);
		let r = 50.0;
		let den_re = 2.0 * r;
		let den_im = x_self;
		let d2 = den_re * den_re + den_im * den_im;
		let s_re = (x_self * den_im) / d2;
		let s_im = (x_self * den_re) / d2;
		close(m.z0[0], r, 1e-12, "z0 = R");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
		close(m.r_im[0], x_self, 0.0, "X11");
		close(m.s_re[0], s_re, 1e-12, "S11 re");
		close(m.s_im[0], s_im, 1e-12, "S11 im");
		assert!(m.s_re[0].hypot(m.s_im[0]) > 1e-3, "S11 from leftover X");
		assert!(m.residual < TAU, "residual {}", m.residual);
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
		assert_x_kernel(&m, &x, &y, x_nn, alpha, beta, 0.0, 0.0);
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
		assert_x_kernel(&m, &x, &y, x_nn, 2.0, 0.0, aniso, 0.0);
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
		assert_x_kernel(&m, &x, &y, x_nn, alpha, beta, aniso, 0.0);
		let o = n_vogel;
		assert!(
			(m.r_im[o * n + o + 1] - m.r_im[o * n + o + 2]).abs() > 1e-9,
			"equal-ρ x/y arms split"
		);
	}

	#[test]
	fn propagation_nn_sign_flip_and_common_z0() {
		let p_re = [0.5f32, 0.1, 0.05, 0.1, 0.5, 0.1, 0.05, 0.1, 0.5];
		let p_im = [0.0f32; 9];
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.0, 0.0];
		let x_nn = 10.0;
		let zc = 45.0;
		let m = MatchedS::from_gram_propagation(
			&p_re, &p_im, 3, Z_REF, &x, &y, x_nn, 0.0, 1.0, 1.0, 1.0, zc, 0.0,
		);
		close(m.r_im[1], x_nn, 1e-12, "X01 nn");
		close(m.r_im[3], x_nn, 1e-12, "X10");
		close(m.r_im[2], -x_nn, 1e-12, "X02 flip");
		close(m.r_re[1], 2.0 * Z_REF * 0.1, 1e-6, "R01 gram");
		close(m.z0[0], zc, 0.0, "z0");
		close(m.z0[1], zc, 0.0, "z0 1");
		close(m.z0_im[0], 0.0, 0.0, "z0 im");
	}

	#[test]
	fn propagation_eps_splits_xy_and_decay_shrinks() {
		let n = 3;
		let mut p_re = vec![0.05f32; n * n];
		for p in 0..n {
			p_re[p * n + p] = 0.5;
		}
		let p_im = vec![0.0f32; n * n];
		let x = [0.0f32, 0.5, 0.0];
		let y = [0.0f32, 0.0, 0.5];
		let x_nn = 8.0;
		let m = MatchedS::from_gram_propagation(
			&p_re, &p_im, n, Z_REF, &x, &y, x_nn, 0.0, 1.0, 4.0, 1.0, Z_REF, 0.0,
		);
		close(m.r_im[1], x_nn, 1e-12, "X01 is ref");
		assert!(
			(m.r_im[2] - m.r_im[1]).abs() > 1e-9,
			"εy≠εx splits equal-ρ arms"
		);
		close(m.r_im[2], m.r_im[2 * n], 1e-12, "X02 = X20");

		let md = MatchedS::from_gram_propagation(
			&p_re, &p_im, n, Z_REF, &[0.0, 0.5, 1.0], &[0.0, 0.0, 0.0], x_nn, 2.0, 1.0, 1.0, 1.0,
			Z_REF, 0.0,
		);
		close(md.r_im[1], x_nn, 1e-12, "decay nn");
		let far = x_nn * (-2.0_f64 * (1.0 - 0.5)).exp() * (-1.0);
		close(md.r_im[2], far, 1e-9, "decay far");
		assert!(md.r_im[2].abs() < md.r_im[1].abs(), "farther pair weaker");
	}

	#[test]
	fn propagation_invalid_zc_clamps_to_zref() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let m = MatchedS::from_gram_propagation(
			&p_re, &p_im, 1, Z_REF, &[0.0], &[0.0], 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
		);
		close(m.z0[0], Z_REF, 0.0, "z0 clamp");
	}

	fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
		a.iter()
			.zip(b)
			.map(|(x, y)| (x - y).abs())
			.fold(0.0, f64::max)
	}

	fn assert_st_close(faer: &StResult, legacy: &StResult, tol: f64, label: &str) {
		let ds = max_abs_diff(&faer.s_re, &legacy.s_re).max(max_abs_diff(&faer.s_im, &legacy.s_im));
		let dt = max_abs_diff(&faer.t_re, &legacy.t_re).max(max_abs_diff(&faer.t_im, &legacy.t_im));
		assert!(ds <= tol, "{label} max |S| diff {ds} > {tol}");
		assert!(dt <= tol, "{label} max |T| diff {dt} > {tol}");
	}

	fn rect_xy(nx: usize, ny: usize, dx: f32, dy: f32) -> (Vec<f32>, Vec<f32>) {
		let n = nx * ny;
		let mut x = vec![0.0f32; n];
		let mut y = vec![0.0f32; n];
		let mut k = 0;
		for ix in 0..nx {
			for iy in 0..ny {
				x[k] = dx * ix as f32;
				y[k] = dy * iy as f32;
				k += 1;
			}
		}
		(x, y)
	}

	fn green_z(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>, usize) {
		let (x, y) = rect_xy(nx, ny, 0.5, 0.5);
		let mut s = crate::prad::PradState::new();
		s.fill_green_pec_dipole_z(&x, &y, 1.0, 0.25, 0.1, 0.001);
		(s.r_re, s.r_im, s.n)
	}

	#[test]
	fn faer_agrees_with_legacy_small_herm() {
		let p_re = [0.5f32, 0.125, 0.125, 0.5];
		let n = 2;
		let scale = 2.0 * Z_REF;
		let r_re: Vec<f64> = p_re.iter().map(|v| scale * *v as f64).collect();
		let r_im = vec![0.0f64; 4];
		let z0 = vec![r_re[0].max(EPS_Z), r_re[3].max(EPS_Z)];
		let faer = form_s_and_t_with(&r_re, &r_im, &z0, n, Z_REF, crate::linalg::solve_herm_multi);
		let legacy = form_s_and_t_with(
			&r_re,
			&r_im,
			&z0,
			n,
			Z_REF,
			crate::legacy_linalg::solve_herm_multi,
		);
		assert_st_close(&faer, &legacy, 1e-9, "herm 2-port");
		let mut a = r_re.clone();
		for p in 0..n {
			a[p * n + p] += z0[p];
		}
		let fy = crate::linalg::inverse_diag_spd(&a, n);
		let ly = crate::legacy_linalg::inverse_diag_spd(&a, n);
		assert!(max_abs_diff(&fy, &ly) <= 1e-9, "inverse_diag_spd");
	}

	#[test]
	fn faer_agrees_with_legacy_small_cplx() {
		let z_re = [10.0, 1.0, 1.0, 12.0];
		let z_im = [3.0, 0.5, 0.5, 4.0];
		let n = 2;
		let z0_re = [10.0, 10.0];
		let z0_im = [0.0, 0.0];
		let faer = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::linalg::solve_cplx_multi,
		);
		let legacy = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::legacy_linalg::solve_cplx_multi,
		);
		assert_st_close(&faer, &legacy, 1e-9, "kurokawa 2-port");
		let mut a_re = z_re.to_vec();
		let mut a_im = z_im.to_vec();
		for p in 0..n {
			a_re[p * n + p] += z0_re[p];
			a_im[p * n + p] += z0_im[p];
		}
		let (fy_re, fy_im) = crate::linalg::inverse_diag_cplx(&a_re, &a_im, n);
		let (ly_re, ly_im) = crate::legacy_linalg::inverse_diag_cplx(&a_re, &a_im, n);
		assert!(
			max_abs_diff(&fy_re, &ly_re).max(max_abs_diff(&fy_im, &ly_im)) <= 1e-9,
			"inverse_diag_cplx"
		);
	}

	#[test]
	fn faer_agrees_with_legacy_green_8x8() {
		let (z_re, z_im, n) = green_z(8, 8);
		let z0_re = vec![Z_REF; n];
		let z0_im = vec![0.0; n];
		let faer = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::linalg::solve_cplx_multi,
		);
		let legacy = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::legacy_linalg::solve_cplx_multi,
		);
		assert_st_close(&faer, &legacy, 1e-9, "kurokawa 8x8 Green");
		let m = MatchedS::from_z(&z_re, &z_im, n, Z_REF, Z_REF, 0.0);
		let ds = max_abs_diff(&m.s_re, &faer.s_re).max(max_abs_diff(&m.s_im, &faer.s_im));
		let dt = max_abs_diff(&m.t_re, &faer.t_re).max(max_abs_diff(&m.t_im, &faer.t_im));
		assert!(ds <= 1e-9, "from_z S vs faer {ds}");
		assert!(dt <= 1e-9, "from_z T vs faer {dt}");
	}

	fn median_ms(mut samples: Vec<f64>) -> f64 {
		samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
		samples[samples.len() / 2]
	}

	fn time_ms(runs: usize, mut f: impl FnMut()) -> f64 {
		f();
		let mut samples = Vec::with_capacity(runs);
		for _ in 0..runs {
			let t0 = std::time::Instant::now();
			f();
			samples.push(t0.elapsed().as_secs_f64() * 1000.0);
		}
		median_ms(samples)
	}

	#[test]
	#[ignore]
	fn bench_linalg_32x32_legacy_vs_faer() {
		let (z_re, z_im, n) = green_z(32, 32);
		let z0_re = vec![Z_REF; n];
		let z0_im = vec![0.0; n];
		let faer = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::linalg::solve_cplx_multi,
		);
		let legacy = form_s_and_t_kurokawa_with(
			&z_re,
			&z_im,
			&z0_re,
			&z0_im,
			n,
			Z_REF,
			crate::legacy_linalg::solve_cplx_multi,
		);
		assert_st_close(&faer, &legacy, 1e-8, "kurokawa 32x32 Green");

		let mut a_re = z_re.clone();
		let a_im = z_im.clone();
		for p in 0..n {
			a_re[p * n + p] += Z_REF;
		}
		let g = Z_REF.max(EPS_Z).sqrt();
		let mut b_re = vec![0.0f64; n * n];
		let b_im = vec![0.0f64; n * n];
		for j in 0..n {
			b_re[j * n + j] = g;
		}

		let legacy_solve_ms = time_ms(3, || {
			let _ = crate::legacy_linalg::solve_cplx_multi(&a_re, &a_im, n, &b_re, &b_im);
		});
		let faer_solve_ms = time_ms(3, || {
			let _ = crate::linalg::solve_cplx_multi(&a_re, &a_im, n, &b_re, &b_im);
		});
		let legacy_full_ms = time_ms(3, || {
			let _ = form_s_and_t_kurokawa_with(
				&z_re,
				&z_im,
				&z0_re,
				&z0_im,
				n,
				Z_REF,
				crate::legacy_linalg::solve_cplx_multi,
			);
		});
		let faer_full_ms = time_ms(3, || {
			let _ = form_s_and_t_kurokawa_with(
				&z_re,
				&z_im,
				&z0_re,
				&z0_im,
				n,
				Z_REF,
				crate::linalg::solve_cplx_multi,
			);
		});
		println!(
			"\nN {n}\n  solve   legacy_ms {legacy_solve_ms:.1} | faer_ms {faer_solve_ms:.1} | speedup {solve_sp:.2}x\n  full S,T legacy_ms {legacy_full_ms:.1} | faer_ms {faer_full_ms:.1} | speedup {full_sp:.2}x",
			solve_sp = legacy_solve_ms / faer_solve_ms.max(1e-9),
			full_sp = legacy_full_ms / faer_full_ms.max(1e-9),
		);
	}
}
