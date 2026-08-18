//! Radiative resistance, simultaneous real port match, and power-wave \(S\).
//!
//! §6–8 of `docs/approximate_matched_basis.md`. Operates on an already-formed
//! radiated-power Gram \(P_H\). Internals are f64; the Gram itself stays f32.

pub const Z_REF: f64 = 50.0;
pub const EPS_Z: f64 = 1e-9;
pub const K_MAX: u32 = 200;
pub const TAU: f64 = 1e-3;

const PIVOT_FLOOR: f64 = 1e-18;

/// Result of \(P_H \to R \to z_0 \to S\).
pub struct MatchedS {
	#[allow(dead_code)]
	pub n: usize,
	#[allow(dead_code)]
	pub z_ref: f64,
	pub z0: Vec<f64>,
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
	pub fn from_gram(p_re: &[f32], p_im: &[f32], n: usize, z_ref: f64) -> Self {
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
}

fn fill_rre_plus_d(r_re: &[f64], z0: &[f64], n: usize, out: &mut [f64]) {
	out.copy_from_slice(r_re);
	for p in 0..n {
		out[p * n + p] += z0[p];
	}
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

#[cfg(test)]
mod tests {
	use super::*;

	fn close(a: f64, b: f64, tol: f64, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
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
	fn n1_gram_gives_zref_and_s0() {
		let p_re = [0.5f32];
		let p_im = [0.0f32];
		let m = MatchedS::from_gram(&p_re, &p_im, 1, Z_REF);
		assert_eq!(m.n, 1);
		close(m.z_ref, Z_REF, 0.0, "z_ref");
		close(m.r_re[0], Z_REF, 1e-12, "R11");
		close(m.r_im[0], 0.0, 0.0, "R11 im");
		close(m.z0[0], Z_REF, 1e-9, "z0");
		close(m.s_re[0], 0.0, 1e-12, "S11 re");
		close(m.s_im[0], 0.0, 1e-12, "S11 im");
		close(m.t_re[0], 1.0, 1e-12, "T11");
		close(m.t_im[0], 0.0, 1e-12, "T11 im");
		assert!(m.residual < TAU, "residual {}", m.residual);
	}

	#[test]
	fn two_port_real_r_matches_and_s_symmetric() {
		// Isolated P_H so R = [[50, 10], [10, 50]].
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
}
