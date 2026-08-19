//! WP6 spectral PEC gate (`docs/green_function_plan.md` §12).
//!
//! Sommerfeld \(k_\rho\) integral of the free-space/PEC spectral Green function
//! for the same \(+\hat x\) Hertzian as WP1. Production still uses `green.rs`.
//! No \(\varepsilon_r\), no wasm-bindgen, no faer in the \(k_\rho\) loop.
//!
//! ```text
//! z_pair_pec_dipole_spectral(dx, dy, h, ell, a, freq_scale) -> (re, im)
//! f_iso_pec_dipole_spectral(theta, phi, h, ell, freq_scale) -> (Eθ, Eφ)
//! f_iso_pec_dipole_power_spectral(...) -> |F|²
//! ```
//!
//! Mutual \(Z\) is Hertzian (moment \(I\ell\)), matching WP1: no extra
//! \(\mathrm{sinc}(k_x\ell/2)\) line factor. Self is \(Z_\mathrm{fs}-Z(0,0,2h)\).
//! \(F^\mathrm{iso}\) is the saddle-point (space-wave) sample of the same
//! spectral \(\tilde E\), not a \(k_\rho\) integral.

use crate::green::{z_fs, ETA0};
use crate::quadrature::gauss_legendre;
use std::f64::consts::PI;

const R_FLOOR: f64 = 1e-30;

/// Quadrature for the Sommerfeld \(k_\rho\) integral.
#[derive(Clone, Copy, Debug)]
pub struct SpectralQuadConfig {
	/// Gauss–Legendre nodes on \(\alpha\in[0,\pi/2]\) (propagating, \(k_\rho=k\sin\alpha\)).
	pub n_k_prop: usize,
	/// Gauss–Legendre nodes per evanescent \(k_\rho\) panel.
	pub n_k_evan: usize,
	/// Truncation \(k_{\rho,\max}/k\) for \(z>0\) image terms.
	pub k_evan_max_over_k: f64,
	/// Bessel-lobe count for the coplanar (\(z=0\)) evanescent tail + Wynn.
	pub n_lobes: usize,
}

impl SpectralQuadConfig {
	pub const DEFAULT: Self = Self {
		n_k_prop: 48,
		n_k_evan: 32,
		k_evan_max_over_k: 12.0,
		n_lobes: 40,
	};

	pub const COARSE: Self = Self {
		n_k_prop: 16,
		n_k_evan: 16,
		k_evan_max_over_k: 6.0,
		n_lobes: 16,
	};

	pub const FINE: Self = Self {
		n_k_prop: 64,
		n_k_evan: 48,
		k_evan_max_over_k: 16.0,
		n_lobes: 56,
	};
}

impl Default for SpectralQuadConfig {
	fn default() -> Self {
		Self::DEFAULT
	}
}

#[derive(Clone, Copy)]
struct Cpx {
	re: f64,
	im: f64,
}

impl Cpx {
	fn new(re: f64, im: f64) -> Self {
		Self { re, im }
	}

	fn nan() -> Self {
		Self {
			re: f64::NAN,
			im: f64::NAN,
		}
	}

	fn zero() -> Self {
		Self { re: 0.0, im: 0.0 }
	}

	fn scale(self, s: f64) -> Self {
		Self {
			re: self.re * s,
			im: self.im * s,
		}
	}

	fn add(self, o: Self) -> Self {
		Self {
			re: self.re + o.re,
			im: self.im + o.im,
		}
	}

	fn sub(self, o: Self) -> Self {
		Self {
			re: self.re - o.re,
			im: self.im - o.im,
		}
	}

	fn mul(self, o: Self) -> Self {
		Self {
			re: self.re * o.re - self.im * o.im,
			im: self.re * o.im + self.im * o.re,
		}
	}

	fn inv(self) -> Self {
		let d = self.re * self.re + self.im * self.im;
		Self {
			re: self.re / d,
			im: -self.im / d,
		}
	}

	fn abs(self) -> f64 {
		self.re.hypot(self.im)
	}

	fn tuple(self) -> (f64, f64) {
		(self.re, self.im)
	}
}

fn k_of(freq_scale: f64) -> f64 {
	2.0 * PI * freq_scale
}

fn geometry_ok(h: f64, ell: f64, a: f64, freq_scale: f64) -> bool {
	h > 0.0
		&& a > 0.0
		&& ell > 2.0 * a
		&& freq_scale > 0.0
		&& h.is_finite()
		&& ell.is_finite()
		&& a.is_finite()
		&& freq_scale.is_finite()
}

/// \(J_0(x)\) (even). Series for \(|x|<20\), Hankel otherwise.
pub(crate) fn bessel_j0(x: f64) -> f64 {
	let ax = x.abs();
	if ax < 1e-14 {
		return 1.0;
	}
	if ax < 20.0 {
		let y = 0.25 * x * x;
		let mut term = 1.0;
		let mut sum = 1.0;
		for n in 1..80 {
			let nf = n as f64;
			term *= -y / (nf * nf);
			sum += term;
			if term.abs() <= 1e-18 * sum.abs().max(1.0) {
				break;
			}
		}
		return sum;
	}
	let z = 8.0 / ax;
	let y = z * z;
	let p0 = 1.0
		+ y * (-1.098628627e-3
			+ y * (2.734510407e-5 + y * (-2.073370639e-6 + y * 2.093887211e-7)));
	let q0 = -1.562499995e-2
		+ y * (1.430488765e-4
			+ y * (-6.911147651e-6 + y * (7.621095161e-7 - y * 9.34945152e-8)));
	let chi = ax - 0.7853981633974483;
	(p0 * chi.cos() - z * q0 * chi.sin()) * (0.6366197723675814 / ax).sqrt()
}

/// \(J_1(x)\) (odd).
fn bessel_j1(x: f64) -> f64 {
	let ax = x.abs();
	if ax < 1e-14 {
		return 0.5 * x;
	}
	if ax < 20.0 {
		let y = 0.25 * x * x;
		let mut term = 1.0;
		let mut sum = 1.0;
		for n in 1..80 {
			let nf = n as f64;
			term *= -y / (nf * (nf + 1.0));
			sum += term;
			if term.abs() <= 1e-18 * sum.abs().max(1.0) {
				break;
			}
		}
		return 0.5 * x * sum;
	}
	let z = 8.0 / ax;
	let y = z * z;
	let p1 = 1.0
		+ y * (1.83105e-3
			+ y * (-3.516396496e-5 + y * (2.457520174e-6 + y * (-2.40337019e-7))));
	let q1 = 4.687499995e-2
		+ y * (-2.002690873e-4
			+ y * (8.449199096e-6 + y * (-8.8228987e-7 + y * 1.05787412e-7)));
	let chi = ax - 2.356194490192345;
	let val = (p1 * chi.cos() - z * q1 * chi.sin()) * (0.6366197723675814 / ax).sqrt();
	if x < 0.0 {
		-val
	} else {
		val
	}
}

/// \(J_2(x)=(2/x)J_1(x)-J_0(x)\) (even).
pub(crate) fn bessel_j2(x: f64) -> f64 {
	let ax = x.abs();
	if ax < 1e-14 {
		return 0.0;
	}
	2.0 / x * bessel_j1(x) - bessel_j0(x)
}

/// Angular integral \(\int_0^{2\pi}(1-k_x^2/k^2)e^{-j\mathbf{k}\cdot\boldsymbol{\rho}}\,d\phi\)
/// \(= \pi\bigl[(2-t)J_0(\beta)+t\cos(2\phi_\rho)J_2(\beta)\bigr]\), \(t=k_\rho^2/k^2\).
fn ang_factor(k_rho: f64, k: f64, rho: f64, cos2phi: f64) -> f64 {
	let t = (k_rho / k) * (k_rho / k);
	let beta = k_rho * rho;
	let j0 = bessel_j0(beta);
	let j2 = bessel_j2(beta);
	(2.0 - t) * j0 + t * cos2phi * j2
}

/// \(\exp(-j\,k_z\,|z|)\). \(k_z\) real propagating or \(-j\alpha\) evanescent.
fn exp_neg_j_kz_z(kz_re: f64, kz_im: f64, abs_z: f64) -> Cpx {
	// exp(-j (kz_re + j kz_im) abs_z) = exp(kz_im abs_z) * (cos(kz_re abs_z) - j sin)
	let mag = (kz_im * abs_z).exp();
	let phase = kz_re * abs_z;
	Cpx::new(mag * phase.cos(), -mag * phase.sin())
}

pub(crate) fn map_gl(xi: f64, wi: f64, a: f64, b: f64) -> (f64, f64) {
	let half = 0.5 * (b - a);
	(half * xi + 0.5 * (a + b), half * wi)
}

/// Wynn \(\varepsilon\)-algorithm. Returns the last even-column Shanks estimate.
fn wynn_epsilon(partial: &[Cpx]) -> Cpx {
	let n = partial.len();
	if n == 0 {
		return Cpx::zero();
	}
	if n == 1 {
		return partial[0];
	}
	let mut prev = vec![Cpx::zero(); n];
	let mut curr = partial.to_vec();
	let mut best = curr[n - 1];
	for col in 1..n {
		let mut next = vec![Cpx::zero(); n - col];
		for i in 0..n - col {
			let d = curr[i + 1].sub(curr[i]);
			if d.abs() < 1e-30 {
				next[i] = curr[i + 1];
			} else {
				next[i] = prev[i + 1].add(d.inv());
			}
		}
		if col % 2 == 0 {
			best = next[next.len() - 1];
		}
		prev = curr;
		curr = next;
	}
	if best.re.is_finite() && best.im.is_finite() {
		best
	} else {
		partial[n - 1]
	}
}

/// Evanescent panel: \(k_\rho\in[a,b]\) with \(k\le a<b\), \(\beta=\mathrm{acosh}(k_\rho/k)\).
fn integrate_evan_beta(
	k: f64,
	a_krho: f64,
	b_krho: f64,
	abs_z: f64,
	rho: f64,
	cos2phi: f64,
	xi: &[f64],
	wi: &[f64],
) -> Cpx {
	let a = (a_krho / k).max(1.0).acosh();
	let b = (b_krho / k).max(1.0).acosh();
	if !(b > a) {
		return Cpx::zero();
	}
	let mut acc = Cpx::zero();
	for i in 0..xi.len() {
		let (beta, w_b) = map_gl(xi[i], wi[i], a, b);
		let sh = beta.sinh();
		let ch = beta.cosh();
		let k_rho = k * ch;
		let e = exp_neg_j_kz_z(0.0, -k * sh, abs_z);
		let ang = ang_factor(k_rho, k, rho, cos2phi);
		acc = acc.add(e.scale(w_b * k * ch * ang).mul(Cpx::new(0.0, 1.0)));
	}
	acc
}

/// Hertzian \(Z(\Delta x,\Delta y,\Delta z)=-E_x\ell/I\) from the Weyl/Sommerfeld
/// integral of the same moment \(I\ell\) as WP1. \(R>0\).
fn z_hertzian_spectral(
	dx: f64,
	dy: f64,
	dz: f64,
	ell: f64,
	freq_scale: f64,
	cfg: SpectralQuadConfig,
) -> Cpx {
	let abs_z = dz.abs();
	let rho = (dx * dx + dy * dy).sqrt();
	if !(rho.is_finite() && abs_z.is_finite()) {
		return Cpx::nan();
	}
	if rho < R_FLOOR && abs_z < R_FLOOR {
		return Cpx::nan();
	}
	let k = k_of(freq_scale);
	if !(k > 0.0 && k.is_finite() && ell.is_finite()) {
		return Cpx::nan();
	}
	let cos2phi = if rho < R_FLOOR {
		0.0
	} else {
		(dx * dx - dy * dy) / (rho * rho)
	};

	let pref = ETA0 * k * ell * ell / (8.0 * PI);
	let mut acc = Cpx::zero();

	let (xi_p, wi_p) = gauss_legendre(cfg.n_k_prop.max(2));
	for i in 0..xi_p.len() {
		let (alpha, w_a) = map_gl(xi_p[i], wi_p[i], 0.0, 0.5 * PI);
		let s = alpha.sin();
		let c = alpha.cos();
		let k_rho = k * s;
		let kz_re = k * c;
		let e = exp_neg_j_kz_z(kz_re, 0.0, abs_z);
		let ang = ang_factor(k_rho, k, rho, cos2phi);
		acc = acc.add(e.scale(w_a * k * s * ang));
	}

	let (xi_e, wi_e) = gauss_legendre(cfg.n_k_evan.max(2));
	let k_over = cfg.k_evan_max_over_k.max(1.001);
	if abs_z > 1e-12 {
		let beta_max = k_over.acosh();
		let n_pan = if rho < 1e-12 {
			1
		} else {
			(((k_over - 1.0) * k * rho / PI / 8.0).ceil() as usize)
				.max(1)
				.min(32)
		};
		let db = beta_max / n_pan as f64;
		for p in 0..n_pan {
			let ba = db * p as f64;
			let bb = ba + db;
			let a_k = k * ba.cosh();
			let b_k = k * bb.cosh();
			acc = acc.add(integrate_evan_beta(
				k, a_k, b_k, abs_z, rho, cos2phi, &xi_e, &wi_e,
			));
		}
	} else {
		let width = PI / rho.max(R_FLOOR);
		let n_lobes = cfg.n_lobes.max(4);
		let k_end = (k + n_lobes as f64 * width).max(k * k_over);
		let mut a = k;
		let mut sum = Cpx::zero();
		let mut partials: Vec<Cpx> = Vec::with_capacity(n_lobes);
		while a < k_end * (1.0 - 1e-14) {
			let b = (a + width).min(k_end);
			sum = sum.add(integrate_evan_beta(
				k, a, b, abs_z, rho, cos2phi, &xi_e, &wi_e,
			));
			partials.push(sum);
			a = b;
		}
		acc = acc.add(wynn_epsilon(&partials));
	}

	acc.scale(pref)
}

pub(crate) fn z_hertzian_spectral_tuple(
	dx: f64,
	dy: f64,
	dz: f64,
	ell: f64,
	freq_scale: f64,
	cfg: SpectralQuadConfig,
) -> (f64, f64) {
	z_hertzian_spectral(dx, dy, dz, ell, freq_scale, cfg).tuple()
}

/// Mutual (or self when `dx==0 && dy==0`) spectral PEC-dipole \(Z\). Ohms.
pub fn z_pair_pec_dipole_spectral(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
) -> (f64, f64) {
	z_pair_pec_dipole_spectral_cfg(dx, dy, h, ell, a, freq_scale, SpectralQuadConfig::DEFAULT)
}

pub fn z_pair_pec_dipole_spectral_cfg(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
	cfg: SpectralQuadConfig,
) -> (f64, f64) {
	if !geometry_ok(h, ell, a, freq_scale) || !dx.is_finite() || !dy.is_finite() {
		return Cpx::nan().tuple();
	}
	let z_img = z_hertzian_spectral(dx, dy, 2.0 * h, ell, freq_scale, cfg);
	let z = if dx == 0.0 && dy == 0.0 {
		let (r_fs, x_fs) = z_fs(ell, a, freq_scale);
		Cpx::new(r_fs, x_fs).sub(z_img)
	} else {
		z_hertzian_spectral(dx, dy, 0.0, ell, freq_scale, cfg).sub(z_img)
	};
	z.tuple()
}

/// Isolated space-wave \((E_\theta,E_\phi)\) from the saddle point of the same
/// spectral \(\tilde E\) used for \(Z\). Zero for \(\theta>\pi/2\).
pub fn f_iso_pec_dipole_spectral(
	theta: f64,
	phi: f64,
	h: f64,
	ell: f64,
	freq_scale: f64,
) -> (f64, f64, f64, f64) {
	if !(theta.is_finite()
		&& phi.is_finite()
		&& h.is_finite()
		&& ell.is_finite()
		&& freq_scale.is_finite())
	{
		return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
	}
	if theta > PI * 0.5 {
		return (0.0, 0.0, 0.0, 0.0);
	}
	let k = k_of(freq_scale);
	if !(k > 0.0 && k.is_finite()) {
		return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
	}
	let ct = theta.cos();
	let st = theta.sin();
	let cp = phi.cos();
	let sp = phi.sin();
	if ct.abs() < 1e-18 {
		// Horizon: Jacobian \(k\cos\theta\to 0\) and image \(\sin(kh\cos\theta)\to 0\).
		return (0.0, 0.0, 0.0, 0.0);
	}

	let kx = k * st * cp;
	let ky = k * st * sp;
	let kz = k * ct;
	let k2 = k * k;
	// Origin Hertzian spectral amplitude (1/(2π)² Weyl convention), I=1.
	let inv_kz = 1.0 / kz;
	let half_eta_k_ell = 0.5 * ETA0 * k * ell;
	let ex_o = -half_eta_k_ell * (1.0 - kx * kx / k2) * inv_kz;
	let ey_o = half_eta_k_ell * (kx * ky / k2) * inv_kz;
	let ez_o = half_eta_k_ell * (kx / k2);
	// PEC image: (e^{j kz h} - e^{-j kz h}) = 2j sin(kz h).
	let s_img = (kz * h).sin();
	let two_j_s = Cpx::new(0.0, 2.0 * s_img);
	let ex_t = two_j_s.scale(ex_o);
	let ey_t = two_j_s.scale(ey_o);
	let ez_t = two_j_s.scale(ez_o);
	// Far-field stripped: E = (-j k cosθ / 2π) Ẽ.
	let jac = Cpx::new(0.0, -k * ct / (2.0 * PI));
	let ex = jac.mul(ex_t);
	let ey = jac.mul(ey_t);
	let ez = jac.mul(ez_t);
	let e_th = Cpx::new(
		ex.re * ct * cp + ey.re * ct * sp - ez.re * st,
		ex.im * ct * cp + ey.im * ct * sp - ez.im * st,
	);
	let e_ph = Cpx::new(
		-ex.re * sp + ey.re * cp,
		-ex.im * sp + ey.im * cp,
	);
	(e_th.re, e_th.im, e_ph.re, e_ph.im)
}

pub fn f_iso_pec_dipole_power_spectral(
	theta: f64,
	phi: f64,
	h: f64,
	ell: f64,
	freq_scale: f64,
) -> f64 {
	let (et_re, et_im, ep_re, ep_im) = f_iso_pec_dipole_spectral(theta, phi, h, ell, freq_scale);
	et_re * et_re + et_im * et_im + ep_re * ep_re + ep_im * ep_im
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::green::{
		f_iso_pec_dipole, f_iso_pec_dipole_power, z_pair_pec_dipole, DEFAULT_A, DEFAULT_ELL,
		DEFAULT_H,
	};

	fn close(a: f64, b: f64, tol: f64, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	fn z_err(dx: f64, dy: f64, cfg: SpectralQuadConfig) -> (f64, f64, f64) {
		let (cr, ci) = z_pair_pec_dipole(dx, dy, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		let (sr, si) = z_pair_pec_dipole_spectral_cfg(
			dx,
			dy,
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			cfg,
		);
		let e = (sr - cr).hypot(si - ci);
		(cr.hypot(ci).max(1.0), e, e / cr.hypot(ci).max(1e-12))
	}

	#[test]
	fn spectral_matches_closed_self() {
		let (cr, ci) = z_pair_pec_dipole(0.0, 0.0, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		let (sr, si) = z_pair_pec_dipole_spectral(
			0.0,
			0.0,
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
		);
		assert!(sr.is_finite() && si.is_finite(), "spectral self NaN");
		assert!(sr > 0.0, "Re Z_self={sr}");
		let e = (sr - cr).hypot(si - ci);
		assert!(e < 1e-4, "self |ΔZ|={e} closed=({cr},{ci}) spec=({sr},{si})");
	}

	#[test]
	fn spectral_matches_closed_mutual_grid() {
		let lags = [
			(0.1, 0.0),
			(0.0, 0.1),
			(0.5, 0.0),
			(0.0, 0.5),
			(0.5, 0.5),
			(0.3, 0.1),
			(0.1, 0.3),
			(1.0, 0.0),
			(0.0, 2.0),
			(5.0, 0.0),
		];
		let mut max_abs = 0.0f64;
		let mut max_rel = 0.0f64;
		for (dx, dy) in lags {
			let (scale, e, rel) = z_err(dx, dy, SpectralQuadConfig::DEFAULT);
			assert!(
				e < 1e-4 || rel < 1e-3,
				"lag ({dx},{dy}) |ΔZ|={e} rel={rel} scale={scale}"
			);
			max_abs = max_abs.max(e);
			max_rel = max_rel.max(rel);
		}
		assert!(max_abs.is_finite());
		let _ = max_rel;
	}

	#[test]
	fn spectral_matches_closed_reciprocity() {
		let (re, im) = z_pair_pec_dipole_spectral(0.4, -0.25, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		let (re_n, im_n) =
			z_pair_pec_dipole_spectral(-0.4, 0.25, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		close(re, re_n, 1e-9, "Z(Δ)=Z(-Δ) re");
		close(im, im_n, 1e-9, "Z(Δ)=Z(-Δ) im");
		let (re_xy, im_xy) =
			z_pair_pec_dipole_spectral(0.5, 0.1, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		let (re_yx, im_yx) =
			z_pair_pec_dipole_spectral(0.1, 0.5, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		let d = (re_xy - re_yx).hypot(im_xy - im_yx);
		assert!(d > 0.05, "x-dipole is not isotropic: |Δ|={d}");
	}

	#[test]
	fn spectral_matches_closed_pattern() {
		let h = 0.25;
		let ell = DEFAULT_ELL;
		let fs = 1.0;
		for &(th, ph) in &[
			(0.0, 0.0),
			(PI * 0.25, 0.0),
			(PI * 0.25, 0.5 * PI),
			(0.4, 0.7),
		] {
			let (ctr, cti, cpr, cpi) = f_iso_pec_dipole(th, ph, h, ell, fs);
			let (str, sti, spr, spi) = f_iso_pec_dipole_spectral(th, ph, h, ell, fs);
			close(str, ctr, 1e-12, "Eθ re");
			close(sti, cti, 1e-12, "Eθ im");
			close(spr, cpr, 1e-12, "Eφ re");
			close(spi, cpi, 1e-12, "Eφ im");
			let pc = f_iso_pec_dipole_power(th, ph, h, ell, fs);
			let ps = f_iso_pec_dipole_power_spectral(th, ph, h, ell, fs);
			close(ps, pc, 1e-12, "|F|²");
		}
		let (et_re, et_im, ep_re, ep_im) =
			f_iso_pec_dipole_spectral(PI * 0.5 + 0.2, 0.3, h, ell, fs);
		assert_eq!((et_re, et_im, ep_re, ep_im), (0.0, 0.0, 0.0, 0.0));
		assert_eq!(
			f_iso_pec_dipole_power_spectral(2.0, 0.0, h, ell, fs),
			0.0
		);
	}

	#[test]
	fn invalid_geometry_is_nan() {
		let (re, im) = z_pair_pec_dipole_spectral(0.0, 0.0, 0.0, DEFAULT_ELL, DEFAULT_A, 1.0);
		assert!(re.is_nan() && im.is_nan(), "h=0");
		let (re, im) = z_pair_pec_dipole_spectral(0.2, 0.0, DEFAULT_H, 0.001, 0.001, 1.0);
		assert!(re.is_nan() && im.is_nan(), "ell <= 2a");
	}

	#[test]
	fn quadrature_convergence() {
		let lags = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)];
		for (dx, dy) in lags {
			let (_, e_c, rel_c) = z_err(dx, dy, SpectralQuadConfig::COARSE);
			let (_, e_d, rel_d) = z_err(dx, dy, SpectralQuadConfig::DEFAULT);
			let (_, e_f, rel_f) = z_err(dx, dy, SpectralQuadConfig::FINE);
			assert!(
				e_d <= e_c * 1.05 + 1e-6 || rel_d <= rel_c,
				"default should not be worse than coarse at ({dx},{dy}): {e_c} vs {e_d}"
			);
			assert!(
				e_f < 1e-4 || rel_f < 1e-3,
				"fine lag ({dx},{dy}) |ΔZ|={e_f} rel={rel_f}"
			);
		}
	}

	#[test]
	fn bessel_j0_j2_known_values() {
		close(bessel_j0(0.0), 1.0, 1e-15, "J0(0)");
		close(bessel_j0(2.404825557695773), 0.0, 1e-10, "J0 first zero");
		close(bessel_j2(0.0), 0.0, 1e-15, "J2(0)");
		close(bessel_j1(0.0), 0.0, 1e-15, "J1(0)");
		close(bessel_j0(1.0), 0.7651976865579666, 1e-12, "J0(1)");
		close(bessel_j2(1.0), 0.1149034849319005, 1e-12, "J2(1)");
		close(bessel_j0(25.0), bessel_j0(-25.0), 0.0, "J0 even");
		close(bessel_j2(25.0), bessel_j2(-25.0), 1e-14, "J2 even");
	}

	#[ignore]
	#[test]
	fn bench_spectral_vs_closed_per_lag() {
		use std::time::Instant;
		let lags = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (5.0, 0.0)];
		for (dx, dy) in lags {
			let t0 = Instant::now();
			let mut c = (0.0, 0.0);
			for _ in 0..200 {
				c = z_pair_pec_dipole(dx, dy, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
			}
			let closed_us = t0.elapsed().as_secs_f64() * 1e6 / 200.0;
			let t1 = Instant::now();
			let mut s = (0.0, 0.0);
			for _ in 0..8 {
				s = z_pair_pec_dipole_spectral(dx, dy, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
			}
			let spec_us = t1.elapsed().as_secs_f64() * 1e6 / 8.0;
			let e = (s.0 - c.0).hypot(s.1 - c.1);
			eprintln!(
				"lag ({dx},{dy}): closed {closed_us:.3} µs | spectral {spec_us:.1} µs | |ΔZ|={e:.3e} Z_c=({:.4},{:.4}) Z_s=({:.4},{:.4})",
				c.0, c.1, s.0, s.1
			);
		}
	}
}
