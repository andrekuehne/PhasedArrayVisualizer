//! WP1 contract (`docs/green_function_plan.md` §5, §7). Finite short dipole
//! \(+\hat x\) over infinite PEC. No wasm-bindgen.
//!
//! ```text
//! z_pair_pec_dipole(dx, dy, h, ell, a, freq_scale) -> (re, im)  // ohms
//! f_iso_pec_dipole(theta, phi, h, ell, freq_scale)
//!     -> (e_th_re, e_th_im, e_ph_re, e_ph_im)
//! f_iso_pec_dipole_power(...) -> |Eθ|² + |Eφ|²
//! ```
//!
//! Positions and \((h,\ell,a)\) in wavelengths at \(f_0\).
//! \(k=2\pi\cdot\texttt{freq_scale}\), \(\eta_0=120\pi\).
//! Time convention \(e^{j\omega t}\), outgoing \(e^{-jkR}/R\), \(Z=R+jX\).
//! \(\theta=0\) is \(+\hat z\); \(\phi=0\) is \(+\hat x\).

use std::f64::consts::PI;

pub const ETA0: f64 = 120.0 * PI;
#[allow(dead_code)]
pub const DEFAULT_H: f64 = 0.25;
#[allow(dead_code)]
pub const DEFAULT_ELL: f64 = 0.1;
#[allow(dead_code)]
pub const DEFAULT_A: f64 = 0.001;

const R_FLOOR: f64 = 1e-30;

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

	fn scale(self, s: f64) -> Self {
		Self {
			re: self.re * s,
			im: self.im * s,
		}
	}

	fn mul(self, o: Self) -> Self {
		Self {
			re: self.re * o.re - self.im * o.im,
			im: self.re * o.im + self.im * o.re,
		}
	}

	fn sub(self, o: Self) -> Self {
		Self {
			re: self.re - o.re,
			im: self.im - o.im,
		}
	}

	fn tuple(self) -> (f64, f64) {
		(self.re, self.im)
	}

	#[allow(dead_code)]
	fn is_finite(self) -> bool {
		self.re.is_finite() && self.im.is_finite()
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

/// Free-space short-dipole \(Z_\mathrm{fs}=R_\mathrm{fs}+j X_\mathrm{fs}\) (no image).
pub(crate) fn z_fs(ell: f64, a: f64, freq_scale: f64) -> (f64, f64) {
	z_fs_cpx(ell, a, freq_scale).tuple()
}

fn z_fs_cpx(ell: f64, a: f64, freq_scale: f64) -> Cpx {
	let k = k_of(freq_scale);
	let kell = k * ell;
	let r_fs = ETA0 * kell * kell / (6.0 * PI);
	let x_fs = -ETA0 / (PI * kell) * ((ell / (2.0 * a)).ln() - 1.0);
	Cpx::new(r_fs, x_fs)
}

/// Hertzian \(Z(\mathbf{R})=-E_x\ell/I\) of an \(+\hat x\) moment \(I\ell\) at the
/// origin, \(I=1\). Full Balanis \(E_R,E_\theta\) (not far-field only). \(R>0\).
fn z_hertzian_x(dx: f64, dy: f64, dz: f64, ell: f64, freq_scale: f64) -> Cpx {
	let r2 = dx * dx + dy * dy + dz * dz;
	if !(r2.is_finite() && r2 > R_FLOOR * R_FLOOR) || !ell.is_finite() || !(freq_scale > 0.0) {
		return Cpx::nan();
	}
	let r = r2.sqrt();
	let k = k_of(freq_scale);
	let kr = k * r;
	if !(kr > 0.0 && kr.is_finite()) {
		return Cpx::nan();
	}

	let cos_th = dx / r;
	let sin_th = (dy * dy + dz * dz).sqrt() / r;
	let e_jkr = Cpx::new(kr.cos(), -kr.sin());
	// 1 + 1/(j kr) = 1 - j/kr
	let near = Cpx::new(1.0, -1.0 / kr);
	// 1 + 1/(j kr) - 1/(k²R²)
	let mid = Cpx::new(1.0 - 1.0 / (kr * kr), -1.0 / kr);
	// j * mid
	let j_mid = Cpx::new(-mid.im, mid.re);

	let e_r = near.mul(e_jkr).scale(ETA0 * ell / (2.0 * PI) * cos_th / (r * r));
	let e_th = j_mid
		.mul(e_jkr)
		.scale(ETA0 * k * ell / (4.0 * PI) * sin_th / r);
	let e_x = Cpx::new(
		e_r.re * cos_th - e_th.re * sin_th,
		e_r.im * cos_th - e_th.im * sin_th,
	);
	e_x.scale(-ell)
}

/// Mutual (or self when `dx == 0 && dy == 0`) impedance of two co-aligned
/// \(+\hat x\) short dipoles at the same height \(h\) over PEC. Ohms.
pub fn z_pair_pec_dipole(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
) -> (f64, f64) {
	if !geometry_ok(h, ell, a, freq_scale) || !dx.is_finite() || !dy.is_finite() {
		return Cpx::nan().tuple();
	}
	let z_img = z_hertzian_x(dx, dy, 2.0 * h, ell, freq_scale);
	let z = if dx == 0.0 && dy == 0.0 {
		z_fs_cpx(ell, a, freq_scale).sub(z_img)
	} else {
		z_hertzian_x(dx, dy, 0.0, ell, freq_scale).sub(z_img)
	};
	z.tuple()
}

/// Isolated space-wave \((E_\theta,E_\phi)\) of a unit-current short dipole plus
/// PEC image. \(e^{-jkr}/r\) stripped. Zero in the back hemisphere \(\theta>\pi/2\).
pub fn f_iso_pec_dipole(
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
	let ct = theta.cos();
	let amp = -ETA0 * k * ell / (2.0 * PI) * (k * h * ct).sin();
	let e_th = amp * ct * phi.cos();
	let e_ph = -amp * phi.sin();
	(e_th, 0.0, e_ph, 0.0)
}

/// \(\lvert F\rvert^2=\lvert E_\theta\rvert^2+\lvert E_\phi\rvert^2\).
pub fn f_iso_pec_dipole_power(theta: f64, phi: f64, h: f64, ell: f64, freq_scale: f64) -> f64 {
	let (et_re, et_im, ep_re, ep_im) = f_iso_pec_dipole(theta, phi, h, ell, freq_scale);
	et_re * et_re + et_im * et_im + ep_re * ep_re + ep_im * ep_im
}

#[cfg(test)]
fn z_fs_only(ell: f64, a: f64, freq_scale: f64) -> (f64, f64) {
	z_fs(ell, a, freq_scale)
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

	fn wrap_pi(x: f64) -> f64 {
		let t = 2.0 * PI;
		let mut y = (x + PI) % t;
		if y < 0.0 {
			y += t;
		}
		y - PI
	}

	#[test]
	fn coincident_uses_self_not_hertzian_origin() {
		let h = DEFAULT_H;
		let ell = DEFAULT_ELL;
		let a = DEFAULT_A;
		let fs = 1.0;
		let (re, im) = z_pair_pec_dipole(0.0, 0.0, h, ell, a, fs);
		assert!(re.is_finite() && im.is_finite(), "self NaN: {re} {im}");
		assert!(re > 0.0, "Re Z_self={re}");

		let origin = z_hertzian_x(0.0, 0.0, 0.0, ell, fs);
		assert!(!origin.is_finite(), "Hertzian at R=0 must not be used");

		let (r_fs, x_fs) = z_fs_only(ell, a, fs);
		let img = z_hertzian_x(0.0, 0.0, 2.0 * h, ell, fs);
		close(re, r_fs - img.re, 1e-12, "self re = R_fs - Re Z(0,0,2h)");
		close(im, x_fs - img.im, 1e-12, "self im = X_fs - Im Z(0,0,2h)");
	}

	#[test]
	fn z_fs_radiation_resistance() {
		let ell = DEFAULT_ELL;
		let a = DEFAULT_A;
		let fs = 1.0;
		let k = k_of(fs);
		let (r_fs, x_fs) = z_fs_only(ell, a, fs);
		close(
			r_fs,
			ETA0 * (k * ell) * (k * ell) / (6.0 * PI),
			1e-14,
			"R_fs",
		);
		assert!(x_fs.is_finite());
		assert!(x_fs < 0.0, "short-dipole X_fs capacitive, got {x_fs}");
	}

	#[test]
	fn reciprocity_and_anisotropy() {
		let h = DEFAULT_H;
		let ell = DEFAULT_ELL;
		let a = DEFAULT_A;
		let fs = 1.0;
		let (re, im) = z_pair_pec_dipole(0.4, -0.25, h, ell, a, fs);
		let (re_n, im_n) = z_pair_pec_dipole(-0.4, 0.25, h, ell, a, fs);
		close(re, re_n, 1e-12, "Z(Δ)=Z(-Δ) re");
		close(im, im_n, 1e-12, "Z(Δ)=Z(-Δ) im");
		assert!(re.is_finite() && im.is_finite());

		let (re_xy, im_xy) = z_pair_pec_dipole(0.5, 0.1, h, ell, a, fs);
		let (re_yx, im_yx) = z_pair_pec_dipole(0.1, 0.5, h, ell, a, fs);
		let d = (re_xy - re_yx).hypot(im_xy - im_yx);
		assert!(d > 0.05, "x-dipole is not isotropic: |Z(Δx,Δy)-Z(Δy,Δx)|={d}");

		let (re_yp, im_yp) = z_pair_pec_dipole(0.4, 0.25, h, ell, a, fs);
		let (re_yn, im_yn) = z_pair_pec_dipole(0.4, -0.25, h, ell, a, fs);
		let (re_xn, im_xn) = z_pair_pec_dipole(-0.4, 0.25, h, ell, a, fs);
		close(re_yp, re_yn, 1e-12, "even in Δy re");
		close(im_yp, im_yn, 1e-12, "even in Δy im");
		close(re_yp, re_xn, 1e-12, "even in Δx re");
		close(im_yp, im_xn, 1e-12, "even in Δx im");
	}

	#[test]
	fn invalid_geometry_is_nan() {
		let (re, im) = z_pair_pec_dipole(0.0, 0.0, 0.0, DEFAULT_ELL, DEFAULT_A, 1.0);
		assert!(re.is_nan() && im.is_nan(), "h=0");
		let (re, im) = z_pair_pec_dipole(0.2, 0.0, DEFAULT_H, 0.001, 0.001, 1.0);
		assert!(re.is_nan() && im.is_nan(), "ell <= 2a");
	}

	#[test]
	fn pattern_vanishes_as_h_to_zero() {
		let ell = DEFAULT_ELL;
		let fs = 1.0;
		let p = f_iso_pec_dipole_power(0.3, 0.4, 0.0, ell, fs);
		assert!(p.abs() < 1e-30, "|F|^2 at h=0: {p}");
		let p_small = f_iso_pec_dipole_power(0.4, 1.0, 1e-12, ell, fs);
		assert!(p_small.abs() < 1e-18, "|F|^2 at tiny h: {p_small}");
	}

	#[test]
	fn back_hemisphere_is_zero() {
		let (et_re, et_im, ep_re, ep_im) =
			f_iso_pec_dipole(PI * 0.5 + 0.2, 0.3, DEFAULT_H, DEFAULT_ELL, 1.0);
		assert_eq!((et_re, et_im, ep_re, ep_im), (0.0, 0.0, 0.0, 0.0));
		assert_eq!(
			f_iso_pec_dipole_power(2.0, 0.0, DEFAULT_H, DEFAULT_ELL, 1.0),
			0.0
		);
	}

	#[test]
	fn large_rho_phase_tracks_exp_mjk_rho() {
		let h = DEFAULT_H;
		let ell = DEFAULT_ELL;
		let a = DEFAULT_A;
		let fs = 1.0;
		let k = k_of(fs);
		let rho1 = 20.0;
		let rho2 = 22.0;
		// Broadside (ŷ): not on the dipole axis.
		let (re1, im1) = z_pair_pec_dipole(0.0, rho1, h, ell, a, fs);
		let (re2, im2) = z_pair_pec_dipole(0.0, rho2, h, ell, a, fs);
		assert!(re1.is_finite() && im1.is_finite() && re2.is_finite() && im2.is_finite());
		let arg1 = im1.atan2(re1);
		let arg2 = im2.atan2(re2);
		let d_arg = wrap_pi(arg2 - arg1);
		let expected = wrap_pi(-k * (rho2 - rho1));
		close(d_arg, expected, 0.05, "Δarg vs -k Δρ");
	}

	#[test]
	fn eh_planes_quarter_wave_height() {
		let h = 0.25;
		let ell = DEFAULT_ELL;
		let fs = 1.0;
		let k = k_of(fs);
		let pre = ETA0 * k * ell / (2.0 * PI);
		let th = PI * 0.25;

		let p_bore = f_iso_pec_dipole_power(0.0, 0.0, h, ell, fs);
		let p_e = f_iso_pec_dipole_power(th, 0.0, h, ell, fs);
		let p_h = f_iso_pec_dipole_power(th, PI * 0.5, h, ell, fs);

		let img0 = (k * h).sin();
		let img = (k * h * th.cos()).sin();
		close(p_bore, (pre * img0).powi(2), 1e-12, "boresight |F|^2");
		close(p_e, (pre * img * th.cos()).powi(2), 1e-12, "E-plane φ=0");
		close(p_h, (pre * img).powi(2), 1e-12, "H-plane φ=π/2");
		assert!(p_h > p_e, "H-plane lacks the extra cosθ of E-plane");
		assert!(p_bore > p_h, "λ/4 image factor peaks at boresight");

		let (et_re, et_im, ep_re, ep_im) = f_iso_pec_dipole(th, 0.0, h, ell, fs);
		close(et_im, 0.0, 0.0, "Eθ im (real far field)");
		close(ep_im, 0.0, 0.0, "Eφ im");
		close(ep_re, 0.0, 1e-15, "Eφ at φ=0");
		assert!(et_re.abs() > 0.0);
	}
}
