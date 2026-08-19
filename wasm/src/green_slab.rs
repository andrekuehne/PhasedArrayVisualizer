//! WP7 grounded-slab Green function (`docs/green_function_plan.md` §13).
//!
//! Horizontal short dipole (\(+\hat x\)) in air at height \(h\) over a grounded
//! dielectric slab. Same \(Z_\mathrm{pair}\)/\(F^\mathrm{iso}\) interface as WP1.
//! Bindgen lives in `element.rs` / `lib.rs`. Closed-form PEC stays `green.rs`.
//!
//! Stack: air \(z>0\); dipole at \(z=h\); dielectric \(-h_\mathrm{sub}<z<0\)
//! with \(\varepsilon=\varepsilon_r(1-j\tan\delta)\); PEC at \(z=-h_\mathrm{sub}\).
//!
//! ```text
//! z_pair_slab_dipole(dx, dy, h, ell, a, freq_scale, env) -> (re, im)
//! f_iso_slab_dipole(theta, phi, h, ell, freq_scale, env) -> (Eθ, Eφ)
//! f_iso_slab_dipole_power(...) -> |F|²
//! slab_dipole_power_budget(h, ell, a, freq_scale, env) -> SlabPowerBudget
//! ```
//!
//! Mutual \(Z\) is Hertzian (moment \(I\ell\)) plus slab reflection. Self is
//! \(Z_\mathrm{fs}+Z_\mathrm{refl}\). \(F^\mathrm{iso}\) is the space-wave saddle
//! of the same spectral \(\tilde E\) (no surface-wave term). Power split:
//! \(\operatorname{Re}Z_\mathrm{self}=P_\mathrm{rad}+P_\mathrm{sw}+P_\mathrm{diss}\)
//! at \(|I|=1\,\mathrm{A}\) (ohm ≡ watt in this current basis). EIRP stays on
//! \(P_\mathrm{stimulated}\) (WP7b).

use crate::green::{z_fs, ETA0};
use crate::green_spectral::{
	bessel_j0, bessel_j2, map_gl, z_hertzian_spectral_tuple, SpectralQuadConfig,
};
use crate::quadrature::gauss_legendre;
use std::f64::consts::PI;

const R_FLOOR: f64 = 1e-30;
const KP_IMAG_FLOOR: f64 = 1e-14;

/// Grounded-slab environment. Lengths in wavelengths at \(f_0\).
#[derive(Clone, Copy, Debug)]
pub struct SlabEnv {
	pub eps_r: f64,
	pub h_sub: f64,
	pub tan_delta: f64,
}

impl SlabEnv {
	pub const DEFAULT: Self = Self {
		eps_r: 10.0,
		h_sub: 0.05,
		tan_delta: 0.0,
	};

	/// Thin slab used as a PEC-recovery check.
	pub const PEC_LIMIT: Self = Self {
		eps_r: 12.0,
		h_sub: 1e-6,
		tan_delta: 0.0,
	};

	fn eps_c(self) -> Cpx {
		Cpx::new(self.eps_r, -self.eps_r * self.tan_delta.max(0.0))
	}

	pub fn ok(self) -> bool {
		self.eps_r >= 1.0
			&& self.h_sub > 0.0
			&& self.tan_delta >= 0.0
			&& self.eps_r.is_finite()
			&& self.h_sub.is_finite()
			&& self.tan_delta.is_finite()
	}
}

/// Isolated-element power split at unit current, ohms \(\equiv\) watts.
#[derive(Clone, Copy, Debug)]
pub struct SlabPowerBudget {
	pub re_z_self: f64,
	pub p_rad: f64,
	pub p_sw: f64,
	pub p_diss: f64,
	pub closure_residual: f64,
}

#[derive(Clone, Copy, Debug)]
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

	fn from_real(re: f64) -> Self {
		Self { re, im: 0.0 }
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

	fn mul_j(self) -> Self {
		Self {
			re: -self.im,
			im: self.re,
		}
	}

	fn div(self, o: Self) -> Self {
		let d = o.re * o.re + o.im * o.im;
		Self {
			re: (self.re * o.re + self.im * o.im) / d,
			im: (self.im * o.re - self.re * o.im) / d,
		}
	}

	fn abs(self) -> f64 {
		self.re.hypot(self.im)
	}

	fn is_finite(self) -> bool {
		self.re.is_finite() && self.im.is_finite()
	}

	/// Principal square root, then flipped so \(\operatorname{Im}\le 0\).
	fn sqrt_kz(self) -> Self {
		let r = self.abs();
		let t = ((r + self.re) * 0.5).max(0.0).sqrt();
		let u = ((r - self.re) * 0.5).max(0.0).sqrt();
		let mut s = if self.im >= 0.0 {
			Self { re: t, im: u }
		} else {
			Self { re: t, im: -u }
		};
		if s.im > 0.0 || (s.im.abs() <= 1e-18 && s.re < 0.0) {
			s = s.scale(-1.0);
		}
		s
	}

	fn sin(self) -> Self {
		Self {
			re: self.re.sin() * self.im.cosh(),
			im: self.re.cos() * self.im.sinh(),
		}
	}

	fn cos(self) -> Self {
		Self {
			re: self.re.cos() * self.im.cosh(),
			im: -self.re.sin() * self.im.sinh(),
		}
	}

	fn log(self) -> Self {
		Self {
			re: self.abs().ln(),
			im: self.im.atan2(self.re),
		}
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

fn kz_air(k_rho: f64, k: f64) -> Cpx {
	let disc = k * k - k_rho * k_rho;
	if disc >= 0.0 {
		Cpx::new(disc.sqrt(), 0.0)
	} else {
		Cpx::new(0.0, -(-disc).sqrt())
	}
}

fn kz_dielectric(k_rho: f64, k: f64, eps_c: Cpx) -> Cpx {
	let disc = eps_c.scale(k * k).sub(Cpx::from_real(k_rho * k_rho));
	disc.sqrt_kz()
}

fn kz_air_cpx(k_rho: Cpx, k: f64) -> Cpx {
	Cpx::from_real(k * k).sub(k_rho.mul(k_rho)).sqrt_kz()
}

fn kz_dielectric_cpx(k_rho: Cpx, k: f64, eps_c: Cpx) -> Cpx {
	eps_c.scale(k * k).sub(k_rho.mul(k_rho)).sqrt_kz()
}

/// \(\Gamma_\mathrm{TE}=(j k_{z0}\sin\theta-k_{zd}\cos\theta)/(j k_{z0}\sin\theta+k_{zd}\cos\theta)\).
fn gamma_te(kz_a: Cpx, kz_d: Cpx, h_sub: f64) -> Cpx {
	let th = kz_d.scale(h_sub);
	let s = th.sin();
	let c = th.cos();
	let js = kz_a.mul(s).mul_j();
	let kc = kz_d.mul(c);
	let den = js.add(kc);
	if den.abs() < 1e-30 {
		return Cpx::nan();
	}
	js.sub(kc).div(den)
}

/// \(\Gamma_\mathrm{TM}=(j k_{zd}\sin\theta-k_{z0}\varepsilon\cos\theta)/(\cdots+)\).
fn gamma_tm(kz_a: Cpx, kz_d: Cpx, eps_c: Cpx, h_sub: f64) -> Cpx {
	let th = kz_d.scale(h_sub);
	let s = th.sin();
	let c = th.cos();
	let js = kz_d.mul(s).mul_j();
	let ec = kz_a.mul(eps_c).mul(c);
	let den = js.add(ec);
	if den.abs() < 1e-30 {
		return Cpx::nan();
	}
	js.sub(ec).div(den)
}

fn tm_den(k_rho: Cpx, k: f64, eps_c: Cpx, h_sub: f64) -> Cpx {
	let kz_a = kz_air_cpx(k_rho, k);
	let kz_d = kz_dielectric_cpx(k_rho, k, eps_c);
	let th = kz_d.scale(h_sub);
	kz_d.mul(th.sin()).mul_j().add(kz_a.mul(eps_c).mul(th.cos()))
}

fn tm_num(k_rho: Cpx, k: f64, eps_c: Cpx, h_sub: f64) -> Cpx {
	let kz_a = kz_air_cpx(k_rho, k);
	let kz_d = kz_dielectric_cpx(k_rho, k, eps_c);
	let th = kz_d.scale(h_sub);
	kz_d.mul(th.sin()).mul_j().sub(kz_a.mul(eps_c).mul(th.cos()))
}

fn tm0_f_lossless(k_rho: f64, k: f64, eps_r: f64, h_sub: f64) -> f64 {
	let alpha = (k_rho * k_rho - k * k).max(0.0).sqrt();
	let kz_d = (eps_r * k * k - k_rho * k_rho).max(0.0).sqrt();
	let theta = kz_d * h_sub;
	kz_d * theta.tan() - eps_r * alpha
}

/// Real-axis TM0 root in \((k,k\sqrt{\varepsilon_r})\) with \(\theta<\pi/2\).
fn tm0_pole_real(k: f64, env: SlabEnv) -> Option<f64> {
	if !(env.eps_r > 1.0 + 1e-12 && env.h_sub > 0.0 && k > 0.0) {
		return None;
	}
	let k_hi = k * env.eps_r.sqrt();
	let kz_cap = 0.5 * PI / env.h_sub;
	let kd_max2 = env.eps_r * k * k;
	let k_lo = if kz_cap * kz_cap < kd_max2 {
		(kd_max2 - kz_cap * kz_cap).max(k * k).sqrt()
	} else {
		k
	};
	let lo0 = (k_lo + 1e-12 * k).min(k_hi * (1.0 - 1e-12));
	let hi0 = k_hi * (1.0 - 1e-12);
	if !(hi0 > lo0) {
		return None;
	}
	let f_lo = tm0_f_lossless(lo0, k, env.eps_r, env.h_sub);
	let f_hi = tm0_f_lossless(hi0, k, env.eps_r, env.h_sub);
	if !f_lo.is_finite() || !f_hi.is_finite() {
		return None;
	}
	// Root from the dielectric-light-line side (TM0, \(\theta\to 0\)).
	let mut lo = lo0;
	let mut hi = hi0;
	if f_lo * f_hi > 0.0 {
		// Scan for a sign change (thin slabs still have one).
		let n = 80;
		let mut found = false;
		let mut a = hi0;
		for i in 0..n {
			let t = (i as f64) / n as f64;
			let b = hi0 + (lo0 - hi0) * t;
			let fa = tm0_f_lossless(a, k, env.eps_r, env.h_sub);
			let fb = tm0_f_lossless(b, k, env.eps_r, env.h_sub);
			if fa.is_finite() && fb.is_finite() && fa * fb <= 0.0 {
				lo = a.min(b);
				hi = a.max(b);
				found = true;
				break;
			}
			a = b;
		}
		if !found {
			return None;
		}
	} else if f_hi.signum() == f_lo.signum() {
		return None;
	}
	for _ in 0..80 {
		let mid = 0.5 * (lo + hi);
		let fm = tm0_f_lossless(mid, k, env.eps_r, env.h_sub);
		if !fm.is_finite() {
			break;
		}
		let fl = tm0_f_lossless(lo, k, env.eps_r, env.h_sub);
		if fl * fm <= 0.0 {
			hi = mid;
		} else {
			lo = mid;
		}
	}
	Some(0.5 * (lo + hi))
}

fn tm0_pole(k: f64, env: SlabEnv) -> Option<Cpx> {
	let k0 = tm0_pole_real(k, env)?;
	let mut kp = Cpx::new(k0, -KP_IMAG_FLOOR);
	let eps_c = env.eps_c();
	if env.tan_delta > 1e-16 {
		for _ in 0..24 {
			let d = tm_den(kp, k, eps_c, env.h_sub);
			let h = (1e-8 * k).max(1e-10);
			let d_p = tm_den(Cpx::new(kp.re + h, kp.im), k, eps_c, env.h_sub);
			let d_m = tm_den(Cpx::new(kp.re - h, kp.im), k, eps_c, env.h_sub);
			let dp = d_p.sub(d_m).scale(0.5 / h);
			if dp.abs() < 1e-30 {
				break;
			}
			let nxt = kp.sub(d.div(dp));
			if !nxt.is_finite() {
				break;
			}
			kp = nxt;
			if d.abs() < 1e-12 * k {
				break;
			}
		}
	}
	if kp.re > k && kp.re < k * env.eps_r.sqrt() * 1.01 && kp.is_finite() {
		if kp.im > 0.0 {
			kp.im = -kp.im.abs();
		}
		Some(kp)
	} else {
		None
	}
}

fn gamma_tm_residue(kp: Cpx, k: f64, env: SlabEnv) -> Cpx {
	let eps_c = env.eps_c();
	let n = tm_num(kp, k, eps_c, env.h_sub);
	let h = (1e-8 * k).max(1e-10);
	let d_p = tm_den(Cpx::new(kp.re + h, kp.im), k, eps_c, env.h_sub);
	let d_m = tm_den(Cpx::new(kp.re - h, kp.im), k, eps_c, env.h_sub);
	let dp = d_p.sub(d_m).scale(0.5 / h);
	n.div(dp)
}

fn a_te(k_rho: f64, rho: f64, cos2phi: f64) -> f64 {
	let j0 = bessel_j0(k_rho * rho);
	let j2 = bessel_j2(k_rho * rho);
	j0 + cos2phi * j2
}

fn a_tm(k_rho: f64, rho: f64, cos2phi: f64) -> f64 {
	let j0 = bessel_j0(k_rho * rho);
	let j2 = bessel_j2(k_rho * rho);
	j0 - cos2phi * j2
}

fn exp_neg_j_kz_z(kz: Cpx, abs_z: f64) -> Cpx {
	// exp(-j kz z) = exp(kz.im * z) * (cos(kz.re z) - j sin)
	let mag = (kz.im * abs_z).exp();
	let phase = kz.re * abs_z;
	Cpx::new(mag * phase.cos(), -mag * phase.sin())
}

/// Spectral density \(G(k_\rho)\) so \(Z_\mathrm{refl}=\mathrm{pref}\int G\,dk_\rho\).
fn refl_g(k_rho: f64, k: f64, h: f64, rho: f64, cos2phi: f64, env: SlabEnv) -> Cpx {
	if !(k_rho >= 0.0 && k_rho.is_finite()) {
		return Cpx::zero();
	}
	let kz_a = kz_air(k_rho, k);
	let kz_d = kz_dielectric(k_rho, k, env.eps_c());
	if kz_a.abs() < 1e-18 {
		return Cpx::zero();
	}
	let g_te = gamma_te(kz_a, kz_d, env.h_sub);
	let g_tm = gamma_tm(kz_a, kz_d, env.eps_c(), env.h_sub);
	if !g_te.is_finite() || !g_tm.is_finite() {
		return Cpx::zero();
	}
	let kz2_over_k2 = kz_a.mul(kz_a).scale(1.0 / (k * k));
	let br = g_te
		.scale(a_te(k_rho, rho, cos2phi))
		.add(g_tm.mul(kz2_over_k2).scale(a_tm(k_rho, rho, cos2phi)));
	let e = exp_neg_j_kz_z(kz_a, 2.0 * h);
	br.mul(e).scale(k_rho).div(kz_a)
}

fn g_residue(kp: Cpx, k: f64, h: f64, rho: f64, cos2phi: f64, env: SlabEnv) -> Cpx {
	let k_rho = kp.re;
	let kz_a = kz_air(k_rho, k);
	let kz2_over_k2 = kz_a.mul(kz_a).scale(1.0 / (k * k));
	let res_g = gamma_tm_residue(kp, k, env);
	let e = exp_neg_j_kz_z(kz_a, 2.0 * h);
	kz2_over_k2
		.scale(a_tm(k_rho, rho, cos2phi))
		.mul(res_g)
		.mul(e)
		.scale(k_rho)
		.div(kz_a)
}

#[derive(Clone, Copy)]
struct ReflSplit {
	total: Cpx,
	prop: Cpx,
	sw_jump: Cpx,
}

fn integrate_refl(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	freq_scale: f64,
	env: SlabEnv,
	cfg: SpectralQuadConfig,
) -> ReflSplit {
	let nan = ReflSplit {
		total: Cpx::nan(),
		prop: Cpx::nan(),
		sw_jump: Cpx::zero(),
	};
	let rho = (dx * dx + dy * dy).sqrt();
	if !(rho.is_finite() && h > 0.0 && env.ok()) {
		return nan;
	}
	let k = k_of(freq_scale);
	if !(k > 0.0 && k.is_finite() && ell.is_finite()) {
		return nan;
	}
	let cos2phi = if rho < R_FLOOR {
		0.0
	} else {
		(dx * dx - dy * dy) / (rho * rho)
	};
	let pref = ETA0 * k * ell * ell / (8.0 * PI);
	let mut kp = tm0_pole(k, env);
	if let Some(p) = kp {
		if p.re - k < 1e-3 * k {
			kp = None;
		}
	}
	let res_g = kp.map(|p| g_residue(p, k, h, rho, cos2phi, env));

	let mut prop = Cpx::zero();
	let (xi_p, wi_p) = gauss_legendre(cfg.n_k_prop.max(2));
	for i in 0..xi_p.len() {
		let (alpha, w_a) = map_gl(xi_p[i], wi_p[i], 0.0, 0.5 * PI);
		let k_rho = k * alpha.sin();
		let g = refl_g(k_rho, k, h, rho, cos2phi, env);
		// G dkρ = G * kz dα; G already has kρ/kz, so G*kz dα = bracket exp kρ dα.
		let kz = kz_air(k_rho, k);
		prop = prop.add(g.mul(kz).scale(w_a));
	}

	let mut evan = Cpx::zero();
	let (xi_e, wi_e) = gauss_legendre(cfg.n_k_evan.max(2));
	let k_over = cfg.k_evan_max_over_k.max(1.001);
	let k_end = k * k_over;
	let k_diel = k * env.eps_r.max(1.001).sqrt();
	let mut edges = vec![k, k_end];
	if k_diel > k && k_diel < k_end {
		edges.push(k_diel);
	}
	if let Some(p) = kp {
		if p.re > k && p.re < k_end {
			edges.push(p.re);
		}
	}
	edges.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
	edges.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * k);

	for w in edges.windows(2) {
		let a_k = w[0];
		let b_k = w[1];
		if !(b_k > a_k) {
			continue;
		}
		let n_pan = ((((b_k - a_k) / k) * k * (2.0 * h) / PI * 4.0).ceil() as usize)
			.max(1)
			.min(24);
		let da = (b_k - a_k) / n_pan as f64;
		for p in 0..n_pan {
			let aa = a_k + da * p as f64;
			let bb = aa + da;
			let ba = (aa / k).max(1.0).acosh();
			let bb_b = (bb / k).max(1.0).acosh();
			if !(bb_b > ba) {
				continue;
			}
			for i in 0..xi_e.len() {
				let (beta, w_b) = map_gl(xi_e[i], wi_e[i], ba, bb_b);
				let k_rho = k * beta.cosh();
				let g = refl_g(k_rho, k, h, rho, cos2phi, env);
				let dkr_db = k * beta.sinh();
				let mut val = g.scale(dkr_db);
				if let (Some(p0), Some(rg)) = (kp, res_g) {
					let den = Cpx::new(k_rho - p0.re, -p0.im);
					if den.abs() > 1e-18 {
						val = val.sub(rg.div(den).scale(dkr_db));
					}
				}
				evan = evan.add(val.scale(w_b));
			}
		}
	}

	let mut residue = Cpx::zero();
	let mut sw_jump = Cpx::zero();
	if let (Some(p0), Some(rg)) = (kp, res_g) {
		let num = Cpx::new(k_end - p0.re, -p0.im);
		let den = Cpx::new(k - p0.re, -p0.im);
		if num.abs() > 0.0 && den.abs() > 0.0 {
			residue = rg.mul(num.div(den).log());
		}
		// Indentation jump \(-j\pi\operatorname{Res}\) (SW launched power).
		sw_jump = rg.mul_j().scale(-PI);
	}

	let total = prop.add(evan).add(residue);
	ReflSplit {
		total: total.scale(pref),
		prop: prop.scale(pref),
		sw_jump: sw_jump.scale(pref),
	}
}

pub fn z_pair_slab_dipole(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
	env: SlabEnv,
) -> (f64, f64) {
	z_pair_slab_dipole_cfg(
		dx,
		dy,
		h,
		ell,
		a,
		freq_scale,
		env,
		SpectralQuadConfig::DEFAULT,
	)
}

pub fn z_pair_slab_dipole_cfg(
	dx: f64,
	dy: f64,
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
	env: SlabEnv,
	cfg: SpectralQuadConfig,
) -> (f64, f64) {
	if !geometry_ok(h, ell, a, freq_scale) || !dx.is_finite() || !dy.is_finite() || !env.ok() {
		return Cpx::nan().tuple();
	}
	let z_r = integrate_refl(dx, dy, h, ell, freq_scale, env, cfg).total;
	let z = if dx == 0.0 && dy == 0.0 {
		let (r_fs, x_fs) = z_fs(ell, a, freq_scale);
		Cpx::new(r_fs, x_fs).add(z_r)
	} else {
		let (dr, di) = z_hertzian_spectral_tuple(dx, dy, 0.0, ell, freq_scale, cfg);
		Cpx::new(dr, di).add(z_r)
	};
	z.tuple()
}

/// Isolated space-wave \((E_\theta,E_\phi)\) from the saddle of the slab spectral
/// \(\tilde E\). Zero for \(\theta>\pi/2\). No TM0 residue.
pub fn f_iso_slab_dipole(
	theta: f64,
	phi: f64,
	h: f64,
	ell: f64,
	freq_scale: f64,
	env: SlabEnv,
) -> (f64, f64, f64, f64) {
	if !(theta.is_finite()
		&& phi.is_finite()
		&& h.is_finite()
		&& ell.is_finite()
		&& freq_scale.is_finite()
		&& env.ok())
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
		return (0.0, 0.0, 0.0, 0.0);
	}
	let kx = k * st * cp;
	let ky = k * st * sp;
	let kz = k * ct;
	let k2 = k * k;
	let k_rho = (kx * kx + ky * ky).sqrt();
	let kz_a = Cpx::from_real(kz);
	let kz_d = kz_dielectric(k_rho, k, env.eps_c());
	let g_te = gamma_te(kz_a, kz_d, env.h_sub);
	let g_tm = gamma_tm(kz_a, kz_d, env.eps_c(), env.h_sub);
	let e_up = Cpx::new((kz * h).cos(), (kz * h).sin());
	let e_dn = Cpx::new((kz * h).cos(), -(kz * h).sin());
	let u_te = e_up.add(g_te.mul(e_dn));
	let u_tm = e_up.add(g_tm.mul(e_dn));
	let half = 0.5 * ETA0 * k * ell;
	let inv_kz = 1.0 / kz;
	let (ex_o, ey_o, ez_o) = if k_rho < 1e-14 {
		let u = u_te;
		(
			u.scale(-half * inv_kz),
			Cpx::zero(),
			Cpx::zero(),
		)
	} else {
		let inv_kr2 = 1.0 / (k_rho * k_rho);
		let w_te_x = ky * ky * inv_kr2;
		let w_tm_x = kx * kx * inv_kr2 * (kz * kz / k2);
		let w_te_y = -kx * ky * inv_kr2;
		let w_tm_y = kx * ky * inv_kr2 * (kz * kz / k2);
		let ex = u_te.scale(-half * inv_kz * w_te_x).add(u_tm.scale(-half * inv_kz * w_tm_x));
		let ey = u_te.scale(-half * inv_kz * w_te_y).add(u_tm.scale(-half * inv_kz * w_tm_y));
		// G_TM horizontal used (kz/k²) not (1/kz)*(kz²/k²) wait: -half/kz * (kz²/k²) = -half kz/k².
		// w_tm already includes kz²/k², times 1/kz. Good.
		let ez = u_tm.scale(half * kx / k2);
		(ex, ey, ez)
	};
	let jac = Cpx::new(0.0, -k * ct / (2.0 * PI));
	let ex = jac.mul(ex_o);
	let ey = jac.mul(ey_o);
	let ez = jac.mul(ez_o);
	let e_th = Cpx::new(
		ex.re * ct * cp + ey.re * ct * sp - ez.re * st,
		ex.im * ct * cp + ey.im * ct * sp - ez.im * st,
	);
	let e_ph = Cpx::new(-ex.re * sp + ey.re * cp, -ex.im * sp + ey.im * cp);
	(e_th.re, e_th.im, e_ph.re, e_ph.im)
}

pub fn f_iso_slab_dipole_power(
	theta: f64,
	phi: f64,
	h: f64,
	ell: f64,
	freq_scale: f64,
	env: SlabEnv,
) -> f64 {
	let (et_re, et_im, ep_re, ep_im) = f_iso_slab_dipole(theta, phi, h, ell, freq_scale, env);
	et_re * et_re + et_im * et_im + ep_re * ep_re + ep_im * ep_im
}

pub fn slab_dipole_power_budget(
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
	env: SlabEnv,
) -> SlabPowerBudget {
	slab_dipole_power_budget_cfg(h, ell, a, freq_scale, env, SpectralQuadConfig::DEFAULT)
}

pub fn slab_dipole_power_budget_cfg(
	h: f64,
	ell: f64,
	a: f64,
	freq_scale: f64,
	env: SlabEnv,
	cfg: SpectralQuadConfig,
) -> SlabPowerBudget {
	let nan = SlabPowerBudget {
		re_z_self: f64::NAN,
		p_rad: f64::NAN,
		p_sw: f64::NAN,
		p_diss: f64::NAN,
		closure_residual: f64::NAN,
	};
	if !geometry_ok(h, ell, a, freq_scale) || !env.ok() {
		return nan;
	}
	let split = integrate_refl(0.0, 0.0, h, ell, freq_scale, env, cfg);
	let (r_fs, _) = z_fs(ell, a, freq_scale);
	let z = Cpx::from_real(r_fs).add(split.total);
	let p_rad = r_fs + split.prop.re;
	let p_sw = split.sw_jump.re;
	let p_diss = z.re - p_rad - p_sw;
	SlabPowerBudget {
		re_z_self: z.re,
		p_rad,
		p_sw,
		p_diss,
		closure_residual: z.re - p_rad - p_sw - p_diss,
	}
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

	#[test]
	fn tm0_pole_in_bound_window() {
		let k = k_of(1.0);
		let env = SlabEnv::DEFAULT;
		let kp = tm0_pole(k, env).expect("TM0 pole");
		assert!(kp.re > k && kp.re < k * env.eps_r.sqrt(), "k_sw={}", kp.re);
		let d = tm_den(kp, k, env.eps_c(), env.h_sub);
		assert!(d.abs() < 1e-6 * k, "|D(k_sw)|={}", d.abs());
	}

	#[test]
	fn pec_limit_z_matches_closed() {
		let env = SlabEnv::PEC_LIMIT;
		let lags = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5), (0.3, 0.1)];
		for (dx, dy) in lags {
			let (cr, ci) = z_pair_pec_dipole(dx, dy, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
			let (sr, si) = z_pair_slab_dipole(
				dx,
				dy,
				DEFAULT_H,
				DEFAULT_ELL,
				DEFAULT_A,
				1.0,
				env,
			);
			assert!(sr.is_finite() && si.is_finite(), "slab NaN at ({dx},{dy})");
			let e = (sr - cr).hypot(si - ci);
			assert!(
				e < 1e-3,
				"PEC-limit |ΔZ|={e} at ({dx},{dy}) pec=({cr},{ci}) slab=({sr},{si})"
			);
		}
	}

	#[test]
	fn slab_self_positive_finite() {
		let (re, im) = z_pair_slab_dipole(
			0.0,
			0.0,
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			SlabEnv::DEFAULT,
		);
		assert!(re.is_finite() && im.is_finite());
		assert!(re > 0.0, "Re Z_self={re}");
	}

	#[test]
	fn slab_reciprocity() {
		let env = SlabEnv::DEFAULT;
		let (re, im) =
			z_pair_slab_dipole(0.4, -0.25, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		let (re_n, im_n) =
			z_pair_slab_dipole(-0.4, 0.25, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		close(re, re_n, 1e-9, "Z(Δ)=Z(-Δ) re");
		close(im, im_n, 1e-9, "Z(Δ)=Z(-Δ) im");
		let (re_xy, im_xy) =
			z_pair_slab_dipole(0.5, 0.1, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		let (re_yx, im_yx) =
			z_pair_slab_dipole(0.1, 0.5, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		let d = (re_xy - re_yx).hypot(im_xy - im_yx);
		assert!(d > 0.02, "x-dipole is not isotropic: |Δ|={d}");
	}

	#[test]
	fn pec_limit_pattern_matches() {
		let env = SlabEnv::PEC_LIMIT;
		for &(th, ph) in &[
			(0.0, 0.0),
			(PI * 0.25, 0.0),
			(PI * 0.25, 0.5 * PI),
			(0.4, 0.7),
		] {
			let (ctr, cti, cpr, cpi) = f_iso_pec_dipole(th, ph, DEFAULT_H, DEFAULT_ELL, 1.0);
			let (str, sti, spr, spi) =
				f_iso_slab_dipole(th, ph, DEFAULT_H, DEFAULT_ELL, 1.0, env);
			close(str, ctr, 5e-4, "Eθ re");
			close(sti, cti, 5e-4, "Eθ im");
			close(spr, cpr, 5e-4, "Eφ re");
			close(spi, cpi, 5e-4, "Eφ im");
			let pc = f_iso_pec_dipole_power(th, ph, DEFAULT_H, DEFAULT_ELL, 1.0);
			let ps = f_iso_slab_dipole_power(th, ph, DEFAULT_H, DEFAULT_ELL, 1.0, env);
			close(ps, pc, 0.05, "|F|²");
		}
		let (et_re, et_im, ep_re, ep_im) =
			f_iso_slab_dipole(PI * 0.5 + 0.2, 0.3, DEFAULT_H, DEFAULT_ELL, 1.0, env);
		assert_eq!((et_re, et_im, ep_re, ep_im), (0.0, 0.0, 0.0, 0.0));
	}

	#[test]
	fn slab_pattern_back_hemisphere_zero() {
		let p = f_iso_slab_dipole_power(2.0, 0.0, DEFAULT_H, DEFAULT_ELL, 1.0, SlabEnv::DEFAULT);
		assert_eq!(p, 0.0);
	}

	#[test]
	fn closure_lossless() {
		let b = slab_dipole_power_budget(
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			SlabEnv::DEFAULT,
		);
		assert!(b.re_z_self > 0.0 && b.p_rad.is_finite());
		assert!(b.p_sw > 0.0, "TM0 should carry power, P_sw={}", b.p_sw);
		let rel = b.p_diss.abs() / b.re_z_self.max(1e-12);
		assert!(
			rel < 1e-3,
			"lossless P_diss should vanish: P_diss={} ReZ={} Prad={} Psw={}",
			b.p_diss,
			b.re_z_self,
			b.p_rad,
			b.p_sw
		);
	}

	#[test]
	fn closure_lossy() {
		let lossless = SlabEnv::DEFAULT;
		let lossy = SlabEnv {
			eps_r: 10.0,
			h_sub: 0.05,
			tan_delta: 0.2,
		};
		let b0 = slab_dipole_power_budget(DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, lossless);
		let b = slab_dipole_power_budget(DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, lossy);
		assert!(
			b.p_diss > b0.p_diss + 1e-4 || b.p_diss > 0.0,
			"lossy P_diss={} lossless P_diss={}",
			b.p_diss,
			b0.p_diss
		);
		assert!(
			b.closure_residual.abs() <= 1e-12 * b.re_z_self.abs().max(1.0),
			"algebraic closure {}",
			b.closure_residual
		);
		assert!(
			(b.re_z_self - b0.re_z_self).abs() > 1e-3,
			"tanδ should change Re(Z): lossless={} lossy={}",
			b0.re_z_self,
			b.re_z_self
		);
	}

	#[test]
	fn pec_limit_budget_no_sw() {
		let b = slab_dipole_power_budget(
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			SlabEnv::PEC_LIMIT,
		);
		let (cr, _) = z_pair_pec_dipole(0.0, 0.0, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0);
		assert!(b.p_sw.abs() < 0.05 * cr.abs().max(1.0), "P_sw={}", b.p_sw);
		close(b.p_rad, cr, 0.05 * cr.abs().max(1.0) + 0.05, "P_rad vs Re Z_pec");
	}

	#[test]
	fn invalid_geometry_is_nan() {
		let (re, im) = z_pair_slab_dipole(
			0.0,
			0.0,
			0.0,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			SlabEnv::DEFAULT,
		);
		assert!(re.is_nan() && im.is_nan(), "h=0");
		let env = SlabEnv {
			eps_r: 10.0,
			h_sub: 0.0,
			tan_delta: 0.0,
		};
		let (re, im) =
			z_pair_slab_dipole(0.2, 0.0, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		assert!(re.is_nan() && im.is_nan(), "h_sub=0");
	}

	#[test]
	fn budget_quadrature_stable() {
		let env = SlabEnv::DEFAULT;
		let coarse = slab_dipole_power_budget_cfg(
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			env,
			SpectralQuadConfig::COARSE,
		);
		let fine = slab_dipole_power_budget_cfg(
			DEFAULT_H,
			DEFAULT_ELL,
			DEFAULT_A,
			1.0,
			env,
			SpectralQuadConfig::FINE,
		);
		let scale = fine.re_z_self.abs().max(1.0);
		assert!(
			(coarse.re_z_self - fine.re_z_self).abs() < 0.05 * scale,
			"ReZ coarse={} fine={}",
			coarse.re_z_self,
			fine.re_z_self
		);
		assert!((coarse.p_rad - fine.p_rad).abs() < 0.05 * scale);
		assert!((coarse.p_sw - fine.p_sw).abs() < 0.05 * scale);
	}

	#[ignore]
	#[test]
	fn bench_slab_per_lag() {
		use std::time::Instant;
		let env = SlabEnv::DEFAULT;
		let lags = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5)];
		for (dx, dy) in lags {
			let t0 = Instant::now();
			let mut z = (0.0, 0.0);
			let n = 4;
			for _ in 0..n {
				z = z_pair_slab_dipole(dx, dy, DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
			}
			let us = t0.elapsed().as_secs_f64() * 1e6 / n as f64;
			eprintln!("slab lag ({dx},{dy}): {us:.1} µs  Z=({:.4},{:.4})", z.0, z.1);
		}
		let t1 = Instant::now();
		let b = slab_dipole_power_budget(DEFAULT_H, DEFAULT_ELL, DEFAULT_A, 1.0, env);
		eprintln!(
			"budget {:.1} µs  ReZ={:.4} Prad={:.4} Psw={:.4} Pdiss={:.4}",
			t1.elapsed().as_secs_f64() * 1e6,
			b.re_z_self,
			b.p_rad,
			b.p_sw,
			b.p_diss
		);
	}
}
