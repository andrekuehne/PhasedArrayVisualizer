//! Radiated-power Gram \(P_H\) on a Gauss-μ × φ hemisphere quadrature.
//!
//! Isolated scalar fields from planar geometry (wavelengths) times a
//! power-conserving element directivity, then \(P_H = (P_0 / 4\pi) A A^H\).

use crate::bessel::j0f;
use crate::element::{PATTERN_COS_N, PATTERN_ISOTROPIC};
use crate::match_s::MatchedS;
use crate::quadrature::HemisphereQuad;
use crate::sincos::{load4, sincos_f32x4, store4};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use wide::f32x4;

pub const P0: f32 = 0.5;

const SAMPLE_BLOCK: usize = 64;

/// Deterministic hasher so the unique-ρ cache works on wasm32 (no getrandom).
#[derive(Default, Clone)]
struct U32Hasher(u64);

impl Hasher for U32Hasher {
	fn finish(&self) -> u64 {
		self.0
	}
	fn write(&mut self, bytes: &[u8]) {
		let mut h = self.0;
		for &b in bytes {
			h = h.wrapping_mul(0x100000001b3).wrapping_add(b as u64);
		}
		self.0 = h;
	}
	fn write_u32(&mut self, i: u32) {
		self.0 = i as u64;
	}
}

type RhoCache = HashMap<u32, f32, BuildHasherDefault<U32Hasher>>;

pub struct PradState {
	pub quad: Option<HemisphereQuad>,
	pub n: usize,
	pub m: usize,
	pub amp: Vec<f32>,
	pub a_re: Vec<f32>,
	pub a_im: Vec<f32>,
	pub p_re: Vec<f32>,
	pub p_im: Vec<f32>,
	pub n_unique_rho: usize,
	pub z0: Vec<f64>,
	pub r_re: Vec<f64>,
	pub r_im: Vec<f64>,
	pub s_re: Vec<f64>,
	pub s_im: Vec<f64>,
	pub match_iterations: u32,
	pub match_residual: f64,
}

impl PradState {
	pub fn new() -> Self {
		Self {
			quad: None,
			n: 0,
			m: 0,
			amp: Vec::new(),
			a_re: Vec::new(),
			a_im: Vec::new(),
			p_re: Vec::new(),
			p_im: Vec::new(),
			n_unique_rho: 0,
			z0: Vec::new(),
			r_re: Vec::new(),
			r_im: Vec::new(),
			s_re: Vec::new(),
			s_im: Vec::new(),
			match_iterations: 0,
			match_residual: 0.0,
		}
	}

	pub fn set_quadrature(&mut self, n_mu: u32, n_phi: u32) {
		let n_mu = n_mu.max(1) as usize;
		let n_phi = n_phi.max(1) as usize;
		self.quad = Some(HemisphereQuad::new(n_mu, n_phi));
		self.m = n_mu * n_phi;
		self.amp.clear();
		self.a_re.clear();
		self.a_im.clear();
	}

	pub fn n_samples(&self) -> usize {
		self.quad.as_ref().map(|q| q.n_samples()).unwrap_or(self.m)
	}

	pub fn n_elements(&self) -> usize {
		self.n
	}

	pub fn fill_isolated(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		self.fill_isolated_range(x, y, frequency_scale, element_kind, element_n, 0, 0);
	}

	/// Fill isolated steering columns for `sample0 .. sample0+sample_count`.
	/// `sample_count == 0` means the remainder of the quadrature grid.
	pub fn fill_isolated_range(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
		sample0: u32,
		sample_count: u32,
	) {
		if self.quad.is_none() || x.len() != y.len() || x.is_empty() {
			self.n = 0;
			return;
		}
		let m_full = self.quad.as_ref().unwrap().n_samples();
		let s0 = (sample0 as usize).min(m_full);
		let sc = if sample_count == 0 {
			m_full.saturating_sub(s0)
		} else {
			(sample_count as usize).min(m_full.saturating_sub(s0))
		};
		if sc == 0 {
			self.n = 0;
			self.m = 0;
			return;
		}
		self.n = x.len();
		self.m = sc;
		let nm = self.n * self.m;
		self.amp.resize(self.m, 0.0);
		self.a_re.resize(nm, 0.0);
		self.a_im.resize(nm, 0.0);
		let k = std::f32::consts::TAU * frequency_scale;
		let quad = self.quad.as_ref().unwrap();
		fill_amp(
			&quad.mu[s0..s0 + sc],
			&quad.omega[s0..s0 + sc],
			element_kind,
			element_n,
			&mut self.amp,
		);
		fill_steering(
			x,
			y,
			k,
			&quad.u[s0..s0 + sc],
			&quad.v[s0..s0 + sc],
			&self.amp,
			&mut self.a_re,
			&mut self.a_im,
		);
	}

	pub fn form_gram(&mut self) {
		if self.n == 0 || self.m == 0 || self.a_re.len() != self.n * self.m {
			self.p_re.clear();
			self.p_im.clear();
			return;
		}
		let nn = self.n * self.n;
		self.p_re.clear();
		self.p_re.resize(nn, 0.0);
		self.p_im.clear();
		self.p_im.resize(nn, 0.0);
		gram_hermitian(
			&self.a_re,
			&self.a_im,
			self.n,
			self.m,
			&mut self.p_re,
			&mut self.p_im,
		);
		let scale = P0 / (4.0 * std::f32::consts::PI);
		for v in self.p_re.iter_mut() {
			*v *= scale;
		}
		for v in self.p_im.iter_mut() {
			*v *= scale;
		}
		hermitianize(self.n, &mut self.p_re, &mut self.p_im);
	}

	pub fn compute(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		self.fill_isolated(x, y, frequency_scale, element_kind, element_n);
		self.form_gram();
	}

	/// Axisymmetric planar fast path: Gauss-μ of \(D(\mu) J_0(k\rho\sqrt{1-\mu^2})\).
	/// Uses `n_mu` from the last `set_quadrature`; does not fill `A`.
	pub fn compute_j0(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		if self.quad.is_none() || x.len() != y.len() || x.is_empty() {
			self.n = 0;
			self.n_unique_rho = 0;
			self.p_re.clear();
			self.p_im.clear();
			return;
		}
		let n = x.len();
		self.n = n;
		let nn = n * n;
		self.p_re.clear();
		self.p_re.resize(nn, 0.0);
		self.p_im.clear();
		self.p_im.resize(nn, 0.0);
		let k = std::f32::consts::TAU * frequency_scale;
		let quad = self.quad.as_ref().unwrap();
		let n_mu = quad.n_mu;
		let mut coeff = vec![0.0f64; n_mu];
		let mut s_mu = vec![0.0f32; n_mu];
		for i in 0..n_mu {
			let mu = quad.mu1d[i];
			let d = element_directivity(mu, element_kind, element_n);
			coeff[i] = (P0 as f64 * 0.5) * (quad.w_mu[i] as f64) * (d as f64);
			s_mu[i] = (1.0 - mu * mu).max(0.0).sqrt();
		}
		let mut cache = RhoCache::default();
		for p in 0..n {
			for q in 0..=p {
				let dx = x[p] - x[q];
				let dy = y[p] - y[q];
				let rho2 = dx.mul_add(dx, dy * dy);
				let val = *cache.entry(rho2.to_bits()).or_insert_with(|| {
					radial_integral(rho2, k, &coeff, &s_mu)
				});
				self.p_re[p * n + q] = val;
				self.p_re[q * n + p] = val;
			}
		}
		self.n_unique_rho = cache.len();
	}

	fn clear_matched(&mut self) {
		self.z0.clear();
		self.r_re.clear();
		self.r_im.clear();
		self.s_re.clear();
		self.s_im.clear();
		self.match_iterations = 0;
		self.match_residual = 0.0;
	}

	/// \(R = 2 Z_\mathrm{ref} P_H\), simultaneous real match, power-wave \(S\).
	/// Uses the current Gram; does not fill \(A\).
	pub fn form_matched_s(&mut self, z_ref: f32) {
		if self.n == 0 || self.p_re.len() != self.n * self.n {
			self.clear_matched();
			return;
		}
		let m = MatchedS::from_gram(&self.p_re, &self.p_im, self.n, z_ref as f64);
		self.z0 = m.z0;
		self.r_re = m.r_re;
		self.r_im = m.r_im;
		self.s_re = m.s_re;
		self.s_im = m.s_im;
		self.match_iterations = m.iterations;
		self.match_residual = m.residual;
	}
}

fn radial_integral(rho2: f32, k: f32, coeff: &[f64], s_mu: &[f32]) -> f32 {
	if !(rho2 > 0.0) {
		return coeff.iter().sum::<f64>() as f32;
	}
	let kr = k * rho2.sqrt();
	let mut acc = 0.0f64;
	for i in 0..coeff.len() {
		acc += coeff[i] * (j0f(kr * s_mu[i]) as f64);
	}
	acc as f32
}

fn element_directivity(mu: f32, kind: u32, n: f32) -> f32 {
	if mu < 0.0 {
		return 0.0;
	}
	match kind {
		PATTERN_COS_N => {
			let n = if n.is_finite() { n.max(0.0) } else { 0.0 };
			if n == 0.0 {
				2.0
			} else {
				2.0 * (n + 1.0) * mu.powf(n)
			}
		}
		PATTERN_ISOTROPIC | _ => 2.0,
	}
}

fn fill_amp(mu: &[f32], omega: &[f32], kind: u32, n: f32, amp: &mut [f32]) {
	let m = mu.len().min(omega.len()).min(amp.len());
	for s in 0..m {
		let d = element_directivity(mu[s], kind, n);
		let w = omega[s] * d;
		amp[s] = if w > 0.0 { w.sqrt() } else { 0.0 };
	}
}

fn fill_steering(
	x: &[f32],
	y: &[f32],
	k: f32,
	u: &[f32],
	v: &[f32],
	amp: &[f32],
	a_re: &mut [f32],
	a_im: &mut [f32],
) {
	let n = x.len();
	let m = u.len();
	for p in 0..n {
		let xp = x[p] * k;
		let yp = y[p] * k;
		let off = p * m;
		let are = &mut a_re[off..off + m];
		let aim = &mut a_im[off..off + m];
		let mut s = 0;
		let xp4 = f32x4::splat(xp);
		let yp4 = f32x4::splat(yp);
		while s + 4 <= m {
			let phase = xp4.mul_add(load4(u, s), yp4 * load4(v, s));
			let (sn, cs) = sincos_f32x4(phase);
			let g = load4(amp, s);
			store4(are, s, g * cs);
			store4(aim, s, g * (-sn));
			s += 4;
		}
		while s < m {
			let phase = xp.mul_add(u[s], yp * v[s]);
			let (sn, cs) = phase.sin_cos();
			are[s] = amp[s] * cs;
			aim[s] = -amp[s] * sn;
			s += 1;
		}
	}
}

fn gram_hermitian(
	a_re: &[f32],
	a_im: &[f32],
	n: usize,
	m: usize,
	p_re: &mut [f32],
	p_im: &mut [f32],
) {
	for s0 in (0..m).step_by(SAMPLE_BLOCK) {
		let s1 = (s0 + SAMPLE_BLOCK).min(m);
		for p in 0..n {
			let pr = &a_re[p * m + s0..p * m + s1];
			let pi = &a_im[p * m + s0..p * m + s1];
			for q in 0..=p {
				let qr = &a_re[q * m + s0..q * m + s1];
				let qi = &a_im[q * m + s0..q * m + s1];
				let (re, im) = dot_conj(pr, pi, qr, qi);
				p_re[p * n + q] += re;
				p_im[p * n + q] += im;
			}
		}
	}
	for p in 0..n {
		for q in 0..p {
			p_re[q * n + p] = p_re[p * n + q];
			p_im[q * n + p] = -p_im[p * n + q];
		}
	}
}

fn dot_conj(pr: &[f32], pi: &[f32], qr: &[f32], qi: &[f32]) -> (f32, f32) {
	let len = pr.len();
	let mut re = f32x4::splat(0.0);
	let mut im = f32x4::splat(0.0);
	let mut s = 0;
	while s + 4 <= len {
		let pr4 = load4(pr, s);
		let pi4 = load4(pi, s);
		let qr4 = load4(qr, s);
		let qi4 = load4(qi, s);
		re = pr4.mul_add(qr4, re);
		re = pi4.mul_add(qi4, re);
		im = pi4.mul_add(qr4, im);
		im = (-pr4).mul_add(qi4, im);
		s += 4;
	}
	let ra = re.to_array();
	let ia = im.to_array();
	let mut re_s = ra[0] + ra[1] + ra[2] + ra[3];
	let mut im_s = ia[0] + ia[1] + ia[2] + ia[3];
	while s < len {
		re_s += pr[s] * qr[s] + pi[s] * qi[s];
		im_s += pi[s] * qr[s] - pr[s] * qi[s];
		s += 1;
	}
	(re_s, im_s)
}

fn hermitianize(n: usize, p_re: &mut [f32], p_im: &mut [f32]) {
	for p in 0..n {
		for q in 0..=p {
			let i = p * n + q;
			let j = q * n + p;
			let re = 0.5 * (p_re[i] + p_re[j]);
			let im = 0.5 * (p_im[i] - p_im[j]);
			p_re[i] = re;
			p_re[j] = re;
			p_im[i] = im;
			p_im[j] = -im;
		}
		p_im[p * n + p] = 0.0;
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	fn close(a: f32, b: f32, tol: f32, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	#[test]
	fn n1_isotropic_radiates_p0() {
		let mut s = PradState::new();
		s.set_quadrature(8, 16);
		s.compute(&[0.0], &[0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		assert_eq!(s.n_elements(), 1);
		close(s.p_re[0], P0, 2e-4, "P11 re");
		close(s.p_im[0], 0.0, 1e-5, "P11 im");
	}

	#[test]
	fn n1_cos_n_radiates_p0() {
		let mut s = PradState::new();
		s.set_quadrature(12, 16);
		let n = 10f32.powf(0.5) * 0.5 - 1.0;
		s.compute(&[0.1], &[-0.2], 1.0, PATTERN_COS_N, n);
		close(s.p_re[0], P0, 5e-4, "P11 cos^n");
		close(s.p_im[0], 0.0, 1e-5, "P11 im");
	}

	#[test]
	fn n1_cos_n_zero_matches_isotropic() {
		let mut iso = PradState::new();
		let mut cos0 = PradState::new();
		iso.set_quadrature(8, 12);
		cos0.set_quadrature(8, 12);
		iso.compute(&[0.0], &[0.0], 1.0, PATTERN_ISOTROPIC, 0.8);
		cos0.compute(&[0.0], &[0.0], 1.0, PATTERN_COS_N, 0.0);
		close(iso.p_re[0], cos0.p_re[0], 1e-5, "n=0 vs isotropic");
	}

	#[test]
	fn two_coincident_elements() {
		let mut s = PradState::new();
		s.set_quadrature(8, 16);
		s.compute(&[0.3, 0.3], &[0.1, 0.1], 1.0, PATTERN_ISOTROPIC, 0.0);
		close(s.p_re[0], P0, 2e-4, "P11");
		close(s.p_re[3], P0, 2e-4, "P22");
		close(s.p_re[1], P0, 2e-4, "P12");
		close(s.p_re[2], P0, 2e-4, "P21");
		close(s.p_im[1], 0.0, 2e-4, "P12 im");
	}

	#[test]
	fn large_separation_offdiag_small() {
		let mut s = PradState::new();
		s.set_quadrature(48, 96);
		s.compute(&[0.0, 20.0], &[0.0, 0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		close(s.p_re[0], P0, 1e-3, "P11");
		let mag = (s.p_re[1] * s.p_re[1] + s.p_im[1] * s.p_im[1]).sqrt();
		assert!(mag < 0.05, "far |P12|={mag}");
	}

	#[test]
	fn hermitian_and_real_for_planar() {
		let mut s = PradState::new();
		s.set_quadrature(16, 32);
		s.compute(&[0.0, 0.5, 1.0], &[0.0, 0.25, -0.25], 1.2, PATTERN_COS_N, 1.0);
		let n = 3;
		for p in 0..n {
			for q in 0..n {
				let a = p * n + q;
				let b = q * n + p;
				close(s.p_re[a], s.p_re[b], 1e-5, "re Hermitian");
				close(s.p_im[a], -s.p_im[b], 1e-5, "im Hermitian");
				close(s.p_im[a], 0.0, 2e-4, "im ~ 0");
			}
			close(s.p_im[p * n + p], 0.0, 1e-8, "diag im");
		}
	}

	#[test]
	fn quadrature_convergence_short_baseline() {
		let x = [0.0f32, 0.5];
		let y = [0.0f32, 0.0];
		let mut coarse = PradState::new();
		let mut fine = PradState::new();
		coarse.set_quadrature(12, 24);
		fine.set_quadrature(24, 48);
		coarse.compute(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		fine.compute(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		close(fine.p_re[0], P0, 1e-4, "fine P11");
		let d11 = (coarse.p_re[0] - fine.p_re[0]).abs();
		let d12 = (coarse.p_re[1] - fine.p_re[1]).abs();
		assert!(d11 < 5e-4, "P11 Δ={d11}");
		assert!(d12 < 5e-3, "P12 Δ={d12}");
	}

	#[test]
	fn gram_matches_naive() {
		let mut s = PradState::new();
		s.set_quadrature(4, 6);
		s.fill_isolated(&[0.0, 0.4], &[0.1, -0.2], 0.9, PATTERN_ISOTROPIC, 0.0);
		let n = s.n;
		let m = s.m;
		let mut naive_re = vec![0.0f32; n * n];
		let mut naive_im = vec![0.0f32; n * n];
		for p in 0..n {
			for q in 0..n {
				let mut re = 0.0f32;
				let mut im = 0.0f32;
				for k in 0..m {
					let ar = s.a_re[p * m + k];
					let ai = s.a_im[p * m + k];
					let br = s.a_re[q * m + k];
					let bi = s.a_im[q * m + k];
					re += ar * br + ai * bi;
					im += ai * br - ar * bi;
				}
				naive_re[p * n + q] = re;
				naive_im[p * n + q] = im;
			}
		}
		s.form_gram();
		let scale = P0 / (4.0 * std::f32::consts::PI);
		for i in 0..n * n {
			close(s.p_re[i], naive_re[i] * scale, 2e-5, "gram re");
			close(s.p_im[i], naive_im[i] * scale, 2e-5, "gram im");
		}
	}

	#[test]
	fn sample_panels_sum_to_full_gram() {
		let x = [0.0f32, 0.4, 0.9, 1.3];
		let y = [0.1f32, -0.2, 0.05, 0.3];
		let mut full = PradState::new();
		full.set_quadrature(6, 10);
		full.compute(&x, &y, 1.15, PATTERN_COS_N, 0.8);
		let m = full.n_samples();
		let n = full.n_elements();
		let ranges = [(0u32, 7u32), (7, 11), (18, 0)];
		let mut acc_re = vec![0.0f32; n * n];
		let mut acc_im = vec![0.0f32; n * n];
		let mut covered = 0usize;
		let mut panel = PradState::new();
		panel.set_quadrature(6, 10);
		for (s0, sc) in ranges {
			panel.fill_isolated_range(&x, &y, 1.15, PATTERN_COS_N, 0.8, s0, sc);
			covered += panel.m;
			panel.form_gram();
			for i in 0..n * n {
				acc_re[i] += panel.p_re[i];
				acc_im[i] += panel.p_im[i];
			}
		}
		assert_eq!(covered, m);
		for i in 0..n * n {
			close(acc_re[i], full.p_re[i], 3e-5, "panel re");
			close(acc_im[i], full.p_im[i], 3e-5, "panel im");
		}
	}

	fn max_abs(a: &[f32], b: &[f32]) -> f32 {
		a.iter()
			.zip(b)
			.map(|(x, y)| (x - y).abs())
			.fold(0.0f32, f32::max)
	}

	fn rect_xy(nx: usize, ny: usize, dx: f32, dy: f32) -> (Vec<f32>, Vec<f32>) {
		let mut x = Vec::with_capacity(nx * ny);
		let mut y = Vec::with_capacity(nx * ny);
		for ix in 0..nx {
			for iy in 0..ny {
				x.push(dx * ix as f32);
				y.push(dy * iy as f32);
			}
		}
		(x, y)
	}

	#[test]
	fn j0_n1_isotropic_radiates_p0() {
		let mut s = PradState::new();
		s.set_quadrature(8, 16);
		s.compute_j0(&[0.0], &[0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		assert_eq!(s.n_elements(), 1);
		assert_eq!(s.n_unique_rho, 1);
		close(s.p_re[0], P0, 2e-6, "J0 P11 re");
		close(s.p_im[0], 0.0, 0.0, "J0 P11 im");
	}

	#[test]
	fn j0_n1_cos_n_radiates_p0() {
		let mut s = PradState::new();
		s.set_quadrature(12, 16);
		let n = 10f32.powf(0.5) * 0.5 - 1.0;
		s.compute_j0(&[0.1], &[-0.2], 1.0, PATTERN_COS_N, n);
		close(s.p_re[0], P0, 5e-4, "J0 P11 cos^n");
		close(s.p_im[0], 0.0, 0.0, "J0 P11 im");
	}

	#[test]
	fn j0_two_coincident_elements() {
		let mut s = PradState::new();
		s.set_quadrature(8, 16);
		s.compute_j0(&[0.3, 0.3], &[0.1, 0.1], 1.0, PATTERN_ISOTROPIC, 0.0);
		close(s.p_re[0], P0, 2e-6, "J0 P11");
		close(s.p_re[3], P0, 2e-6, "J0 P22");
		close(s.p_re[1], P0, 2e-6, "J0 P12");
		close(s.p_re[2], P0, 2e-6, "J0 P21");
		close(s.p_im[1], 0.0, 0.0, "J0 P12 im");
	}

	#[test]
	fn j0_matches_product_two_element() {
		let x = [0.0f32, 0.5];
		let y = [0.0f32, 0.0];
		let mut j0 = PradState::new();
		let mut prod = PradState::new();
		j0.set_quadrature(16, 2);
		prod.set_quadrature(16, 64);
		j0.compute_j0(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		prod.compute(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		close(j0.p_re[0], P0, 2e-6, "J0 P11");
		close(prod.p_re[0], P0, 1e-4, "prod P11");
		let dre = max_abs(&j0.p_re, &prod.p_re);
		let dim = j0.p_im.iter().fold(0.0f32, |m, v| m.max(v.abs()));
		assert!(dre < 1e-3, "J0 vs product max|Δre|={dre}");
		assert!(dim == 0.0, "J0 imag {dim}");
	}

	#[test]
	fn j0_matches_product_irregular_cos_n() {
		let x = [0.0f32, 0.4, 0.9, 1.3];
		let y = [0.1f32, -0.2, 0.05, 0.3];
		let mut j0 = PradState::new();
		let mut prod = PradState::new();
		j0.set_quadrature(24, 2);
		prod.set_quadrature(24, 96);
		j0.compute_j0(&x, &y, 1.15, PATTERN_COS_N, 0.8);
		prod.compute(&x, &y, 1.15, PATTERN_COS_N, 0.8);
		let n = 4;
		for p in 0..n {
			for q in 0..n {
				let a = p * n + q;
				let b = q * n + p;
				close(j0.p_re[a], j0.p_re[b], 0.0, "J0 re symmetric");
				close(j0.p_im[a], 0.0, 0.0, "J0 im");
			}
		}
		let dre = max_abs(&j0.p_re, &prod.p_re);
		assert!(dre < 1e-3, "J0 vs product max|Δre|={dre}");
	}

	#[test]
	fn j0_large_separation_offdiag_small() {
		let mut s = PradState::new();
		s.set_quadrature(48, 2);
		s.compute_j0(&[0.0, 20.0], &[0.0, 0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		close(s.p_re[0], P0, 1e-5, "J0 P11");
		assert!(s.p_re[1].abs() < 0.05, "J0 far |P12|={}", s.p_re[1]);
	}

	#[test]
	fn j0_rect_lattice_unique_rho_and_product() {
		let (x, y) = rect_xy(8, 8, 0.5, 0.5);
		let n = x.len();
		let mut j0 = PradState::new();
		let mut prod = PradState::new();
		j0.set_quadrature(16, 2);
		prod.set_quadrature(16, 64);
		j0.compute_j0(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		prod.compute(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		let pairs = n * (n + 1) / 2;
		assert!(
			j0.n_unique_rho < pairs / 8,
			"unique ρ={} vs pairs={}",
			j0.n_unique_rho,
			pairs
		);
		close(j0.p_re[0], P0, 2e-6, "J0 P11 lattice");
		let dre = max_abs(&j0.p_re, &prod.p_re);
		assert!(dre < 2e-3, "8x8 J0 vs product max|Δre|={dre}");
	}

	fn close64(a: f64, b: f64, tol: f64, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	#[test]
	fn j0_n1_matched_s_is_open() {
		use crate::match_s::{TAU, Z_REF};
		let mut s = PradState::new();
		s.set_quadrature(8, 2);
		s.compute_j0(&[0.0], &[0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32);
		close64(s.r_re[0], Z_REF, 1e-5, "R11");
		close64(s.z0[0], Z_REF, 1e-6, "z0");
		close64(s.s_re[0], 0.0, 1e-9, "S11 re");
		close64(s.s_im[0], 0.0, 1e-12, "S11 im");
		assert!(s.match_residual < TAU, "residual {}", s.match_residual);
	}

	#[test]
	fn j0_far_pair_weak_coupling() {
		use crate::match_s::{TAU, Z_REF};
		let mut s = PradState::new();
		s.set_quadrature(48, 2);
		s.compute_j0(&[0.0, 20.0], &[0.0, 0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32);
		close64(s.z0[0], Z_REF, 1.0, "z0_1 ~ 50");
		close64(s.z0[1], Z_REF, 1.0, "z0_2 ~ 50");
		assert!(s.match_residual < TAU, "residual {}", s.match_residual);
		let mag12 = (s.s_re[1] * s.s_re[1] + s.s_im[1] * s.s_im[1]).sqrt();
		assert!(mag12 < 0.1, "|S12|={mag12}");
		close64(s.s_re[0], 0.0, 2e-3, "S11");
		close64(s.s_re[3], 0.0, 2e-3, "S22");
	}

	#[test]
	fn j0_three_element_s_symmetric_and_matched() {
		use crate::match_s::{TAU, Z_REF};
		let mut s = PradState::new();
		s.set_quadrature(24, 2);
		s.compute_j0(&[0.0, 0.5, 1.0], &[0.0, 0.25, -0.25], 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32);
		let n = 3;
		assert!(s.match_residual < TAU, "residual {}", s.match_residual);
		for p in 0..n {
			let mag_ii = (s.s_re[p * n + p] * s.s_re[p * n + p]
				+ s.s_im[p * n + p] * s.s_im[p * n + p])
				.sqrt();
			assert!(mag_ii < 2e-3, "|S{p}{p}|={mag_ii}");
			for q in 0..n {
				close64(s.s_re[p * n + q], s.s_re[q * n + p], 1e-9, "S re symmetric");
				close64(s.s_im[p * n + q], s.s_im[q * n + p], 1e-9, "S im symmetric");
				close64(s.s_im[p * n + q], 0.0, 1e-8, "S im ~ 0");
			}
		}
	}
}
