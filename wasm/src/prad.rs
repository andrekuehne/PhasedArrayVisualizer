//! Radiated-power Gram \(P_H\) on a Gauss-μ × φ hemisphere quadrature.
//!
//! Isolated scalar fields from planar geometry (wavelengths) times a
//! power-conserving element directivity, then \(P_H = (P_0 / 4\pi) A A^H\).

use crate::element::{PATTERN_COS_N, PATTERN_ISOTROPIC};
use crate::quadrature::HemisphereQuad;
use crate::sincos::{load4, sincos_f32x4, store4};
use wide::f32x4;

pub const P0: f32 = 0.5;

const SAMPLE_BLOCK: usize = 64;

pub struct PradState {
	pub quad: Option<HemisphereQuad>,
	pub n: usize,
	pub m: usize,
	pub amp: Vec<f32>,
	pub a_re: Vec<f32>,
	pub a_im: Vec<f32>,
	pub p_re: Vec<f32>,
	pub p_im: Vec<f32>,
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
}
