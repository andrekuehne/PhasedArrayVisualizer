//! Radiated-power Gram \(P_H\) on a Gauss-μ × φ hemisphere quadrature.
//!
//! Isolated scalar fields from planar geometry (wavelengths) times a
//! power-conserving element directivity, then \(P_H = (P_0 / 4\pi) A A^H\).

use crate::bessel::j0f;
use crate::element::{PATTERN_COS_N, PATTERN_ISOTROPIC};
use crate::green::z_pair_pec_dipole;
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

	fn write_u64(&mut self, i: u64) {
		self.0 = self.0.wrapping_mul(0x100000001b3).wrapping_add(i);
	}
}

type RhoCache = HashMap<u32, f32, BuildHasherDefault<U32Hasher>>;
type LagCache = HashMap<(u64, u64), (f64, f64), BuildHasherDefault<U32Hasher>>;

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
	pub n_unique_lag: usize,
	pub z0: Vec<f64>,
	pub z0_im: Vec<f64>,
	pub r_re: Vec<f64>,
	pub r_im: Vec<f64>,
	pub s_re: Vec<f64>,
	pub s_im: Vec<f64>,
	pub t_re: Vec<f64>,
	pub t_im: Vec<f64>,
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
			n_unique_lag: 0,
			z0: Vec::new(),
			z0_im: Vec::new(),
			r_re: Vec::new(),
			r_im: Vec::new(),
			s_re: Vec::new(),
			s_im: Vec::new(),
			t_re: Vec::new(),
			t_im: Vec::new(),
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
		self.z0_im.clear();
		self.r_re.clear();
		self.r_im.clear();
		self.s_re.clear();
		self.s_im.clear();
		self.t_re.clear();
		self.t_im.clear();
		self.match_iterations = 0;
		self.match_residual = 0.0;
	}

	fn apply_matched(&mut self, m: MatchedS) {
		self.z0 = m.z0;
		self.z0_im = m.z0_im;
		self.r_re = m.r_re;
		self.r_im = m.r_im;
		self.s_re = m.s_re;
		self.s_im = m.s_im;
		self.t_re = m.t_re;
		self.t_im = m.t_im;
		self.match_iterations = m.iterations;
		self.match_residual = m.residual;
	}

	/// \(R = 2 Z_\mathrm{ref} P_H\), optional \(jX(\Delta x,\Delta y)+jX_\mathrm{self}I\),
	/// real \(z_0\), power-wave \(S\). Uses the current Gram; does not fill \(A\).
	/// Finite \(z_\mathrm{common,re}>0\) skips the per-port solver.
	pub fn form_matched_s(
		&mut self,
		z_ref: f32,
		x: &[f32],
		y: &[f32],
		x_nn: f32,
		alpha: f32,
		beta: f32,
		aniso: f32,
		z_common_re: f32,
		x_self: f32,
	) {
		if self.n == 0 || self.p_re.len() != self.n * self.n {
			self.clear_matched();
			return;
		}
		let m = MatchedS::from_gram_coupled(
			&self.p_re,
			&self.p_im,
			self.n,
			z_ref as f64,
			x,
			y,
			x_nn as f64,
			alpha as f64,
			beta as f64,
			aniso as f64,
			z_common_re as f64,
			x_self as f64,
		);
		self.apply_matched(m);
	}

	/// Like [`Self::form_matched_s`] with the propagation \(X\) overlay.
	/// Always a common real \(z_c\) (non-positive \(z_c\) becomes \(z_\mathrm{ref}\)).
	pub fn form_matched_s_propagation(
		&mut self,
		z_ref: f32,
		x: &[f32],
		y: &[f32],
		x_nn: f32,
		att: f32,
		eps_x: f32,
		eps_y: f32,
		freq: f32,
		z_common_re: f32,
		x_self: f32,
	) {
		if self.n == 0 || self.p_re.len() != self.n * self.n {
			self.clear_matched();
			return;
		}
		let m = MatchedS::from_gram_propagation(
			&self.p_re,
			&self.p_im,
			self.n,
			z_ref as f64,
			x,
			y,
			x_nn as f64,
			att as f64,
			eps_x as f64,
			eps_y as f64,
			freq as f64,
			z_common_re as f64,
			x_self as f64,
		);
		self.apply_matched(m);
	}

	/// Unique-lag PEC-dipole \(Z\) into `r_re`/`r_im`. No Gram, no `from_z`.
	pub fn fill_green_pec_dipole_z(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
	) {
		let h = h as f64;
		let ell = ell as f64;
		let a = a as f64;
		let fs = frequency_scale as f64;
		self.fill_green_pec_z_with(x, y, |dx, dy| z_pair_pec_dipole(dx, dy, h, ell, a, fs));
	}

	fn fill_green_pec_z_with<F>(&mut self, x: &[f32], y: &[f32], mut z_pair: F)
	where
		F: FnMut(f64, f64) -> (f64, f64),
	{
		if x.len() != y.len() || x.is_empty() {
			self.n = 0;
			self.n_unique_lag = 0;
			self.r_re.clear();
			self.r_im.clear();
			return;
		}
		let n = x.len();
		self.n = n;
		let nn = n * n;
		self.r_re.clear();
		self.r_re.resize(nn, 0.0);
		self.r_im.clear();
		self.r_im.resize(nn, 0.0);
		if let Some((nx, ny, ix, iy, xs, ys)) = uniform_product_lattice(x, y) {
			let mut tre = vec![0.0f64; nx * ny];
			let mut tim = vec![0.0f64; nx * ny];
			for di in 0..nx {
				for dj in 0..ny {
					let dx = xs[di] as f64 - xs[0] as f64;
					let dy = ys[dj] as f64 - ys[0] as f64;
					let (re, im) = z_pair(dx, dy);
					let k = di * ny + dj;
					tre[k] = re;
					tim[k] = im;
				}
			}
			for p in 0..n {
				for q in 0..=p {
					let di = ix[p].abs_diff(ix[q]);
					let dj = iy[p].abs_diff(iy[q]);
					let k = di * ny + dj;
					let re = tre[k];
					let im = tim[k];
					self.r_re[p * n + q] = re;
					self.r_im[p * n + q] = im;
					self.r_re[q * n + p] = re;
					self.r_im[q * n + p] = im;
				}
			}
			self.n_unique_lag = nx * ny;
		} else {
			let mut cache = LagCache::default();
			for p in 0..n {
				for q in 0..=p {
					let dx = (x[p] as f64 - x[q] as f64).abs();
					let dy = (y[p] as f64 - y[q] as f64).abs();
					let key = (dx.to_bits(), dy.to_bits());
					let (re, im) = *cache.entry(key).or_insert_with(|| z_pair(dx, dy));
					self.r_re[p * n + q] = re;
					self.r_im[p * n + q] = im;
					self.r_re[q * n + p] = re;
					self.r_im[q * n + p] = im;
				}
			}
			self.n_unique_lag = cache.len();
		}
	}

	#[cfg(test)]
	fn fill_green_pec_dipole_z_spectral(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
	) {
		let h = h as f64;
		let ell = ell as f64;
		let a = a as f64;
		let fs = frequency_scale as f64;
		self.fill_green_pec_z_with(x, y, |dx, dy| {
			crate::green_spectral::z_pair_pec_dipole_spectral(dx, dy, h, ell, a, fs)
		});
	}

	#[cfg(test)]
	fn fill_green_slab_dipole_z(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
		env: crate::green_slab::SlabEnv,
	) {
		let h = h as f64;
		let ell = ell as f64;
		let a = a as f64;
		let fs = frequency_scale as f64;
		self.fill_green_pec_z_with(x, y, |dx, dy| {
			crate::green_slab::z_pair_slab_dipole(dx, dy, h, ell, a, fs, env)
		});
	}

	#[cfg(test)]
	fn form_green_slab_dipole(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
		env: crate::green_slab::SlabEnv,
		z_ref: f32,
		z_common_re: f32,
		x_self: f64,
	) {
		self.fill_green_slab_dipole_z(x, y, frequency_scale, h, ell, a, env);
		self.form_from_z(z_ref, z_common_re, x_self);
	}

	/// [`MatchedS::from_z`] on the current Green \(Z\) (`r_re`/`r_im`).
	/// \(x_\mathrm{self}\) is added to \(\mathrm{diag}(\Im Z)\).
	pub fn form_from_z(&mut self, z_ref: f32, z_common_re: f32, x_self: f64) {
		if self.n == 0 || self.r_re.len() != self.n * self.n {
			self.clear_matched();
			return;
		}
		let m = MatchedS::from_z(
			&self.r_re,
			&self.r_im,
			self.n,
			z_ref as f64,
			z_common_re as f64,
			x_self,
		);
		self.apply_matched(m);
	}

	/// Unique-lag PEC-dipole \(Z\), then [`MatchedS::from_z`]. No Gram.
	pub fn form_green_pec_dipole(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
		z_ref: f32,
		z_common_re: f32,
		x_self: f64,
	) {
		self.fill_green_pec_dipole_z(x, y, frequency_scale, h, ell, a);
		self.form_from_z(z_ref, z_common_re, x_self);
	}
}

#[cfg(test)]
fn fill_green_pec_dipole_z_naive(
	s: &mut PradState,
	x: &[f32],
	y: &[f32],
	frequency_scale: f32,
	h: f32,
	ell: f32,
	a: f32,
) {
	if x.len() != y.len() || x.is_empty() {
		s.n = 0;
		s.r_re.clear();
		s.r_im.clear();
		return;
	}
	let n = x.len();
	s.n = n;
	let nn = n * n;
	s.r_re.clear();
	s.r_re.resize(nn, 0.0);
	s.r_im.clear();
	s.r_im.resize(nn, 0.0);
	let h = h as f64;
	let ell = ell as f64;
	let a = a as f64;
	let fs = frequency_scale as f64;
	for p in 0..n {
		for q in 0..n {
			let dx = x[p] as f64 - x[q] as f64;
			let dy = y[p] as f64 - y[q] as f64;
			let (re, im) = z_pair_pec_dipole(dx, dy, h, ell, a, fs);
			s.r_re[p * n + q] = re;
			s.r_im[p * n + q] = im;
		}
	}
}

fn unique_sorted_bits(v: &[f32]) -> Vec<f32> {
	let mut bits: Vec<u32> = v.iter().map(|x| x.to_bits()).collect();
	bits.sort_unstable();
	bits.dedup();
	bits.into_iter().map(f32::from_bits).collect()
}

fn is_uniform_sorted(v: &[f32]) -> bool {
	if v.len() <= 2 {
		return true;
	}
	let d0 = v[1] - v[0];
	if !d0.is_finite() || d0 == 0.0 {
		return false;
	}
	let tol = 1e-5 * d0.abs().max(1e-12);
	for i in 2..v.len() {
		let d = v[i] - v[i - 1];
		if !d.is_finite() || (d - d0).abs() > tol {
			return false;
		}
	}
	true
}

fn uniform_product_lattice(
	x: &[f32],
	y: &[f32],
) -> Option<(usize, usize, Vec<usize>, Vec<usize>, Vec<f32>, Vec<f32>)> {
	let n = x.len();
	let xs = unique_sorted_bits(x);
	let ys = unique_sorted_bits(y);
	let nx = xs.len();
	let ny = ys.len();
	if nx.checked_mul(ny) != Some(n) {
		return None;
	}
	if !is_uniform_sorted(&xs) || !is_uniform_sorted(&ys) {
		return None;
	}
	let mut ix = vec![0usize; n];
	let mut iy = vec![0usize; n];
	let mut seen = vec![false; n];
	for p in 0..n {
		let i = xs
			.binary_search_by(|a| a.to_bits().cmp(&x[p].to_bits()))
			.ok()?;
		let j = ys
			.binary_search_by(|a| a.to_bits().cmp(&y[p].to_bits()))
			.ok()?;
		let slot = i * ny + j;
		if seen[slot] {
			return None;
		}
		seen[slot] = true;
		ix[p] = i;
		iy[p] = j;
	}
	Some((nx, ny, ix, iy, xs, ys))
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
		s.form_matched_s(Z_REF as f32, &[0.0], &[0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		close64(s.r_re[0], Z_REF, 1e-5, "R11");
		close64(s.z0[0], Z_REF, 1e-6, "z0");
		close64(s.s_re[0], 0.0, 1e-9, "S11 re");
		close64(s.s_im[0], 0.0, 1e-12, "S11 im");
		close64(s.t_re[0], 1.0, 1e-6, "T11");
		close64(s.t_im[0], 0.0, 1e-12, "T11 im");
		assert!(s.match_residual < TAU, "residual {}", s.match_residual);
	}

	#[test]
	fn j0_far_pair_weak_coupling() {
		use crate::match_s::{TAU, Z_REF};
		let mut s = PradState::new();
		s.set_quadrature(48, 2);
		s.compute_j0(&[0.0, 20.0], &[0.0, 0.0], 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32, &[0.0, 20.0], &[0.0, 0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
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
		s.form_matched_s(Z_REF as f32, &[0.0, 0.5, 1.0], &[0.0, 0.25, -0.25], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
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

	#[test]
	fn j0_irregular_reactance_is_real_z0() {
		use crate::match_s::{TAU, Z_REF};
		let mut s = PradState::new();
		s.set_quadrature(24, 2);
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.25, -0.25];
		s.compute_j0(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32, &x, &y, 10.0, 2.0, 0.0, 0.0, 0.0, 0.0);
		assert!(s.match_residual < TAU, "residual {}", s.match_residual);
		assert!(s.z0_im.iter().all(|v| v.abs() == 0.0), "z0 real");
		assert!(s.r_im[1].abs() > 1.0, "X01");
		close64(s.r_im[1], s.r_im[3], 1e-12, "X01 = X10");
		close64(s.r_im[0], 0.0, 0.0, "X00");
		let mag00 = (s.s_re[0] * s.s_re[0] + s.s_im[0] * s.s_im[0]).sqrt();
		assert!(mag00 > 1e-3, "|S00| leftover reactance {mag00}");
	}

	#[test]
	fn j0_common_z0_is_flat_and_mismatched() {
		use crate::match_s::Z_REF;
		let mut s = PradState::new();
		s.set_quadrature(24, 2);
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.25, -0.25];
		s.compute_j0(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s(Z_REF as f32, &x, &y, 0.0, 0.0, 0.0, 0.0, Z_REF as f32, 0.0);
		assert_eq!(s.match_iterations, 0);
		close64(s.z0[0], Z_REF, 0.0, "z0_0");
		close64(s.z0[1], Z_REF, 0.0, "z0_1");
		close64(s.z0[2], Z_REF, 0.0, "z0_2");
		let mut max_sii = 0.0f64;
		for p in 0..3 {
			let mag_ii = (s.s_re[p * 3 + p] * s.s_re[p * 3 + p]
				+ s.s_im[p * 3 + p] * s.s_im[p * 3 + p])
				.sqrt();
			max_sii = max_sii.max(mag_ii);
		}
		assert!(max_sii > 0.01, "common |Sii|={max_sii}");
	}

	#[test]
	fn j0_propagation_nn_sign_flip() {
		use crate::match_s::Z_REF;
		let mut s = PradState::new();
		s.set_quadrature(24, 2);
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.0, 0.0];
		s.compute_j0(&x, &y, 1.0, PATTERN_ISOTROPIC, 0.0);
		s.form_matched_s_propagation(
			Z_REF as f32, &x, &y, 10.0, 0.0, 1.0, 1.0, 1.0, Z_REF as f32, 0.0,
		);
		close64(s.r_im[1], 10.0, 1e-9, "X01 nn");
		close64(s.r_im[2], -10.0, 1e-9, "X02 flip");
		close64(s.r_im[1], s.r_im[3], 1e-12, "X01 = X10");
		close64(s.z0[0], Z_REF, 0.0, "z0");
		close64(s.r_re[0], 2.0 * Z_REF as f64 * s.p_re[0] as f64, 1e-6, "R00 gram");
	}

	#[test]
	fn green_n1_matches_z_pair_self_no_gram() {
		use crate::green::{z_pair_pec_dipole, DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::match_s::Z_REF;
		let h = DEFAULT_H as f32 as f64;
		let ell = DEFAULT_ELL as f32 as f64;
		let a = DEFAULT_A as f32 as f64;
		let mut s = PradState::new();
		let p_re_before = s.p_re.len();
		s.form_green_pec_dipole(
			&[0.0],
			&[0.0],
			1.0,
			h as f32,
			ell as f32,
			a as f32,
			Z_REF as f32,
			0.0,
			0.0,
		);
		assert_eq!(s.n_elements(), 1);
		assert_eq!(s.p_re.len(), p_re_before);
		let (z_re, z_im) = z_pair_pec_dipole(0.0, 0.0, h, ell, a, 1.0);
		close64(s.r_re[0], z_re, 1e-12, "Z11 re");
		close64(s.r_im[0], z_im, 1e-12, "Z11 im");
		assert!(z_re > 0.0 && z_im.is_finite());
		let zc = z_re;
		s.form_green_pec_dipole(
			&[0.0],
			&[0.0],
			1.0,
			h as f32,
			ell as f32,
			a as f32,
			Z_REF as f32,
			zc as f32,
			0.0,
		);
		close64(s.z0[0], zc as f32 as f64, 1e-6, "z0 = Re Z11");
		assert_eq!(s.s_re.len(), 1);
		assert_eq!(s.t_re.len(), 1);
		assert!(s.s_re[0].is_finite() && s.s_im[0].is_finite());
		assert!(s.t_re[0].is_finite() && s.t_im[0].is_finite());
		let mag_s = s.s_re[0].hypot(s.s_im[0]);
		assert!(mag_s > 1e-3, "|S11| leftover X {mag_s}");
	}

	#[test]
	fn green_n1_cancelled_xself_is_open() {
		use crate::green::{z_pair_pec_dipole, DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::match_s::Z_REF;
		let h = DEFAULT_H as f32 as f64;
		let ell = DEFAULT_ELL as f32 as f64;
		let a = DEFAULT_A as f32 as f64;
		let (z_re, z_im) = z_pair_pec_dipole(0.0, 0.0, h, ell, a, 1.0);
		let mut s = PradState::new();
		s.form_green_pec_dipole(
			&[0.0],
			&[0.0],
			1.0,
			h as f32,
			ell as f32,
			a as f32,
			Z_REF as f32,
			z_re as f32,
			-z_im,
		);
		close64(s.r_im[0], 0.0, 1e-12, "X11 cancelled");
		close64(s.z0[0], z_re as f32 as f64, 1e-6, "z0 = Re Z11");
		let mag_s = s.s_re[0].hypot(s.s_im[0]);
		assert!(mag_s < 1e-6, "|S11| after cancel {mag_s}");
		let zc = s.z0[0];
		let t_want = (Z_REF / zc).sqrt();
		close64(s.t_re[0], t_want, 1e-6, "T11 Kurokawa");
		close64(s.t_im[0], 0.0, 1e-6, "T11 im");
	}

	#[test]
	fn green_n2_reciprocal_matches_z_pair() {
		use crate::green::{z_pair_pec_dipole, DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::match_s::Z_REF;
		let h = DEFAULT_H as f32 as f64;
		let ell = DEFAULT_ELL as f32 as f64;
		let a = DEFAULT_A as f32 as f64;
		let x = [0.0f32, 0.5];
		let y = [0.0f32, 0.0];
		let mut s = PradState::new();
		s.form_green_pec_dipole(
			&x,
			&y,
			1.0,
			h as f32,
			ell as f32,
			a as f32,
			Z_REF as f32,
			Z_REF as f32,
			0.0,
		);
		assert_eq!(s.n_elements(), 2);
		assert_eq!(s.s_re.len(), 4);
		assert_eq!(s.s_im.len(), 4);
		assert_eq!(s.t_re.len(), 4);
		assert_eq!(s.t_im.len(), 4);
		close64(s.r_re[1], s.r_re[2], 1e-12, "Z12 re = Z21 re");
		close64(s.r_im[1], s.r_im[2], 1e-12, "Z12 im = Z21 im");
		close64(s.s_re[1], s.s_re[2], 1e-12, "S12 re = S21 re");
		close64(s.s_im[1], s.s_im[2], 1e-12, "S12 im = S21 im");
		let (z12_re, z12_im) = z_pair_pec_dipole(
			x[0] as f64 - x[1] as f64,
			y[0] as f64 - y[1] as f64,
			h,
			ell,
			a,
			1.0,
		);
		close64(s.r_re[1], z12_re, 1e-12, "Z12 vs z_pair");
		close64(s.r_im[1], z12_im, 1e-12, "Z12 im vs z_pair");
		assert!(s.p_re.is_empty(), "Green path must not fill Gram");
	}

	#[test]
	fn green_4x4_unique_lag_matches_naive() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let (x, y) = rect_xy(4, 4, 0.5, 0.5);
		let mut uniq = PradState::new();
		let mut naive = PradState::new();
		uniq.fill_green_pec_dipole_z(&x, &y, 1.0, h, ell, a);
		fill_green_pec_dipole_z_naive(&mut naive, &x, &y, 1.0, h, ell, a);
		assert_eq!(uniq.n_unique_lag, 16, "4×4 lattice U = nx ny");
		assert_eq!(uniq.r_re.len(), 16 * 16);
		let mut max_d = 0.0f64;
		for i in 0..uniq.r_re.len() {
			max_d = max_d.max((uniq.r_re[i] - naive.r_re[i]).abs());
			max_d = max_d.max((uniq.r_im[i] - naive.r_im[i]).abs());
		}
		assert!(max_d < 1e-12, "unique vs naïve max|Δ|={max_d}");
	}

	#[test]
	fn green_irregular_triple_matches_naive() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let x = [0.0f32, 0.37, 0.9];
		let y = [0.0f32, 0.21, -0.4];
		let mut uniq = PradState::new();
		let mut naive = PradState::new();
		uniq.fill_green_pec_dipole_z(&x, &y, 1.0, h, ell, a);
		fill_green_pec_dipole_z_naive(&mut naive, &x, &y, 1.0, h, ell, a);
		assert!(
			uniq.n_unique_lag < 9,
			"unique lags {} should collapse vs N²",
			uniq.n_unique_lag
		);
		let mut max_d = 0.0f64;
		for i in 0..9 {
			max_d = max_d.max((uniq.r_re[i] - naive.r_re[i]).abs());
			max_d = max_d.max((uniq.r_im[i] - naive.r_im[i]).abs());
		}
		assert!(max_d < 1e-12, "irregular unique vs naïve max|Δ|={max_d}");
	}

	fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
		a.iter()
			.zip(b.iter())
			.map(|(x, y)| (x - y).abs())
			.fold(0.0f64, f64::max)
	}

	#[test]
	fn spectral_lag_fill_matches_closed_4x4() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let (x, y) = rect_xy(4, 4, 0.5, 0.5);
		let mut closed = PradState::new();
		let mut spectral = PradState::new();
		closed.fill_green_pec_dipole_z(&x, &y, 1.0, h, ell, a);
		spectral.fill_green_pec_dipole_z_spectral(&x, &y, 1.0, h, ell, a);
		assert_eq!(closed.n_unique_lag, 16);
		assert_eq!(spectral.n_unique_lag, 16);
		let mut max_d = 0.0f64;
		for i in 0..closed.r_re.len() {
			max_d = max_d.max((closed.r_re[i] - spectral.r_re[i]).abs());
			max_d = max_d.max((closed.r_im[i] - spectral.r_im[i]).abs());
		}
		assert!(max_d < 1e-3, "spectral vs closed unique-lag max|ΔZ|={max_d}");
	}

	#[test]
	fn spectral_z_from_z_agrees_closed_4x4() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::match_s::Z_REF;
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let (x, y) = rect_xy(4, 4, 0.5, 0.5);
		let mut closed = PradState::new();
		let mut spectral = PradState::new();
		closed.fill_green_pec_dipole_z(&x, &y, 1.0, h, ell, a);
		spectral.fill_green_pec_dipole_z_spectral(&x, &y, 1.0, h, ell, a);
		closed.form_from_z(Z_REF as f32, Z_REF as f32, 0.0);
		spectral.form_from_z(Z_REF as f32, Z_REF as f32, 0.0);
		let ds = max_abs_diff(&closed.s_re, &spectral.s_re)
			.max(max_abs_diff(&closed.s_im, &spectral.s_im));
		let dt = max_abs_diff(&closed.t_re, &spectral.t_re)
			.max(max_abs_diff(&closed.t_im, &spectral.t_im));
		assert!(ds < 1e-3, "from_z S spectral vs closed max|Δ|={ds}");
		assert!(dt < 1e-3, "from_z T spectral vs closed max|Δ|={dt}");
		assert!(
			closed.p_re.is_empty() && spectral.p_re.is_empty(),
			"Green path must not fill Gram"
		);
	}

	fn fill_green_slab_dipole_z_naive(
		s: &mut PradState,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		h: f32,
		ell: f32,
		a: f32,
		env: crate::green_slab::SlabEnv,
	) {
		if x.len() != y.len() || x.is_empty() {
			s.n = 0;
			s.r_re.clear();
			s.r_im.clear();
			return;
		}
		let n = x.len();
		s.n = n;
		let nn = n * n;
		s.r_re.clear();
		s.r_re.resize(nn, 0.0);
		s.r_im.clear();
		s.r_im.resize(nn, 0.0);
		let h = h as f64;
		let ell = ell as f64;
		let a = a as f64;
		let fs = frequency_scale as f64;
		for p in 0..n {
			for q in 0..n {
				let dx = x[p] as f64 - x[q] as f64;
				let dy = y[p] as f64 - y[q] as f64;
				let (re, im) = crate::green_slab::z_pair_slab_dipole(dx, dy, h, ell, a, fs, env);
				s.r_re[p * n + q] = re;
				s.r_im[p * n + q] = im;
			}
		}
	}

	#[test]
	fn slab_4x4_unique_lag_matches_naive() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::green_slab::SlabEnv;
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let env = SlabEnv::DEFAULT;
		let (x, y) = rect_xy(4, 4, 0.5, 0.5);
		let mut uniq = PradState::new();
		let mut naive = PradState::new();
		uniq.fill_green_slab_dipole_z(&x, &y, 1.0, h, ell, a, env);
		fill_green_slab_dipole_z_naive(&mut naive, &x, &y, 1.0, h, ell, a, env);
		assert_eq!(uniq.n_unique_lag, 16, "4×4 lattice U = nx ny");
		let mut max_d = 0.0f64;
		for i in 0..uniq.r_re.len() {
			max_d = max_d.max((uniq.r_re[i] - naive.r_re[i]).abs());
			max_d = max_d.max((uniq.r_im[i] - naive.r_im[i]).abs());
		}
		assert!(max_d < 1e-10, "slab unique vs naïve max|Δ|={max_d}");
		let n = 16;
		for p in 0..n {
			for q in 0..n {
				let a = uniq.r_re[p * n + q].hypot(uniq.r_im[p * n + q]);
				let b = uniq.r_re[q * n + p].hypot(uniq.r_im[q * n + p]);
				assert!((uniq.r_re[p * n + q] - uniq.r_re[q * n + p]).abs() < 1e-12);
				assert!((uniq.r_im[p * n + q] - uniq.r_im[q * n + p]).abs() < 1e-12);
				assert!(a.is_finite() && b.is_finite());
			}
		}
	}

	#[test]
	fn slab_4x4_from_z_s_finite_no_gram() {
		use crate::green::{DEFAULT_A, DEFAULT_ELL, DEFAULT_H};
		use crate::green_slab::SlabEnv;
		use crate::match_s::Z_REF;
		let h = DEFAULT_H as f32;
		let ell = DEFAULT_ELL as f32;
		let a = DEFAULT_A as f32;
		let (x, y) = rect_xy(4, 4, 0.5, 0.5);
		let mut s = PradState::new();
		s.form_green_slab_dipole(
			&x,
			&y,
			1.0,
			h,
			ell,
			a,
			SlabEnv::DEFAULT,
			Z_REF as f32,
			Z_REF as f32,
			0.0,
		);
		assert_eq!(s.s_re.len(), 16 * 16);
		assert!(s.s_re.iter().chain(s.s_im.iter()).all(|v| v.is_finite()));
		assert!(s.t_re.iter().chain(s.t_im.iter()).all(|v| v.is_finite()));
		assert!(s.p_re.is_empty(), "slab Green path must not fill Gram");
		for i in 0..16 {
			assert!(s.s_re[i * 16 + i].hypot(s.s_im[i * 16 + i]).is_finite());
		}
	}
}
