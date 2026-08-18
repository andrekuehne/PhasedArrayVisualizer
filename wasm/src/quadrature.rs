//! Gauss–Legendre μ × trapezoid-φ product quadrature on the front hemisphere.
//!
//! μ = cos θ ∈ [0, 1], dΩ = dμ dφ (no sin θ). φ ∈ [0, 2π) with equal spacing.

const PI: f64 = std::f64::consts::PI;

/// n-point Gauss–Legendre nodes and weights on [-1, 1].
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
	assert!(n >= 1, "Gauss–Legendre order must be >= 1");
	let mut x = vec![0.0f64; n];
	let mut w = vec![0.0f64; n];
	let m = (n + 1) / 2;
	let nf = n as f64;
	for i in 0..m {
		let mut z = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
		let mut p1;
		let mut p2;
		let mut pp;
		loop {
			p1 = 1.0;
			p2 = 0.0;
			for j in 1..=n {
				let jf = j as f64;
				let p3 = p2;
				p2 = p1;
				p1 = ((2.0 * jf - 1.0) * z * p2 - (jf - 1.0) * p3) / jf;
			}
			pp = nf * (z * p1 - p2) / (z * z - 1.0);
			let z1 = z;
			z -= p1 / pp;
			if (z - z1).abs() < 1e-14 {
				break;
			}
		}
		x[i] = -z;
		x[n - 1 - i] = z;
		let wi = 2.0 / ((1.0 - z * z) * pp * pp);
		w[i] = wi;
		w[n - 1 - i] = wi;
	}
	(x, w)
}

/// Gauss–Legendre on the unit interval [0, 1]: μ = (ξ + 1)/2, w_μ = w_ξ / 2.
pub fn gauss_legendre_01(n: usize) -> (Vec<f64>, Vec<f64>) {
	let (xi, wi) = gauss_legendre(n);
	let mu: Vec<f64> = xi.iter().map(|x| 0.5 * (x + 1.0)).collect();
	let w: Vec<f64> = wi.iter().map(|ww| 0.5 * ww).collect();
	(mu, w)
}

/// Flattened product grid: μ outer, φ inner. Length `n_mu * n_phi`.
pub struct HemisphereQuad {
	#[allow(dead_code)]
	pub n_mu: usize,
	#[allow(dead_code)]
	pub n_phi: usize,
	pub mu: Vec<f32>,
	pub u: Vec<f32>,
	pub v: Vec<f32>,
	pub omega: Vec<f32>,
}

impl HemisphereQuad {
	pub fn new(n_mu: usize, n_phi: usize) -> Self {
		assert!(n_mu >= 1 && n_phi >= 1, "quadrature orders must be >= 1");
		let (mu64, wmu) = gauss_legendre_01(n_mu);
		let dphi = 2.0 * PI / n_phi as f64;
		let m = n_mu * n_phi;
		let mut mu = vec![0.0f32; m];
		let mut u = vec![0.0f32; m];
		let mut v = vec![0.0f32; m];
		let mut omega = vec![0.0f32; m];
		for i in 0..n_mu {
			let mui = mu64[i];
			let rho = (1.0 - mui * mui).max(0.0).sqrt();
			let wi = wmu[i] * dphi;
			for j in 0..n_phi {
				let phi = dphi * j as f64;
				let s = i * n_phi + j;
				mu[s] = mui as f32;
				u[s] = (rho * phi.cos()) as f32;
				v[s] = (rho * phi.sin()) as f32;
				omega[s] = wi as f32;
			}
		}
		Self {
			n_mu,
			n_phi,
			mu,
			u,
			v,
			omega,
		}
	}

	pub fn n_samples(&self) -> usize {
		self.mu.len()
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	fn close(a: f64, b: f64, tol: f64) {
		assert!(
			(a - b).abs() <= tol,
			"{a} vs {b} (tol {tol})"
		);
	}

	#[test]
	fn gl_n1() {
		let (x, w) = gauss_legendre(1);
		close(x[0], 0.0, 1e-14);
		close(w[0], 2.0, 1e-14);
	}

	#[test]
	fn gl_n2() {
		let (x, w) = gauss_legendre(2);
		let r = 1.0 / 3.0f64.sqrt();
		close(x[0], -r, 1e-12);
		close(x[1], r, 1e-12);
		close(w[0], 1.0, 1e-12);
		close(w[1], 1.0, 1e-12);
	}

	#[test]
	fn gl_n3() {
		let (x, w) = gauss_legendre(3);
		let r = (3.0 / 5.0f64).sqrt();
		close(x[0], -r, 1e-12);
		close(x[1], 0.0, 1e-12);
		close(x[2], r, 1e-12);
		close(w[0], 5.0 / 9.0, 1e-12);
		close(w[1], 8.0 / 9.0, 1e-12);
		close(w[2], 5.0 / 9.0, 1e-12);
	}

	#[test]
	fn gl_integrates_monomials() {
		let (x, w) = gauss_legendre(4);
		let mut i0 = 0.0;
		let mut i2 = 0.0;
		let mut i6 = 0.0;
		for k in 0..x.len() {
			i0 += w[k];
			i2 += w[k] * x[k] * x[k];
			i6 += w[k] * x[k].powi(6);
		}
		close(i0, 2.0, 1e-12);
		close(i2, 2.0 / 3.0, 1e-12);
		close(i6, 2.0 / 7.0, 1e-12);
	}

	#[test]
	fn interval_01_weights_sum_to_one() {
		let (_, w) = gauss_legendre_01(8);
		close(w.iter().sum(), 1.0, 1e-12);
	}

	#[test]
	fn hemisphere_solid_angle_is_two_pi() {
		let q = HemisphereQuad::new(8, 16);
		let omega: f64 = q.omega.iter().map(|w| *w as f64).sum();
		close(omega, 2.0 * PI, 1e-5);
		let four_pi: f64 = q.omega.iter().map(|w| 2.0 * *w as f64).sum();
		close(four_pi, 4.0 * PI, 2e-5);
	}

	#[test]
	fn sample_count_is_product() {
		let q = HemisphereQuad::new(5, 12);
		assert_eq!(q.n_samples(), 60);
		assert_eq!(q.n_mu, 5);
		assert_eq!(q.n_phi, 12);
	}

	#[test]
	fn directions_are_unit_hemisphere() {
		let q = HemisphereQuad::new(6, 8);
		for s in 0..q.n_samples() {
			let r2 = q.u[s] * q.u[s] + q.v[s] * q.v[s] + q.mu[s] * q.mu[s];
			assert!((r2 - 1.0).abs() < 2e-6, "r2={r2} s={s}");
			assert!(q.mu[s] >= -1e-6);
		}
	}
}
