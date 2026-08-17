use crate::sincos::accumulate_linear;

pub const DOMAIN_SPHERICAL: u32 = 0;
pub const DOMAIN_UV: u32 = 1;
pub const DOMAIN_LUDWIG3: u32 = 2;

pub struct FarfieldState {
	pub n1: usize,
	pub n2: usize,
	pub re: Vec<f32>,
	pub im: Vec<f32>,
	pub total: Vec<f32>,
	pub x: Vec<f32>,
	pub y: Vec<f32>,
	pub mag: Vec<f32>,
	pub pha: Vec<f32>,
	pub ax1: Vec<f32>,
	pub ax2: Vec<f32>,
	scratch: Vec<f32>,
}

impl FarfieldState {
	pub fn new() -> Self {
		Self {
			n1: 0,
			n2: 0,
			re: Vec::new(),
			im: Vec::new(),
			total: Vec::new(),
			x: Vec::new(),
			y: Vec::new(),
			mag: Vec::new(),
			pha: Vec::new(),
			ax1: Vec::new(),
			ax2: Vec::new(),
			scratch: Vec::new(),
		}
	}

	pub fn prepare(&mut self, n1: u32, n2: u32) {
		self.n1 = n1 as usize;
		self.n2 = n2 as usize;
		let n = self.n1.saturating_mul(self.n2);
		self.re.clear();
		self.re.resize(n, 0.0);
		self.im.clear();
		self.im.resize(n, 0.0);
		self.total.clear();
		self.total.resize(n, 0.0);
	}

	pub fn set_inputs(
		&mut self,
		x: &[f32],
		y: &[f32],
		mag: &[f32],
		pha: &[f32],
		ax1: &[f32],
		ax2: &[f32],
	) {
		self.x.clear();
		self.x.extend_from_slice(x);
		self.y.clear();
		self.y.extend_from_slice(y);
		self.mag.clear();
		self.mag.extend_from_slice(mag);
		self.pha.clear();
		self.pha.extend_from_slice(pha);
		self.ax1.clear();
		self.ax1.extend_from_slice(ax1);
		self.ax2.clear();
		self.ax2.extend_from_slice(ax2);
	}

	fn inputs_ok(&self) -> bool {
		let n_el = self.x.len();
		n_el > 0
			&& n_el == self.y.len()
			&& n_el == self.mag.len()
			&& n_el == self.pha.len()
			&& self.n1 > 0
			&& self.n2 > 0
			&& self.ax1.len() == self.n1
			&& self.ax2.len() == self.n2
			&& self.re.len() == self.n1 * self.n2
	}

	pub fn accumulate_tile(&mut self, domain: u32, frequency_scale: f32, row0: u32, row_count: u32) {
		if !self.inputs_ok() {
			return;
		}
		let row0 = row0 as usize;
		let row_end = (row0 + row_count as usize).min(self.n2);
		if row0 >= row_end {
			return;
		}
		match domain {
			DOMAIN_SPHERICAL => self.accumulate_spherical(frequency_scale, row0, row_end),
			DOMAIN_UV => self.accumulate_uv(row0, row_end),
			DOMAIN_LUDWIG3 => self.accumulate_ludwig3(row0, row_end),
			_ => {}
		}
	}

	fn accumulate_spherical(&mut self, frequency_scale: f32, row0: usize, row_end: usize) {
		let n1 = self.n1;
		let sc = std::f32::consts::TAU * frequency_scale;
		self.scratch.clear();
		self.scratch.resize(n1, 0.0);
		for (i, &theta) in self.ax1.iter().enumerate() {
			self.scratch[i] = sc * theta.sin();
		}
		let n_el = self.x.len();
		for row in row0..row_end {
			let phi = self.ax2[row];
			let cphi = phi.cos();
			let sphi = phi.sin();
			let off = row * n1;
			for i in 0..n_el {
				let scale = self.x[i].mul_add(cphi, self.y[i] * sphi);
				let (re_row, im_row) = split_row(&mut self.re, &mut self.im, off, n1);
				accumulate_linear(
					n1,
					self.mag[i],
					scale,
					self.pha[i],
					&self.scratch,
					re_row,
					im_row,
				);
			}
		}
	}

	fn accumulate_uv(&mut self, row0: usize, row_end: usize) {
		let n1 = self.n1;
		let two_pi = std::f32::consts::TAU;
		let n_el = self.x.len();
		for row in row0..row_end {
			let v = self.ax2[row];
			let off = row * n1;
			for i in 0..n_el {
				let a = self.x[i] * two_pi;
				let b = self.y[i].mul_add(v * two_pi, self.pha[i]);
				let (re_row, im_row) = split_row(&mut self.re, &mut self.im, off, n1);
				accumulate_linear(n1, self.mag[i], a, b, &self.ax1, re_row, im_row);
			}
		}
	}

	fn accumulate_ludwig3(&mut self, row0: usize, row_end: usize) {
		let n1 = self.n1;
		let two_pi = std::f32::consts::TAU;
		self.scratch.clear();
		self.scratch.resize(n1, 0.0);
		for (i, &az) in self.ax1.iter().enumerate() {
			self.scratch[i] = az.sin();
		}
		let n_el = self.x.len();
		for row in row0..row_end {
			let el = self.ax2[row];
			let cel = el.cos();
			let sel = el.sin();
			let off = row * n1;
			for i in 0..n_el {
				let xxv = self.x[i] * cel;
				let yyv = self.y[i] * sel;
				let a = xxv * two_pi;
				let b = yyv.mul_add(two_pi, self.pha[i]);
				let (re_row, im_row) = split_row(&mut self.re, &mut self.im, off, n1);
				accumulate_linear(n1, self.mag[i], a, b, &self.scratch, re_row, im_row);
			}
		}
	}

	pub fn finalize(&mut self, n_elements: u32) -> f32 {
		let n = n_elements.max(1) as f32;
		let mut max_value = f32::NEG_INFINITY;
		for i in 0..self.total.len() {
			let t = (self.re[i] * self.re[i] + self.im[i] * self.im[i]).abs() / n;
			self.total[i] = t;
			if t > max_value {
				max_value = t;
			}
		}
		max_value
	}
}

fn split_row<'a>(
	re: &'a mut [f32],
	im: &'a mut [f32],
	off: usize,
	n1: usize,
) -> (&'a mut [f32], &'a mut [f32]) {
	(&mut re[off..off + n1], &mut im[off..off + n1])
}

#[cfg(test)]
mod tests {
	use super::*;

	fn js_spherical(
		x: &[f32],
		y: &[f32],
		mag: &[f32],
		pha: &[f32],
		theta: &[f32],
		phi: &[f32],
		freq: f32,
	) -> (Vec<f32>, Vec<f32>) {
		let n1 = theta.len();
		let n2 = phi.len();
		let mut re = vec![0.0f32; n1 * n2];
		let mut im = vec![0.0f32; n1 * n2];
		let sc = std::f32::consts::TAU * freq;
		let jk: Vec<f32> = theta.iter().map(|t| sc * t.sin()).collect();
		for i in 0..x.len() {
			for ip in 0..n2 {
				let xxv = x[i] * phi[ip].cos();
				let yyv = y[i] * phi[ip].sin();
				for it in 0..n1 {
					let v = (xxv + yyv) * jk[it] + pha[i];
					let idx = ip * n1 + it;
					re[idx] += mag[i] * v.cos();
					im[idx] += mag[i] * v.sin();
				}
			}
		}
		(re, im)
	}

	#[test]
	fn spherical_matches_reference() {
		let x = [0.0f32, 0.5, 1.0];
		let y = [0.0f32, 0.25, -0.25];
		let mag = [1.0f32, 0.7, 0.4];
		let pha = [0.0f32, 0.3, -1.2];
		let theta: Vec<f32> = (0..7)
			.map(|i| -std::f32::consts::FRAC_PI_2 + i as f32 * std::f32::consts::PI / 6.0)
			.collect();
		let phi: Vec<f32> = (0..5)
			.map(|i| -std::f32::consts::FRAC_PI_2 + i as f32 * std::f32::consts::PI / 4.0)
			.collect();
		let (re_ref, im_ref) = js_spherical(&x, &y, &mag, &pha, &theta, &phi, 1.1);

		let mut s = FarfieldState::new();
		s.prepare(theta.len() as u32, phi.len() as u32);
		s.set_inputs(&x, &y, &mag, &pha, &theta, &phi);
		s.accumulate_tile(DOMAIN_SPHERICAL, 1.1, 0, phi.len() as u32);

		for i in 0..re_ref.len() {
			assert!(
				(s.re[i] - re_ref[i]).abs() < 2e-4,
				"re[{i}] {} vs {}",
				s.re[i],
				re_ref[i]
			);
			assert!(
				(s.im[i] - im_ref[i]).abs() < 2e-4,
				"im[{i}] {} vs {}",
				s.im[i],
				im_ref[i]
			);
		}
	}

	#[test]
	fn uv_matches_reference() {
		let x = [0.0f32, 0.5];
		let y = [0.1f32, -0.2];
		let mag = [1.0f32, 0.5];
		let pha = [0.2f32, -0.4];
		let u = [-1.0f32, -0.5, 0.0, 0.5, 1.0];
		let v = [-1.0f32, 0.0, 1.0];
		let two_pi = std::f32::consts::TAU;
		let mut re_ref = vec![0.0f32; u.len() * v.len()];
		let mut im_ref = vec![0.0f32; u.len() * v.len()];
		for i in 0..x.len() {
			for iv in 0..v.len() {
				let xxv = x[i];
				let yyv = y[i] * v[iv];
				for iu in 0..u.len() {
					let phase = (xxv * u[iu] + yyv) * two_pi + pha[i];
					let idx = iv * u.len() + iu;
					re_ref[idx] += mag[i] * phase.cos();
					im_ref[idx] += mag[i] * phase.sin();
				}
			}
		}
		let mut s = FarfieldState::new();
		s.prepare(u.len() as u32, v.len() as u32);
		s.set_inputs(&x, &y, &mag, &pha, &u, &v);
		s.accumulate_tile(DOMAIN_UV, 1.0, 0, v.len() as u32);
		for i in 0..re_ref.len() {
			assert!((s.re[i] - re_ref[i]).abs() < 2e-4);
			assert!((s.im[i] - im_ref[i]).abs() < 2e-4);
		}
	}

	#[test]
	fn finalize_peak_and_norm() {
		let mut s = FarfieldState::new();
		s.prepare(2, 1);
		s.re = vec![3.0, 0.0];
		s.im = vec![4.0, 0.0];
		s.total = vec![0.0, 0.0];
		let max = s.finalize(1);
		assert!((s.total[0] - 25.0).abs() < 1e-5);
		assert!((max - 25.0).abs() < 1e-5);
	}
}
