mod bessel;
mod element;
mod kernel;
mod match_s;
mod metrics;
mod prad;
mod quadrature;
mod sincos;

use kernel::FarfieldState;
use prad::PradState;
use wasm_bindgen::prelude::*;

pub use kernel::{DOMAIN_LUDWIG3, DOMAIN_SPHERICAL, DOMAIN_UV};

#[wasm_bindgen]
pub struct FarfieldKernel {
	state: FarfieldState,
}

#[wasm_bindgen]
impl FarfieldKernel {
	#[wasm_bindgen(constructor)]
	pub fn new() -> FarfieldKernel {
		FarfieldKernel {
			state: FarfieldState::new(),
		}
	}

	pub fn prepare(&mut self, n1: u32, n2: u32) {
		self.state.prepare(n1, n2);
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
		self.state.set_inputs(x, y, mag, pha, ax1, ax2);
	}

	pub fn accumulate_tile(&mut self, domain: u32, frequency_scale: f32, row0: u32, row_count: u32) {
		self.state
			.accumulate_tile(domain, frequency_scale, row0, row_count);
	}

	pub fn finalize(&mut self, n_elements: u32) -> f32 {
		self.state.finalize(n_elements)
	}

	pub fn take_total(&self) -> Vec<f32> {
		self.state.total.clone()
	}
}

#[wasm_bindgen]
pub struct RadiatedPowerKernel {
	state: PradState,
}

#[wasm_bindgen]
impl RadiatedPowerKernel {
	#[wasm_bindgen(constructor)]
	pub fn new() -> RadiatedPowerKernel {
		RadiatedPowerKernel {
			state: PradState::new(),
		}
	}

	pub fn set_quadrature(&mut self, n_mu: u32, n_phi: u32) {
		self.state.set_quadrature(n_mu, n_phi);
	}

	pub fn fill_isolated(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		self.state
			.fill_isolated(x, y, frequency_scale, element_kind, element_n);
	}

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
		self.state.fill_isolated_range(
			x,
			y,
			frequency_scale,
			element_kind,
			element_n,
			sample0,
			sample_count,
		);
	}

	pub fn form_gram(&mut self) {
		self.state.form_gram();
	}

	pub fn compute(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		self.state
			.compute(x, y, frequency_scale, element_kind, element_n);
	}

	pub fn compute_j0(
		&mut self,
		x: &[f32],
		y: &[f32],
		frequency_scale: f32,
		element_kind: u32,
		element_n: f32,
	) {
		self.state
			.compute_j0(x, y, frequency_scale, element_kind, element_n);
	}

	pub fn take_re(&self) -> Vec<f32> {
		self.state.p_re.clone()
	}

	pub fn take_im(&self) -> Vec<f32> {
		self.state.p_im.clone()
	}

	pub fn form_matched_s(
		&mut self,
		z_ref: f32,
		x: &[f32],
		y: &[f32],
		x_nn: f32,
		alpha: f32,
	) {
		self.state.form_matched_s(z_ref, x, y, x_nn, alpha);
	}

	pub fn take_z0(&self) -> Vec<f64> {
		self.state.z0.clone()
	}

	pub fn take_z0_im(&self) -> Vec<f64> {
		self.state.z0_im.clone()
	}

	pub fn take_z_re(&self) -> Vec<f64> {
		self.state.r_re.clone()
	}

	pub fn take_z_im(&self) -> Vec<f64> {
		self.state.r_im.clone()
	}

	pub fn take_s_re(&self) -> Vec<f64> {
		self.state.s_re.clone()
	}

	pub fn take_s_im(&self) -> Vec<f64> {
		self.state.s_im.clone()
	}

	pub fn take_t_re(&self) -> Vec<f64> {
		self.state.t_re.clone()
	}

	pub fn take_t_im(&self) -> Vec<f64> {
		self.state.t_im.clone()
	}

	pub fn match_iterations(&self) -> u32 {
		self.state.match_iterations
	}

	pub fn match_residual(&self) -> f64 {
		self.state.match_residual
	}

	pub fn n_samples(&self) -> u32 {
		self.state.n_samples() as u32
	}

	pub fn n_elements(&self) -> u32 {
		self.state.n_elements() as u32
	}
}
