mod kernel;
mod metrics;
mod sincos;

use kernel::FarfieldState;
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
