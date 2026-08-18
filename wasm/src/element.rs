use crate::metrics::boresight_cosine;
use wasm_bindgen::prelude::*;

#[allow(dead_code)]
pub const PATTERN_ISOTROPIC: u32 = 0;
pub const PATTERN_COS_N: u32 = 1;

/// Power-conserving cos^n exponent from peak element gain in dBi:
/// `n = 10^(element_gain/10)/2 - 1`, clamped at 0.
#[wasm_bindgen]
pub fn element_exponent_from_peak_dbi(gain_dbi: f32) -> f32 {
	if !gain_dbi.is_finite() {
		return 0.0;
	}
	let g = 10f32.powf(gain_dbi * 0.1);
	(g * 0.5 - 1.0).max(0.0)
}

/// Multiply AF intensity `total` (row-major `i2 * n1 + i1`) by the element
/// power pattern. Returns the new peak. Cos^n uses `[max(w,0)]^n` and is 0
/// for invisible/back directions even when `n == 0` (`0^0` would otherwise be 1).
#[wasm_bindgen]
pub fn apply_element_pattern(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &mut [f32],
	kind: u32,
	n: f32,
) -> f32 {
	let n1 = ax1.len();
	let n2 = ax2.len();
	let apply_cos = kind == PATTERN_COS_N && n.is_finite();
	let n = n.max(0.0);
	let mut max_value = f32::NEG_INFINITY;
	if apply_cos {
		for i2 in 0..n2 {
			let a2 = ax2[i2];
			let off = i2 * n1;
			for i1 in 0..n1 {
				let idx = off + i1;
				if idx >= total.len() {
					break;
				}
				let w = boresight_cosine(domain, ax1[i1], a2);
				let factor = if w <= 0.0 {
					0.0
				} else if n == 0.0 {
					1.0
				} else {
					w.powf(n)
				};
				let t = total[idx] * factor;
				total[idx] = t;
				if t > max_value {
					max_value = t;
				}
			}
		}
	} else {
		for &t in total.iter() {
			if t > max_value {
				max_value = t;
			}
		}
	}
	max_value
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::kernel::{DOMAIN_LUDWIG3, DOMAIN_SPHERICAL, DOMAIN_UV};

	#[test]
	fn exponent_matches_excel_formula() {
		let n5 = element_exponent_from_peak_dbi(5.0);
		assert!((n5 - (10f32.powf(0.5) * 0.5 - 1.0)).abs() < 1e-6);
		assert!((n5 - 0.58113883).abs() < 1e-5);

		let n_iso = element_exponent_from_peak_dbi(10.0 * 2.0f32.log10());
		assert!(n_iso.abs() < 1e-6);

		assert_eq!(element_exponent_from_peak_dbi(3.0), 0.0);
		assert_eq!(element_exponent_from_peak_dbi(f32::NAN), 0.0);
	}

	#[test]
	fn isotropic_leaves_intensity_unchanged() {
		let ax1 = [0.0f32, 0.5];
		let ax2 = [0.0f32];
		let mut total = vec![2.0f32, 4.0];
		let orig = total.clone();
		let peak = apply_element_pattern(
			DOMAIN_SPHERICAL,
			&ax1,
			&ax2,
			&mut total,
			PATTERN_ISOTROPIC,
			0.8,
		);
		assert_eq!(total, orig);
		assert!((peak - 4.0).abs() < 1e-6);
	}

	#[test]
	fn cos_n_zero_is_front_hemisphere() {
		let ax1 = [0.0f32];
		let ax2 = [0.0f32];
		let mut total = vec![3.0f32];
		let peak = apply_element_pattern(DOMAIN_SPHERICAL, &ax1, &ax2, &mut total, PATTERN_COS_N, 0.0);
		assert!((total[0] - 3.0).abs() < 1e-6);
		assert!((peak - 3.0).abs() < 1e-6);
	}

	#[test]
	fn cos_n_zeros_uv_outside_unit_circle_even_at_n_zero() {
		let u = [0.0f32, 1.2];
		let v = [0.0f32];
		let mut total = vec![1.0f32, 1.0];
		apply_element_pattern(DOMAIN_UV, &u, &v, &mut total, PATTERN_COS_N, 0.0);
		assert!((total[0] - 1.0).abs() < 1e-6);
		assert_eq!(total[1], 0.0);
	}

	#[test]
	fn spherical_cos_n_matches_analytic_w() {
		let theta = [0.0f32, std::f32::consts::FRAC_PI_4, std::f32::consts::FRAC_PI_2];
		let phi = [0.0f32];
		let mut total = vec![1.0f32; 3];
		let n = 1.0;
		apply_element_pattern(DOMAIN_SPHERICAL, &theta, &phi, &mut total, PATTERN_COS_N, n);
		for i in 0..3 {
			let w = theta[i].cos().max(0.0);
			assert!((total[i] - w.powf(n)).abs() < 1e-6);
		}
	}

	#[test]
	fn uv_cos_n_matches_analytic_w() {
		let u = [0.0f32, 0.6, 1.0];
		let v = [0.0f32, 0.8];
		let mut total = vec![1.0f32; 6];
		let n = 2.0;
		apply_element_pattern(DOMAIN_UV, &u, &v, &mut total, PATTERN_COS_N, n);
		for iv in 0..2 {
			for iu in 0..3 {
				let r2 = u[iu] * u[iu] + v[iv] * v[iv];
				let w = if r2 >= 1.0 { 0.0 } else { (1.0 - r2).sqrt() };
				let idx = iv * 3 + iu;
				assert!((total[idx] - w.powf(n)).abs() < 1e-5, "idx {idx}");
			}
		}
	}

	#[test]
	fn ludwig3_cos_n_matches_analytic_w() {
		let az = [0.0f32, std::f32::consts::FRAC_PI_4];
		let el = [0.0f32, 0.3];
		let mut total = vec![1.0f32; 4];
		let n = 0.58113883;
		apply_element_pattern(DOMAIN_LUDWIG3, &az, &el, &mut total, PATTERN_COS_N, n);
		for ie in 0..2 {
			for ia in 0..2 {
				let w = (el[ie].cos() * az[ia].cos()).max(0.0);
				let idx = ie * 2 + ia;
				assert!((total[idx] - w.powf(n)).abs() < 1e-5, "idx {idx}");
			}
		}
	}
}
