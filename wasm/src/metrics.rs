use crate::kernel::{DOMAIN_LUDWIG3, DOMAIN_SPHERICAL, DOMAIN_UV};
use wasm_bindgen::prelude::*;

#[derive(Clone, Copy)]
struct Vec3 {
	x: f32,
	y: f32,
	z: f32,
}

impl Vec3 {
	fn dot(self, o: Vec3) -> f32 {
		self.x * o.x + self.y * o.y + self.z * o.z
	}

	fn cross(self, o: Vec3) -> Vec3 {
		Vec3 {
			x: self.y * o.z - self.z * o.y,
			y: self.z * o.x - self.x * o.z,
			z: self.x * o.y - self.y * o.x,
		}
	}

	fn add(self, o: Vec3) -> Vec3 {
		Vec3 {
			x: self.x + o.x,
			y: self.y + o.y,
			z: self.z + o.z,
		}
	}

	fn sub(self, o: Vec3) -> Vec3 {
		Vec3 {
			x: self.x - o.x,
			y: self.y - o.y,
			z: self.z - o.z,
		}
	}

	fn scale(self, s: f32) -> Vec3 {
		Vec3 {
			x: self.x * s,
			y: self.y * s,
			z: self.z * s,
		}
	}

	fn norm(self) -> f32 {
		self.dot(self).sqrt()
	}

	fn normalized(self) -> Option<Vec3> {
		let n = self.norm();
		if n < 1e-20 {
			None
		} else {
			Some(self.scale(1.0 / n))
		}
	}
}

/// Pattern-feature metrics extracted from a computed intensity map.
#[wasm_bindgen]
pub struct PatternMetrics {
	pub peak_i1: u32,
	pub peak_i2: u32,
	pub peak_ax1: f32,
	pub peak_ax2: f32,
	pub hpbw_ax1: f32,
	pub hpbw_ax2: f32,
	pub hpbw_ax1_deg: f32,
	pub hpbw_ax2_deg: f32,
	pub hpbw_ax1_clipped: bool,
	pub hpbw_ax2_clipped: bool,
	pub hpbw_large: f32,
	pub hpbw_small: f32,
	pub hpbw_large_deg: f32,
	pub hpbw_small_deg: f32,
	pub hpbw_large_clipped: bool,
	pub hpbw_small_clipped: bool,
	pub hpbw_large_angle_deg: f32,
	pub hpbw_small_angle_deg: f32,
	pub nearest_sll_db: f32,
	pub largest_sll_db: f32,
	pub nearest_sll_ax1: f32,
	pub nearest_sll_ax2: f32,
	pub largest_sll_ax1: f32,
	pub largest_sll_ax2: f32,
}

fn nan_metrics() -> PatternMetrics {
	PatternMetrics {
		peak_i1: 0,
		peak_i2: 0,
		peak_ax1: f32::NAN,
		peak_ax2: f32::NAN,
		hpbw_ax1: f32::NAN,
		hpbw_ax2: f32::NAN,
		hpbw_ax1_deg: f32::NAN,
		hpbw_ax2_deg: f32::NAN,
		hpbw_ax1_clipped: true,
		hpbw_ax2_clipped: true,
		hpbw_large: f32::NAN,
		hpbw_small: f32::NAN,
		hpbw_large_deg: f32::NAN,
		hpbw_small_deg: f32::NAN,
		hpbw_large_clipped: true,
		hpbw_small_clipped: true,
		hpbw_large_angle_deg: f32::NAN,
		hpbw_small_angle_deg: f32::NAN,
		nearest_sll_db: f32::NAN,
		largest_sll_db: f32::NAN,
		nearest_sll_ax1: f32::NAN,
		nearest_sll_ax2: f32::NAN,
		largest_sll_ax1: f32::NAN,
		largest_sll_ax2: f32::NAN,
	}
}

fn look(domain: u32, a1: f32, a2: f32) -> Option<Vec3> {
	match domain {
		DOMAIN_SPHERICAL => {
			let (st, ct) = a1.sin_cos();
			let (sp, cp) = a2.sin_cos();
			Some(Vec3 {
				x: st * cp,
				y: st * sp,
				z: ct,
			})
		}
		DOMAIN_UV => {
			let r2 = a1 * a1 + a2 * a2;
			if r2 >= 1.0 {
				None
			} else {
				Some(Vec3 {
					x: a1,
					y: a2,
					z: (1.0 - r2).sqrt(),
				})
			}
		}
		DOMAIN_LUDWIG3 => {
			let (sa, ca) = a1.sin_cos();
			let (se, ce) = a2.sin_cos();
			Some(Vec3 {
				x: ce * sa,
				y: se,
				z: ce * ca,
			})
		}
		_ => None,
	}
}

fn adjust_theta_phi(mut theta: f32, mut phi: f32) -> (f32, f32) {
	const HALF: f32 = std::f32::consts::FRAC_PI_2;
	const FULL: f32 = std::f32::consts::PI;
	if phi > HALF {
		phi -= FULL;
		theta = -theta;
	}
	if phi < -HALF {
		phi += FULL;
		theta = -theta;
	}
	(theta, phi)
}

fn from_look(domain: u32, r: Vec3) -> Option<(f32, f32)> {
	let n = r.norm();
	if n < 1e-20 {
		return None;
	}
	let x = r.x / n;
	let y = r.y / n;
	let z = r.z / n;
	match domain {
		DOMAIN_SPHERICAL => {
			let theta = z.clamp(-1.0, 1.0).acos();
			let phi = y.atan2(x);
			Some(adjust_theta_phi(theta, phi))
		}
		DOMAIN_UV => {
			let r2 = x * x + y * y;
			if z < 0.0 || r2 >= 1.0 {
				None
			} else {
				Some((x, y))
			}
		}
		DOMAIN_LUDWIG3 => Some((x.atan2(z), y.clamp(-1.0, 1.0).asin())),
		_ => None,
	}
}

fn sample_valid(domain: u32, ax1: &[f32], ax2: &[f32], i1: usize, i2: usize) -> bool {
	look(domain, ax1[i1], ax2[i2]).is_some()
}

fn angle_rad(a: Vec3, b: Vec3) -> f32 {
	a.dot(b).clamp(-1.0, 1.0).acos()
}

fn angle_deg(a: Vec3, b: Vec3) -> f32 {
	angle_rad(a, b).to_degrees()
}

fn lerp_look(a: Vec3, b: Vec3, t: f32) -> Option<Vec3> {
	let x = a.x + t * (b.x - a.x);
	let y = a.y + t * (b.y - a.y);
	let z = a.z + t * (b.z - a.z);
	let n = (x * x + y * y + z * z).sqrt();
	if n < 1e-20 {
		None
	} else {
		Some(Vec3 {
			x: x / n,
			y: y / n,
			z: z / n,
		})
	}
}

fn to_db(intensity: f32, peak: f32) -> f32 {
	if intensity <= 0.0 || peak <= 0.0 {
		f32::NEG_INFINITY
	} else {
		10.0 * (intensity / peak).log10()
	}
}

fn idx(i1: usize, i2: usize, n1: usize) -> usize {
	i2 * n1 + i1
}

fn nearest_index(ax: &[f32], target: f32) -> usize {
	let mut best = 0;
	let mut best_d = f32::INFINITY;
	for (i, &v) in ax.iter().enumerate() {
		let d = (v - target).abs();
		if d < best_d {
			best_d = d;
			best = i;
		}
	}
	best
}

fn orthogonal_phi_index(ax2: &[f32], phi_peak: f32) -> usize {
	let plus = phi_peak + std::f32::consts::FRAC_PI_2;
	let minus = phi_peak - std::f32::consts::FRAC_PI_2;
	let i_plus = nearest_index(ax2, plus);
	let i_minus = nearest_index(ax2, minus);
	if (ax2[i_plus] - plus).abs() <= (ax2[i_minus] - minus).abs() {
		i_plus
	} else {
		i_minus
	}
}

fn interp_cross(x0: f32, y0: f32, x1: f32, y1: f32, half: f32) -> f32 {
	let dy = y1 - y0;
	if dy.abs() < 1e-20 {
		x0
	} else {
		x0 + (half - y0) / dy * (x1 - x0)
	}
}

struct CutResult {
	native: f32,
	deg: f32,
	clipped: bool,
}

fn clipped_cut() -> CutResult {
	CutResult {
		native: f32::NAN,
		deg: f32::NAN,
		clipped: true,
	}
}

fn finish_cut(
	domain: u32,
	left_ok: bool,
	right_ok: bool,
	left_a1: f32,
	left_a2: f32,
	right_a1: f32,
	right_a2: f32,
	native: f32,
) -> CutResult {
	if !left_ok || !right_ok {
		return clipped_cut();
	}
	let Some(r_l) = look(domain, left_a1, left_a2) else {
		return clipped_cut();
	};
	let Some(r_r) = look(domain, right_a1, right_a2) else {
		return clipped_cut();
	};
	CutResult {
		native,
		deg: angle_deg(r_l, r_r),
		clipped: false,
	}
}

fn hpbw_along_ax1(
	domain: u32,
	ax1: &[f32],
	ax2_fixed: f32,
	y: impl Fn(usize) -> f32,
	valid: impl Fn(usize) -> bool,
	i_peak: usize,
	half: f32,
) -> CutResult {
	let n = ax1.len();
	if n == 0 || i_peak >= n || !valid(i_peak) {
		return clipped_cut();
	}

	let mut left_x = ax1[0];
	let mut left_ok = false;
	let mut i = i_peak;
	while i > 0 {
		let j = i - 1;
		if !valid(j) {
			break;
		}
		if y(j) < half {
			left_x = interp_cross(ax1[j], y(j), ax1[i], y(i), half);
			left_ok = true;
			break;
		}
		i = j;
	}

	let mut right_x = ax1[n - 1];
	let mut right_ok = false;
	i = i_peak;
	while i + 1 < n {
		let j = i + 1;
		if !valid(j) {
			break;
		}
		if y(j) < half {
			right_x = interp_cross(ax1[i], y(i), ax1[j], y(j), half);
			right_ok = true;
			break;
		}
		i = j;
	}

	finish_cut(
		domain,
		left_ok,
		right_ok,
		left_x,
		ax2_fixed,
		right_x,
		ax2_fixed,
		right_x - left_x,
	)
}

fn hpbw_along_ax2(
	domain: u32,
	ax1_fixed: f32,
	ax2: &[f32],
	y: impl Fn(usize) -> f32,
	valid: impl Fn(usize) -> bool,
	i_peak: usize,
	half: f32,
) -> CutResult {
	let n = ax2.len();
	if n == 0 || i_peak >= n || !valid(i_peak) {
		return clipped_cut();
	}

	let mut left_x = ax2[0];
	let mut left_ok = false;
	let mut i = i_peak;
	while i > 0 {
		let j = i - 1;
		if !valid(j) {
			break;
		}
		if y(j) < half {
			left_x = interp_cross(ax2[j], y(j), ax2[i], y(i), half);
			left_ok = true;
			break;
		}
		i = j;
	}

	let mut right_x = ax2[n - 1];
	let mut right_ok = false;
	i = i_peak;
	while i + 1 < n {
		let j = i + 1;
		if !valid(j) {
			break;
		}
		if y(j) < half {
			right_x = interp_cross(ax2[i], y(i), ax2[j], y(j), half);
			right_ok = true;
			break;
		}
		i = j;
	}

	finish_cut(
		domain,
		left_ok,
		right_ok,
		ax1_fixed,
		left_x,
		ax1_fixed,
		right_x,
		right_x - left_x,
	)
}

/// Next spherical (θ, φ) index walking φ by `dir` (±1).
///
/// The chart identifies `(θ, +90°) ≡ (-θ, -90°)` (see `adjust_theta_phi` in JS),
/// so crossing a φ edge continues on the sign-flipped θ row and skips the
/// duplicate sample.
fn step_spherical_phi(
	ax1: &[f32],
	n2: usize,
	i1: usize,
	i2: usize,
	dir: isize,
) -> Option<(usize, usize)> {
	if n2 < 2 || dir == 0 {
		return None;
	}
	let next_i2 = i2 as isize + dir;
	if next_i2 >= 0 && next_i2 < n2 as isize {
		return Some((i1, next_i2 as usize));
	}
	let i1_flip = nearest_index(ax1, -ax1[i1]);
	let i2_edge = if dir > 0 { 0usize } else { n2 - 1 };
	let i2_after = i2_edge as isize + dir;
	if i2_after < 0 || i2_after >= n2 as isize {
		return None;
	}
	Some((i1_flip, i2_after as usize))
}

struct PhiCross {
	unwrap: f32,
	r: Vec3,
}

fn walk_spherical_phi(
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	n1: usize,
	n2: usize,
	i1_peak: usize,
	i2_peak: usize,
	dir: isize,
	half: f32,
) -> Option<PhiCross> {
	let dphi = if n2 > 1 {
		((ax2[n2 - 1] - ax2[0]) / (n2 as f32 - 1.0)).abs()
	} else {
		return None;
	};
	let mut i1 = i1_peak;
	let mut i2 = i2_peak;
	let mut unwrap = 0.0f32;
	let mut y = total[idx(i1, i2, n1)];
	let mut r = look(DOMAIN_SPHERICAL, ax1[i1], ax2[i2])?;
	let mut seen = vec![false; n1 * n2];
	seen[idx(i1, i2, n1)] = true;

	for _ in 0..(2 * n2 + 4) {
		let (j1, j2) = step_spherical_phi(ax1, n2, i1, i2, dir)?;
		let k = idx(j1, j2, n1);
		if seen[k] {
			return None;
		}
		seen[k] = true;
		let y2 = total[k];
		let r2 = look(DOMAIN_SPHERICAL, ax1[j1], ax2[j2])?;
		let d_unwrap = dir as f32 * dphi;
		if y2 < half {
			let dy = y2 - y;
			let t = if dy.abs() < 1e-20 {
				0.0
			} else {
				(half - y) / dy
			};
			let r_x = lerp_look(r, r2, t)?;
			return Some(PhiCross {
				unwrap: unwrap + t * d_unwrap,
				r: r_x,
			});
		}
		unwrap += d_unwrap;
		i1 = j1;
		i2 = j2;
		y = y2;
		r = r2;
	}
	None
}

/// HPBW along φ, continuing across the ±90° chart wrap onto the −θ row.
fn hpbw_along_spherical_phi(
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	n1: usize,
	n2: usize,
	i1_peak: usize,
	i2_peak: usize,
	half: f32,
) -> CutResult {
	let Some(left) = walk_spherical_phi(
		ax1, ax2, total, n1, n2, i1_peak, i2_peak, -1, half,
	) else {
		return clipped_cut();
	};
	let Some(right) = walk_spherical_phi(
		ax1, ax2, total, n1, n2, i1_peak, i2_peak, 1, half,
	) else {
		return clipped_cut();
	};
	CutResult {
		native: right.unwrap - left.unwrap,
		deg: angle_deg(left.r, right.r),
		clipped: false,
	}
}

fn hpbw_orthogonal_meridian(
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	n1: usize,
	_n2: usize,
	i1_peak: usize,
	peak_ax2: f32,
	half: f32,
) -> CutResult {
	let i2_ortho = orthogonal_phi_index(ax2, peak_ax2);
	let mut i1_ortho = i1_peak;
	let mut best_v = f32::NEG_INFINITY;
	for i1 in 0..n1 {
		if !sample_valid(DOMAIN_SPHERICAL, ax1, ax2, i1, i2_ortho) {
			continue;
		}
		let v = total[idx(i1, i2_ortho, n1)];
		if v > best_v {
			best_v = v;
			i1_ortho = i1;
		}
	}
	hpbw_along_ax1(
		DOMAIN_SPHERICAL,
		ax1,
		ax2[i2_ortho],
		|i1| total[idx(i1, i2_ortho, n1)],
		|i1| sample_valid(DOMAIN_SPHERICAL, ax1, ax2, i1, i2_ortho),
		i1_ortho,
		half,
	)
}

fn find_peak(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	n1: usize,
	n2: usize,
) -> Option<(usize, usize, f32)> {
	let mut best_i1 = 0;
	let mut best_i2 = 0;
	let mut best = f32::NEG_INFINITY;
	let mut found = false;
	for i2 in 0..n2 {
		for i1 in 0..n1 {
			if !sample_valid(domain, ax1, ax2, i1, i2) {
				continue;
			}
			let v = total[idx(i1, i2, n1)];
			if v > best {
				best = v;
				best_i1 = i1;
				best_i2 = i2;
				found = true;
			}
		}
	}
	if found && best > 0.0 {
		Some((best_i1, best_i2, best))
	} else {
		None
	}
}

fn flood_main_beam(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	n1: usize,
	n2: usize,
	i1_peak: usize,
	i2_peak: usize,
	half: f32,
) -> Vec<bool> {
	let n = n1 * n2;
	let mut in_beam = vec![false; n];
	let mut stack = vec![(i1_peak, i2_peak)];
	in_beam[idx(i1_peak, i2_peak, n1)] = true;
	while let Some((i1, i2)) = stack.pop() {
		let mut neighbors = vec![
			(i1.wrapping_sub(1), i2),
			(i1 + 1, i2),
			(i1, i2.wrapping_sub(1)),
			(i1, i2 + 1),
		];
		if domain == DOMAIN_SPHERICAL {
			if i2 == 0 {
				if let Some(p) = step_spherical_phi(ax1, n2, i1, i2, -1) {
					neighbors.push(p);
				}
			}
			if i2 + 1 == n2 {
				if let Some(p) = step_spherical_phi(ax1, n2, i1, i2, 1) {
					neighbors.push(p);
				}
			}
		}
		for (j1, j2) in neighbors {
			if j1 >= n1 || j2 >= n2 {
				continue;
			}
			let k = idx(j1, j2, n1);
			if in_beam[k] {
				continue;
			}
			if !sample_valid(domain, ax1, ax2, j1, j2) {
				continue;
			}
			if total[k] >= half {
				in_beam[k] = true;
				stack.push((j1, j2));
			}
		}
	}
	in_beam
}

fn is_local_max(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	in_beam: &[bool],
	n1: usize,
	n2: usize,
	i1: usize,
	i2: usize,
) -> bool {
	if !sample_valid(domain, ax1, ax2, i1, i2) {
		return false;
	}
	let k = idx(i1, i2, n1);
	if in_beam[k] {
		return false;
	}
	let v = total[k];
	for d2 in -1isize..=1 {
		for d1 in -1isize..=1 {
			if d1 == 0 && d2 == 0 {
				continue;
			}
			let j1 = i1 as isize + d1;
			let j2 = i2 as isize + d2;
			if j1 < 0 || j2 < 0 || j1 >= n1 as isize || j2 >= n2 as isize {
				continue;
			}
			let u1 = j1 as usize;
			let u2 = j2 as usize;
			if !sample_valid(domain, ax1, ax2, u1, u2) {
				continue;
			}
			if total[idx(u1, u2, n1)] >= v {
				return false;
			}
		}
	}
	if domain == DOMAIN_SPHERICAL {
		for dir in [-1isize, 1] {
			if let Some((u1, u2)) = step_spherical_phi(ax1, n2, i1, i2, dir) {
				if total[idx(u1, u2, n1)] >= v {
					return false;
				}
			}
		}
	}
	true
}

const HPBW_SWEEP_N: usize = 36;

fn geodesic_point(r0: Vec3, d: Vec3, psi: f32) -> Option<Vec3> {
	let (s, c) = psi.sin_cos();
	r0.scale(c).add(d.scale(s)).normalized()
}

fn axis_frac(ax: &[f32], val: f32) -> Option<(usize, f32)> {
	let n = ax.len();
	if n < 2 {
		return None;
	}
	let lo = ax[0];
	let hi = ax[n - 1];
	let span = hi - lo;
	if span.abs() < 1e-20 {
		return None;
	}
	let mut v = val;
	let tol = span.abs() * 1e-4 + 1e-6;
	if span > 0.0 {
		if v < lo - tol || v > hi + tol {
			return None;
		}
		v = v.clamp(lo, hi);
	} else {
		if v > lo + tol || v < hi - tol {
			return None;
		}
		v = v.clamp(hi, lo);
	}
	let f = (v - lo) / span * (n as f32 - 1.0);
	let i = (f.floor() as usize).min(n - 2);
	let t = (f - i as f32).clamp(0.0, 1.0);
	Some((i, t))
}

fn bilinear_intensity(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	a1: f32,
	a2: f32,
) -> Option<f32> {
	if look(domain, a1, a2).is_none() {
		return None;
	}
	let (i, t1) = axis_frac(ax1, a1)?;
	let (j, t2) = axis_frac(ax2, a2)?;
	let n1 = ax1.len();
	let y00 = total[idx(i, j, n1)];
	let y10 = total[idx(i + 1, j, n1)];
	let y01 = total[idx(i, j + 1, n1)];
	let y11 = total[idx(i + 1, j + 1, n1)];
	let y0 = y00 + t1 * (y10 - y00);
	let y1 = y01 + t1 * (y11 - y01);
	Some(y0 + t2 * (y1 - y0))
}

fn tangent_basis(
	domain: u32,
	ax1: &[f32],
	_ax2: &[f32],
	a1: f32,
	a2: f32,
) -> Option<(Vec3, Vec3, Vec3)> {
	let r0 = look(domain, a1, a2)?;
	let n1 = ax1.len();
	let da = if n1 > 1 {
		(ax1[n1 - 1] - ax1[0]) / (n1 as f32 - 1.0)
	} else {
		1e-3
	};
	let da = if da.abs() < 1e-8 {
		1e-3f32.copysign(if da == 0.0 { 1.0 } else { da })
	} else {
		da
	};
	let (r1, sign) = if let Some(r) = look(domain, a1 + da, a2) {
		(r, 1.0f32)
	} else {
		(look(domain, a1 - da, a2)?, -1.0f32)
	};
	let mut e1 = r1.sub(r0).scale(sign);
	e1 = e1.sub(r0.scale(e1.dot(r0)));
	let e1 = e1.normalized()?;
	let e2 = r0.cross(e1).normalized()?;
	Some((r0, e1, e2))
}

fn local_psi_step(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	i1: usize,
	i2: usize,
	r0: Vec3,
) -> f32 {
	let n1 = ax1.len();
	let n2 = ax2.len();
	let mut best = f32::INFINITY;
	if i1 + 1 < n1 {
		if let Some(r) = look(domain, ax1[i1 + 1], ax2[i2]) {
			best = best.min(angle_rad(r0, r));
		}
	}
	if i1 > 0 {
		if let Some(r) = look(domain, ax1[i1 - 1], ax2[i2]) {
			best = best.min(angle_rad(r0, r));
		}
	}
	if i2 + 1 < n2 {
		if let Some(r) = look(domain, ax1[i1], ax2[i2 + 1]) {
			best = best.min(angle_rad(r0, r));
		}
	}
	if i2 > 0 {
		if let Some(r) = look(domain, ax1[i1], ax2[i2 - 1]) {
			best = best.min(angle_rad(r0, r));
		}
	}
	(0.5 * best).clamp(1e-4, 0.05)
}

struct GeodesicCross {
	psi: f32,
	a1: f32,
	a2: f32,
}

fn walk_geodesic(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	r0: Vec3,
	d: Vec3,
	half: f32,
	y0: f32,
	dpsi: f32,
) -> Option<GeodesicCross> {
	let psi_max = std::f32::consts::PI;
	let mut psi = 0.0f32;
	let mut y = y0;
	let max_steps = ((psi_max / dpsi).ceil() as usize + 2).min(4096);
	for _ in 0..max_steps {
		let psi2 = psi + dpsi;
		if psi2 > psi_max {
			return None;
		}
		let r2 = geodesic_point(r0, d, psi2)?;
		let (a1, a2) = from_look(domain, r2)?;
		let y2 = bilinear_intensity(domain, ax1, ax2, total, a1, a2)?;
		if y2 < half {
			let psi_x = interp_cross(psi, y, psi2, y2, half);
			let r_x = geodesic_point(r0, d, psi_x)?;
			let (a1_x, a2_x) = from_look(domain, r_x)?;
			return Some(GeodesicCross {
				psi: psi_x,
				a1: a1_x,
				a2: a2_x,
			});
		}
		psi = psi2;
		y = y2;
	}
	None
}

struct ExtremaHpbw {
	large: f32,
	small: f32,
	large_deg: f32,
	small_deg: f32,
	large_angle_deg: f32,
	small_angle_deg: f32,
	clipped: bool,
}

fn clipped_extrema() -> ExtremaHpbw {
	ExtremaHpbw {
		large: f32::NAN,
		small: f32::NAN,
		large_deg: f32::NAN,
		small_deg: f32::NAN,
		large_angle_deg: f32::NAN,
		small_angle_deg: f32::NAN,
		clipped: true,
	}
}

fn hpbw_large_small(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
	i1_peak: usize,
	i2_peak: usize,
	peak_ax1: f32,
	peak_ax2: f32,
	half: f32,
) -> ExtremaHpbw {
	let Some((r0, e1, e2)) = tangent_basis(domain, ax1, ax2, peak_ax1, peak_ax2) else {
		return clipped_extrema();
	};
	let y0 = total[idx(i1_peak, i2_peak, ax1.len())];
	let dpsi = local_psi_step(domain, ax1, ax2, i1_peak, i2_peak, r0);
	let mut large_deg = f32::NEG_INFINITY;
	let mut small_deg = f32::INFINITY;
	let mut large = f32::NAN;
	let mut small = f32::NAN;
	let mut large_angle = f32::NAN;
	let mut small_angle = f32::NAN;
	let mut found = false;

	for k in 0..HPBW_SWEEP_N {
		let alpha = k as f32 * std::f32::consts::PI / HPBW_SWEEP_N as f32;
		let (s, c) = alpha.sin_cos();
		let d = e1.scale(c).add(e2.scale(s));
		let Some(left) = walk_geodesic(domain, ax1, ax2, total, r0, d.scale(-1.0), half, y0, dpsi)
		else {
			continue;
		};
		let Some(right) = walk_geodesic(domain, ax1, ax2, total, r0, d, half, y0, dpsi) else {
			continue;
		};
		let deg = (left.psi + right.psi).to_degrees();
		let native = if domain == DOMAIN_UV {
			let du = right.a1 - left.a1;
			let dv = right.a2 - left.a2;
			(du * du + dv * dv).sqrt()
		} else {
			f32::NAN
		};
		let angle_deg = alpha.to_degrees();
		if !found || deg > large_deg {
			large_deg = deg;
			large = native;
			large_angle = angle_deg;
		}
		if !found || deg < small_deg {
			small_deg = deg;
			small = native;
			small_angle = angle_deg;
		}
		found = true;
	}

	if !found {
		return clipped_extrema();
	}
	ExtremaHpbw {
		large,
		small,
		large_deg,
		small_deg,
		large_angle_deg: large_angle,
		small_angle_deg: small_angle,
		clipped: false,
	}
}

/// Extract HPBW and sidelobe metrics from a computed far-field intensity map.
///
/// `total` is row-major with `index = i2 * n1 + i1` (ax2 outer), matching
/// `FarfieldKernel::take_total`.
pub fn compute_pattern_metrics(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
) -> PatternMetrics {
	let n1 = ax1.len();
	let n2 = ax2.len();
	if n1 == 0 || n2 == 0 || total.len() != n1 * n2 {
		return nan_metrics();
	}
	if !matches!(domain, DOMAIN_SPHERICAL | DOMAIN_UV | DOMAIN_LUDWIG3) {
		return nan_metrics();
	}

	let Some((i1_peak, i2_peak, peak)) = find_peak(domain, ax1, ax2, total, n1, n2) else {
		return nan_metrics();
	};
	let half = peak * 0.5;
	let peak_ax1 = ax1[i1_peak];
	let peak_ax2 = ax2[i2_peak];

	let cut_ax1 = hpbw_along_ax1(
		domain,
		ax1,
		peak_ax2,
		|i1| total[idx(i1, i2_peak, n1)],
		|i1| sample_valid(domain, ax1, ax2, i1, i2_peak),
		i1_peak,
		half,
	);

	let d_ax1 = if n1 > 1 {
		(ax1[n1 - 1] - ax1[0]) / (n1 as f32 - 1.0)
	} else {
		0.0
	};
	let pole = domain == DOMAIN_SPHERICAL && peak_ax1.abs() <= d_ax1.abs().max(1e-6);

	let cut_ax2 = if domain == DOMAIN_SPHERICAL {
		if pole {
			hpbw_orthogonal_meridian(ax1, ax2, total, n1, n2, i1_peak, peak_ax2, half)
		} else {
			let wrapped =
				hpbw_along_spherical_phi(ax1, ax2, total, n1, n2, i1_peak, i2_peak, half);
			if wrapped.clipped {
				hpbw_orthogonal_meridian(ax1, ax2, total, n1, n2, i1_peak, peak_ax2, half)
			} else {
				wrapped
			}
		}
	} else {
		hpbw_along_ax2(
			domain,
			peak_ax1,
			ax2,
			|i2| total[idx(i1_peak, i2, n1)],
			|i2| sample_valid(domain, ax1, ax2, i1_peak, i2),
			i2_peak,
			half,
		)
	};

	let extrema = hpbw_large_small(
		domain,
		ax1,
		ax2,
		total,
		i1_peak,
		i2_peak,
		peak_ax1,
		peak_ax2,
		half,
	);

	let in_beam = flood_main_beam(
		domain, ax1, ax2, total, n1, n2, i1_peak, i2_peak, half,
	);

	let r_peak = look(domain, peak_ax1, peak_ax2);
	let mut nearest_i1 = 0usize;
	let mut nearest_i2 = 0usize;
	let mut nearest_ang = f32::INFINITY;
	let mut nearest_i = 0.0f32;
	let mut largest_i1 = 0usize;
	let mut largest_i2 = 0usize;
	let mut largest_i = f32::NEG_INFINITY;
	let mut found_local = false;
	let mut found_unmasked = false;
	let mut fallback_i1 = 0usize;
	let mut fallback_i2 = 0usize;
	let mut fallback_i = f32::NEG_INFINITY;

	for i2 in 0..n2 {
		for i1 in 0..n1 {
			if !sample_valid(domain, ax1, ax2, i1, i2) {
				continue;
			}
			let k = idx(i1, i2, n1);
			if in_beam[k] {
				continue;
			}
			found_unmasked = true;
			let v = total[k];
			if v > fallback_i {
				fallback_i = v;
				fallback_i1 = i1;
				fallback_i2 = i2;
			}
			if !is_local_max(domain, ax1, ax2, total, &in_beam, n1, n2, i1, i2) {
				continue;
			}
			found_local = true;
			if v > largest_i {
				largest_i = v;
				largest_i1 = i1;
				largest_i2 = i2;
			}
			if let Some(rp) = r_peak {
				if let Some(rl) = look(domain, ax1[i1], ax2[i2]) {
					let ang = angle_deg(rp, rl);
					if ang < nearest_ang {
						nearest_ang = ang;
						nearest_i1 = i1;
						nearest_i2 = i2;
						nearest_i = v;
					}
				}
			}
		}
	}

	let (nearest_sll_db, nearest_sll_ax1, nearest_sll_ax2, largest_sll_db, largest_sll_ax1, largest_sll_ax2) =
		if !found_unmasked {
			(
				f32::NAN,
				f32::NAN,
				f32::NAN,
				f32::NAN,
				f32::NAN,
				f32::NAN,
			)
		} else if found_local {
			(
				to_db(nearest_i, peak),
				ax1[nearest_i1],
				ax2[nearest_i2],
				to_db(largest_i, peak),
				ax1[largest_i1],
				ax2[largest_i2],
			)
		} else {
			let db = to_db(fallback_i, peak);
			(
				db,
				ax1[fallback_i1],
				ax2[fallback_i2],
				db,
				ax1[fallback_i1],
				ax2[fallback_i2],
			)
		};

	PatternMetrics {
		peak_i1: i1_peak as u32,
		peak_i2: i2_peak as u32,
		peak_ax1,
		peak_ax2,
		hpbw_ax1: cut_ax1.native,
		hpbw_ax2: cut_ax2.native,
		hpbw_ax1_deg: cut_ax1.deg,
		hpbw_ax2_deg: cut_ax2.deg,
		hpbw_ax1_clipped: cut_ax1.clipped,
		hpbw_ax2_clipped: cut_ax2.clipped,
		hpbw_large: extrema.large,
		hpbw_small: extrema.small,
		hpbw_large_deg: extrema.large_deg,
		hpbw_small_deg: extrema.small_deg,
		hpbw_large_clipped: extrema.clipped,
		hpbw_small_clipped: extrema.clipped,
		hpbw_large_angle_deg: extrema.large_angle_deg,
		hpbw_small_angle_deg: extrema.small_angle_deg,
		nearest_sll_db,
		largest_sll_db,
		nearest_sll_ax1,
		nearest_sll_ax2,
		largest_sll_ax1,
		largest_sll_ax2,
	}
}

#[wasm_bindgen]
pub fn extract_pattern_metrics(
	domain: u32,
	ax1: &[f32],
	ax2: &[f32],
	total: &[f32],
) -> PatternMetrics {
	compute_pattern_metrics(domain, ax1, ax2, total)
}

#[cfg(test)]
mod tests {
	use super::*;

	fn linspace(a: f32, b: f32, n: usize) -> Vec<f32> {
		if n == 1 {
			return vec![a];
		}
		(0..n)
			.map(|i| a + (b - a) * (i as f32) / (n as f32 - 1.0))
			.collect()
	}

	fn add_gaussian(
		total: &mut [f32],
		ax1: &[f32],
		ax2: &[f32],
		c1: f32,
		c2: f32,
		amp: f32,
		sigma: f32,
	) {
		let n1 = ax1.len();
		let n2 = ax2.len();
		let s2 = 2.0 * sigma * sigma;
		for i2 in 0..n2 {
			for i1 in 0..n1 {
				let d1 = ax1[i1] - c1;
				let d2 = ax2[i2] - c2;
				total[i2 * n1 + i1] += amp * (-(d1 * d1 + d2 * d2) / s2).exp();
			}
		}
	}

	fn expected_hpbw(sigma: f32) -> f32 {
		2.0 * sigma * (2.0 * std::f32::consts::LN_2).sqrt()
	}

	#[test]
	fn gaussian_hpbw_ludwig3() {
		let n = 101;
		let ax1 = linspace(-0.4, 0.4, n);
		let ax2 = linspace(-0.4, 0.4, n);
		let mut total = vec![0.0f32; n * n];
		let sigma = 0.05;
		add_gaussian(&mut total, &ax1, &ax2, 0.0, 0.0, 1.0, sigma);
		let m = compute_pattern_metrics(DOMAIN_LUDWIG3, &ax1, &ax2, &total);
		let want = expected_hpbw(sigma);
		assert!(!m.hpbw_ax1_clipped && !m.hpbw_ax2_clipped);
		assert!(
			(m.hpbw_ax1 - want).abs() / want < 0.05,
			"hpbw_ax1 {} vs {}",
			m.hpbw_ax1,
			want
		);
		assert!(
			(m.hpbw_ax2 - want).abs() / want < 0.05,
			"hpbw_ax2 {} vs {}",
			m.hpbw_ax2,
			want
		);
		let want_deg = want.to_degrees();
		assert!((m.hpbw_ax1_deg - want_deg).abs() / want_deg < 0.08);
		assert!((m.hpbw_ax2_deg - want_deg).abs() / want_deg < 0.08);
		assert!(!m.hpbw_large_clipped && !m.hpbw_small_clipped);
		assert!((m.hpbw_large_deg - want_deg).abs() / want_deg < 0.08);
		assert!((m.hpbw_small_deg - want_deg).abs() / want_deg < 0.08);
		assert!(m.peak_ax1.abs() < 1e-5 && m.peak_ax2.abs() < 1e-5);
	}

	#[test]
	fn gaussian_hpbw_uv() {
		let n = 101;
		let ax1 = linspace(-0.4, 0.4, n);
		let ax2 = linspace(-0.4, 0.4, n);
		let mut total = vec![0.0f32; n * n];
		let sigma = 0.05;
		add_gaussian(&mut total, &ax1, &ax2, 0.0, 0.0, 1.0, sigma);
		let m = compute_pattern_metrics(DOMAIN_UV, &ax1, &ax2, &total);
		let want = expected_hpbw(sigma);
		assert!(!m.hpbw_ax1_clipped && !m.hpbw_ax2_clipped);
		assert!((m.hpbw_ax1 - want).abs() / want < 0.05);
		assert!((m.hpbw_ax2 - want).abs() / want < 0.05);
		assert!(!m.hpbw_large_clipped && !m.hpbw_small_clipped);
		let want_deg = want.to_degrees();
		assert!((m.hpbw_large_deg - want_deg).abs() / want_deg < 0.08);
		assert!((m.hpbw_small_deg - want_deg).abs() / want_deg < 0.08);
		assert!((m.hpbw_large - want).abs() / want < 0.08);
		assert!((m.hpbw_small - want).abs() / want < 0.08);
	}

	#[test]
	fn two_gaussians_nearest_and_largest_sll() {
		let n = 121;
		let ax1 = linspace(-0.6, 0.6, n);
		let ax2 = linspace(-0.6, 0.6, n);
		let mut total = vec![0.0f32; n * n];
		add_gaussian(&mut total, &ax1, &ax2, 0.0, 0.0, 1.0, 0.04);
		add_gaussian(&mut total, &ax1, &ax2, 0.18, 0.0, 0.1, 0.03);
		add_gaussian(&mut total, &ax1, &ax2, 0.35, 0.25, 0.25, 0.03);
		let m = compute_pattern_metrics(DOMAIN_LUDWIG3, &ax1, &ax2, &total);
		let nearest_db = 10.0 * 0.1f32.log10();
		let largest_db = 10.0 * 0.25f32.log10();
		assert!(
			(m.nearest_sll_db - nearest_db).abs() < 0.5,
			"nearest {} vs {}",
			m.nearest_sll_db,
			nearest_db
		);
		assert!(
			(m.largest_sll_db - largest_db).abs() < 0.5,
			"largest {} vs {}",
			m.largest_sll_db,
			largest_db
		);
		assert!((m.nearest_sll_ax1 - 0.18).abs() < 0.02);
		assert!(m.nearest_sll_ax2.abs() < 0.02);
		assert!((m.largest_sll_ax1 - 0.35).abs() < 0.02);
		assert!((m.largest_sll_ax2 - 0.25).abs() < 0.02);
		assert!(m.largest_sll_db > m.nearest_sll_db);
	}

	#[test]
	fn spherical_pole_uses_orthogonal_meridian() {
		let n = 101;
		let ax1 = linspace(-std::f32::consts::FRAC_PI_2, std::f32::consts::FRAC_PI_2, n);
		let ax2 = linspace(-std::f32::consts::FRAC_PI_2, std::f32::consts::FRAC_PI_2, n);
		let mut total = vec![0.0f32; n * n];
		let sigma = 0.08;
		for i2 in 0..n {
			for i1 in 0..n {
				let th = ax1[i1];
				total[i2 * n + i1] = (-(th * th) / (2.0 * sigma * sigma)).exp();
			}
		}
		let m = compute_pattern_metrics(DOMAIN_SPHERICAL, &ax1, &ax2, &total);
		assert!(m.peak_ax1.abs() < ax1[1] - ax1[0]);
		assert!(!m.hpbw_ax1_clipped);
		assert!(!m.hpbw_ax2_clipped, "phi cut at pole must not be used");
		let want = expected_hpbw(sigma);
		assert!(
			(m.hpbw_ax1 - want).abs() / want < 0.08,
			"hpbw_ax1 {} vs {}",
			m.hpbw_ax1,
			want
		);
		assert!(
			(m.hpbw_ax2 - want).abs() / want < 0.08,
			"orthogonal meridian {} vs {}",
			m.hpbw_ax2,
			want
		);
		assert!(!m.hpbw_large_clipped && !m.hpbw_small_clipped);
		assert!(
			(m.hpbw_large_deg - m.hpbw_ax1_deg).abs() / m.hpbw_ax1_deg < 0.12,
			"large {} vs ax1 {}",
			m.hpbw_large_deg,
			m.hpbw_ax1_deg
		);
		assert!(
			(m.hpbw_small_deg - m.hpbw_ax1_deg).abs() / m.hpbw_ax1_deg < 0.12,
			"small {} vs ax1 {}",
			m.hpbw_small_deg,
			m.hpbw_ax1_deg
		);
	}

	fn fill_spherical_cap(total: &mut [f32], ax1: &[f32], ax2: &[f32], th0: f32, ph0: f32, sigma: f32) {
		let n1 = ax1.len();
		let n2 = ax2.len();
		let r0 = look(DOMAIN_SPHERICAL, th0, ph0).unwrap();
		for i2 in 0..n2 {
			for i1 in 0..n1 {
				let r = look(DOMAIN_SPHERICAL, ax1[i1], ax2[i2]).unwrap();
				let ang = angle_deg(r0, r).to_radians();
				total[i2 * n1 + i1] = (-ang * ang / (2.0 * sigma * sigma)).exp();
			}
		}
	}

	#[test]
	fn spherical_phi_wrap_at_edge_is_not_clipped() {
		let n = 101;
		let ax1 = linspace(-std::f32::consts::FRAC_PI_2, std::f32::consts::FRAC_PI_2, n);
		let ax2 = linspace(-std::f32::consts::FRAC_PI_2, std::f32::consts::FRAC_PI_2, n);
		let sigma = 0.08;
		let th0 = 0.3;
		let mut total_edge = vec![0.0f32; n * n];
		fill_spherical_cap(&mut total_edge, &ax1, &ax2, th0, std::f32::consts::FRAC_PI_2, sigma);
		let edge = compute_pattern_metrics(DOMAIN_SPHERICAL, &ax1, &ax2, &total_edge);
		assert!(
			!edge.hpbw_ax2_clipped,
			"phi HPBW clipped at wrap: peak φ={}",
			edge.peak_ax2
		);
		assert!(edge.peak_ax2.abs() > 1.0, "peak should sit on a φ edge");

		let mut total_mid = vec![0.0f32; n * n];
		fill_spherical_cap(&mut total_mid, &ax1, &ax2, th0, 0.0, sigma);
		let mid = compute_pattern_metrics(DOMAIN_SPHERICAL, &ax1, &ax2, &total_mid);
		assert!(!mid.hpbw_ax2_clipped);
		assert!(
			(edge.hpbw_ax2_deg - mid.hpbw_ax2_deg).abs() / mid.hpbw_ax2_deg < 0.15,
			"wrap HPBW {} deg vs interior {} deg",
			edge.hpbw_ax2_deg,
			mid.hpbw_ax2_deg
		);
		assert!(!edge.hpbw_large_clipped && !edge.hpbw_small_clipped);
		assert!(!mid.hpbw_large_clipped && !mid.hpbw_small_clipped);
		assert!(
			(mid.hpbw_large_deg - mid.hpbw_small_deg).abs() / mid.hpbw_large_deg < 0.12,
			"circular cap large {} vs small {}",
			mid.hpbw_large_deg,
			mid.hpbw_small_deg
		);
	}

	#[test]
	fn uv_invisible_samples_ignored_for_peak() {
		let n = 5;
		let ax1 = linspace(-1.2, 1.2, n);
		let ax2 = linspace(-1.2, 1.2, n);
		let mut total = vec![0.0f32; n * n];
		for i2 in 0..n {
			for i1 in 0..n {
				let u = ax1[i1];
				let v = ax2[i2];
				if u * u + v * v >= 1.0 {
					total[i2 * n + i1] = 100.0;
				}
			}
		}
		let mid = n / 2;
		total[mid * n + mid] = 1.0;
		let m = compute_pattern_metrics(DOMAIN_UV, &ax1, &ax2, &total);
		assert_eq!(m.peak_i1, mid as u32);
		assert_eq!(m.peak_i2, mid as u32);
		assert!(m.peak_ax1.abs() < 1e-5);
		assert!(m.peak_ax2.abs() < 1e-5);
	}

	#[test]
	fn clipped_beam_returns_nan_and_flags() {
		let n = 11;
		let ax1 = linspace(-0.1, 0.1, n);
		let ax2 = linspace(-0.1, 0.1, n);
		let total = vec![1.0f32; n * n];
		let m = compute_pattern_metrics(DOMAIN_LUDWIG3, &ax1, &ax2, &total);
		assert!(m.hpbw_ax1_clipped && m.hpbw_ax2_clipped);
		assert!(m.hpbw_ax1.is_nan() && m.hpbw_ax2.is_nan());
		assert!(m.hpbw_ax1_deg.is_nan() && m.hpbw_ax2_deg.is_nan());
		assert!(m.hpbw_large_clipped && m.hpbw_small_clipped);
		assert!(m.hpbw_large_deg.is_nan() && m.hpbw_small_deg.is_nan());
		assert!(m.nearest_sll_db.is_nan());
		assert!(m.largest_sll_db.is_nan());
	}

	#[test]
	fn empty_inputs_return_nan() {
		let m = compute_pattern_metrics(DOMAIN_UV, &[], &[], &[]);
		assert!(m.peak_ax1.is_nan());
		let ax1 = [0.0f32];
		let ax2 = [0.0f32];
		let m = compute_pattern_metrics(DOMAIN_UV, &ax1, &ax2, &[1.0, 2.0]);
		assert!(m.peak_ax1.is_nan());
	}

	fn abs_angle_diff_180(a: f32, b: f32) -> f32 {
		let mut d = (a - b) % 180.0;
		if d > 90.0 {
			d -= 180.0;
		}
		if d < -90.0 {
			d += 180.0;
		}
		d.abs()
	}

	fn add_rotated_gaussian(
		total: &mut [f32],
		ax1: &[f32],
		ax2: &[f32],
		sigma1: f32,
		sigma2: f32,
		rot: f32,
	) {
		let n1 = ax1.len();
		let n2 = ax2.len();
		let (s, c) = rot.sin_cos();
		for i2 in 0..n2 {
			for i1 in 0..n1 {
				let d1 = ax1[i1];
				let d2 = ax2[i2];
				let p1 = c * d1 + s * d2;
				let p2 = -s * d1 + c * d2;
				total[i2 * n1 + i1] =
					(-(p1 * p1 / (2.0 * sigma1 * sigma1) + p2 * p2 / (2.0 * sigma2 * sigma2))).exp();
			}
		}
	}

	#[test]
	fn rotated_ellipse_large_small_hpbw_ludwig3() {
		let n = 121;
		let ax1 = linspace(-0.4, 0.4, n);
		let ax2 = linspace(-0.4, 0.4, n);
		let mut total = vec![0.0f32; n * n];
		let sigma1 = 0.08;
		let sigma2 = 0.03;
		let rot = 30.0f32.to_radians();
		add_rotated_gaussian(&mut total, &ax1, &ax2, sigma1, sigma2, rot);
		let m = compute_pattern_metrics(DOMAIN_LUDWIG3, &ax1, &ax2, &total);
		assert!(!m.hpbw_large_clipped && !m.hpbw_small_clipped);
		assert!(!m.hpbw_ax1_clipped && !m.hpbw_ax2_clipped);
		let want_large = expected_hpbw(sigma1).to_degrees();
		let want_small = expected_hpbw(sigma2).to_degrees();
		assert!(
			(m.hpbw_large_deg - want_large).abs() / want_large < 0.08,
			"large {} vs {}",
			m.hpbw_large_deg,
			want_large
		);
		assert!(
			(m.hpbw_small_deg - want_small).abs() / want_small < 0.08,
			"small {} vs {}",
			m.hpbw_small_deg,
			want_small
		);
		let ax_min = m.hpbw_ax1_deg.min(m.hpbw_ax2_deg);
		let ax_max = m.hpbw_ax1_deg.max(m.hpbw_ax2_deg);
		assert!(
			ax_min > m.hpbw_small_deg * 0.9 && ax_max < m.hpbw_large_deg * 1.1,
			"aligned ({}, {}) not between small {} and large {}",
			m.hpbw_ax1_deg,
			m.hpbw_ax2_deg,
			m.hpbw_small_deg,
			m.hpbw_large_deg
		);
		assert!(
			abs_angle_diff_180(m.hpbw_large_angle_deg, 30.0) < 8.0,
			"large angle {}",
			m.hpbw_large_angle_deg
		);
		assert!(
			abs_angle_diff_180(m.hpbw_small_angle_deg, 120.0) < 8.0,
			"small angle {}",
			m.hpbw_small_angle_deg
		);
	}
}
