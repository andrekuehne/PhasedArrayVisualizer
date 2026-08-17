use wide::f32x4;

#[inline(always)]
pub fn sincos_f32x4(x: f32x4) -> (f32x4, f32x4) {
	(x.sin(), x.cos())
}

#[inline(always)]
pub fn load4(buf: &[f32], i: usize) -> f32x4 {
	f32x4::from([buf[i], buf[i + 1], buf[i + 2], buf[i + 3]])
}

#[inline(always)]
pub fn store4(buf: &mut [f32], i: usize, v: f32x4) {
	let a = v.to_array();
	buf[i] = a[0];
	buf[i + 1] = a[1];
	buf[i + 2] = a[2];
	buf[i + 3] = a[3];
}

/// Accumulate `mag * exp(j * (a * t[i] + b))` into a contiguous row of length `n`.
#[inline(always)]
pub fn accumulate_linear(
	n: usize,
	mag: f32,
	a: f32,
	b: f32,
	t: &[f32],
	re: &mut [f32],
	im: &mut [f32],
) {
	let mag_v = f32x4::splat(mag);
	let a_v = f32x4::splat(a);
	let b_v = f32x4::splat(b);
	let mut i = 0;
	while i + 4 <= n {
		let phase = a_v.mul_add(load4(t, i), b_v);
		let (s, c) = sincos_f32x4(phase);
		store4(re, i, load4(re, i) + mag_v * c);
		store4(im, i, load4(im, i) + mag_v * s);
		i += 4;
	}
	while i < n {
		let phase = a.mul_add(t[i], b);
		let (s, c) = phase.sin_cos();
		re[i] += mag * c;
		im[i] += mag * s;
		i += 1;
	}
}
