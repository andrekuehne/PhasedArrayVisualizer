//! f32 Bessel J0 (fdlibm / musl `j0f`).
//!
//! origin: FreeBSD /usr/src/lib/msun/src/e_j0f.c
//! Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
//!
//! ====================================================
//! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//!
//! Developed at SunPro, a Sun Microsystems, Inc. business.
//! Permission to use, copy, modify, and distribute this
//! software is freely granted, provided that this notice
//! is preserved.
//! ====================================================

const INVSQRTPI: f32 = 5.6418961287e-01; /* 0x3f106ebb */

/* R0/S0 on [0, 2] */
const R02: f32 = 1.5625000000e-02; /* 0x3c800000 */
const R03: f32 = -1.8997929874e-04; /* 0xb947352e */
const R04: f32 = 1.8295404516e-06; /* 0x35f58e88 */
const R05: f32 = -4.6183270541e-09; /* 0xb19eaf3c */
const S01: f32 = 1.5619102865e-02; /* 0x3c7fe744 */
const S02: f32 = 1.1692678527e-04; /* 0x38f53697 */
const S03: f32 = 5.1354652442e-07; /* 0x3509daa6 */
const S04: f32 = 1.1661400734e-09; /* 0x30a045e8 */

/// J0(x) for f32. Even; J0(0) = 1.
pub fn j0f(mut x: f32) -> f32 {
	let mut ix = x.to_bits();
	ix &= 0x7fffffff;
	if ix >= 0x7f800000 {
		return 1.0 / (x * x);
	}
	x = x.abs();
	if ix >= 0x40000000 {
		/* |x| >= 2 */
		return j0_large(ix, x);
	}
	if ix >= 0x3a000000 {
		/* |x| >= 2^-11 */
		let z = x * x;
		let r = z * (R02 + z * (R03 + z * (R04 + z * R05)));
		let s = 1.0 + z * (S01 + z * (S02 + z * (S03 + z * S04)));
		return (1.0 + x / 2.0) * (1.0 - x / 2.0) + z * (r / s);
	}
	if ix >= 0x21800000 {
		/* |x| >= 2^-60 */
		x = 0.25 * x * x;
	}
	1.0 - x
}

fn j0_large(ix: u32, x: f32) -> f32 {
	let s = x.sin();
	let c = x.cos();
	let mut cc = s + c;
	if ix < 0x7f000000 {
		let mut ss = s - c;
		let z = -(2.0 * x).cos();
		if s * c < 0.0 {
			cc = z / ss;
		} else {
			ss = z / cc;
		}
		if ix < 0x58800000 {
			cc = pzerof(x) * cc - qzerof(x) * ss;
		}
	}
	INVSQRTPI * cc / x.sqrt()
}

const PR8: [f32; 6] = [
	0.0000000000e+00,
	-7.0312500000e-02,
	-8.0816707611e+00,
	-2.5706311035e+02,
	-2.4852163086e+03,
	-5.2530439453e+03,
];
const PS8: [f32; 5] = [
	1.1653436279e+02,
	3.8337448730e+03,
	4.0597855469e+04,
	1.1675296875e+05,
	4.7627726562e+04,
];
const PR5: [f32; 6] = [
	-1.1412546255e-11,
	-7.0312492549e-02,
	-4.1596107483e+00,
	-6.7674766541e+01,
	-3.3123129272e+02,
	-3.4643338013e+02,
];
const PS5: [f32; 5] = [
	6.0753936768e+01,
	1.0512523193e+03,
	5.9789707031e+03,
	9.6254453125e+03,
	2.4060581055e+03,
];
const PR3: [f32; 6] = [
	-2.5470459075e-09,
	-7.0311963558e-02,
	-2.4090321064e+00,
	-2.1965976715e+01,
	-5.8079170227e+01,
	-3.1447946548e+01,
];
const PS3: [f32; 5] = [
	3.5856033325e+01,
	3.6151397705e+02,
	1.1936077881e+03,
	1.1279968262e+03,
	1.7358093262e+02,
];
const PR2: [f32; 6] = [
	-8.8753431271e-08,
	-7.0303097367e-02,
	-1.4507384300e+00,
	-7.6356959343e+00,
	-1.1193166733e+01,
	-3.2336456776e+00,
];
const PS2: [f32; 5] = [
	2.2220300674e+01,
	1.3620678711e+02,
	2.7047027588e+02,
	1.5387539673e+02,
	1.4657617569e+01,
];

fn pzerof(x: f32) -> f32 {
	let ix = x.to_bits() & 0x7fffffff;
	let (p, q): (&[f32; 6], &[f32; 5]) = if ix >= 0x41000000 {
		(&PR8, &PS8)
	} else if ix >= 0x409173eb {
		(&PR5, &PS5)
	} else if ix >= 0x4036d917 {
		(&PR3, &PS3)
	} else {
		(&PR2, &PS2)
	};
	let z = 1.0 / (x * x);
	let r = p[0] + z * (p[1] + z * (p[2] + z * (p[3] + z * (p[4] + z * p[5]))));
	let s = 1.0 + z * (q[0] + z * (q[1] + z * (q[2] + z * (q[3] + z * q[4]))));
	1.0 + r / s
}

const QR8: [f32; 6] = [
	0.0000000000e+00,
	7.3242187500e-02,
	1.1768206596e+01,
	5.5767340088e+02,
	8.8591972656e+03,
	3.7014625000e+04,
];
const QS8: [f32; 6] = [
	1.6377603149e+02,
	8.0983447266e+03,
	1.4253829688e+05,
	8.0330925000e+05,
	8.4050156250e+05,
	-3.4389928125e+05,
];
const QR5: [f32; 6] = [
	1.8408595828e-11,
	7.3242180049e-02,
	5.8356351852e+00,
	1.3511157227e+02,
	1.0272437744e+03,
	1.9899779053e+03,
];
const QS5: [f32; 6] = [
	8.2776611328e+01,
	2.0778142090e+03,
	1.8847289062e+04,
	5.6751113281e+04,
	3.5976753906e+04,
	-5.3543427734e+03,
];
const QR3: [f32; 6] = [
	4.3774099900e-09,
	7.3241114616e-02,
	3.3442313671e+00,
	4.2621845245e+01,
	1.7080809021e+02,
	1.6673394775e+02,
];
const QS3: [f32; 6] = [
	4.8758872986e+01,
	7.0968920898e+02,
	3.7041481934e+03,
	6.4604252930e+03,
	2.5163337402e+03,
	-1.4924745178e+02,
];
const QR2: [f32; 6] = [
	1.5044444979e-07,
	7.3223426938e-02,
	1.9981917143e+00,
	1.4495602608e+01,
	3.1666231155e+01,
	1.6252708435e+01,
];
const QS2: [f32; 6] = [
	3.0365585327e+01,
	2.6934811401e+02,
	8.4478375244e+02,
	8.8293585205e+02,
	2.1266638184e+02,
	-5.3109550476e+00,
];

fn qzerof(x: f32) -> f32 {
	let ix = x.to_bits() & 0x7fffffff;
	let (p, q): (&[f32; 6], &[f32; 6]) = if ix >= 0x41000000 {
		(&QR8, &QS8)
	} else if ix >= 0x409173eb {
		(&QR5, &QS5)
	} else if ix >= 0x4036d917 {
		(&QR3, &QS3)
	} else {
		(&QR2, &QS2)
	};
	let z = 1.0 / (x * x);
	let r = p[0] + z * (p[1] + z * (p[2] + z * (p[3] + z * (p[4] + z * p[5]))));
	let s = 1.0 + z * (q[0] + z * (q[1] + z * (q[2] + z * (q[3] + z * (q[4] + z * q[5])))));
	(-0.125 + r / s) / x
}

#[cfg(test)]
mod tests {
	use super::j0f;

	fn close(a: f32, b: f32, tol: f32, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	#[test]
	fn j0_at_zero() {
		close(j0f(0.0), 1.0, 0.0, "J0(0)");
		close(j0f(-0.0), 1.0, 0.0, "J0(-0)");
	}

	#[test]
	fn j0_even() {
		for x in [0.3f32, 1.0, 2.5, 8.0, 40.0, 200.0] {
			close(j0f(-x), j0f(x), 0.0, "even");
		}
	}

	#[test]
	fn j0_tabulated() {
		close(j0f(1.0), 0.7651976866, 2e-7, "J0(1)");
		close(j0f(0.5), 0.9384698072, 2e-7, "J0(0.5)");
		close(j0f(2.0), 0.2238907791, 5e-6, "J0(2)");
	}

	#[test]
	fn j0_first_zero() {
		let z0 = 2.4048255577f32;
		close(j0f(z0), 0.0, 3e-4, "first zero");
	}

	#[test]
	fn j0_large_oscillates() {
		let x = 100.0f32;
		let asym = (2.0 / (std::f32::consts::PI * x)).sqrt()
			* (x - std::f32::consts::FRAC_PI_4).cos();
		close(j0f(x), asym, 2e-3, "Hankel |x|=100");
	}
}
