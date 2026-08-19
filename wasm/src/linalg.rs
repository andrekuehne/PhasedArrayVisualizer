//! faer-backed dense solvers for matched \(S\)/\(T\).
//!
//! Inputs are row-major split real/imag (`index = row * n + col`). All calls
//! use sequential parallelism so the same path works in WASM.

use faer::linalg::solvers::Solve;
use faer::{c64, Mat, Par, Side};

const PIVOT_FLOOR: f64 = 1e-18;

fn ensure_seq() {
	static INIT: std::sync::Once = std::sync::Once::new();
	INIT.call_once(|| {
		faer::set_global_parallelism(Par::Seq);
	});
}

fn pack_cplx(re: &[f64], im: &[f64], n: usize) -> Mat<c64> {
	Mat::from_fn(n, n, |i, j| {
		let k = i * n + j;
		c64::new(re[k], im[k])
	})
}

fn pack_real(a: &[f64], n: usize) -> Mat<f64> {
	Mat::from_fn(n, n, |i, j| a[i * n + j])
}

fn unpack_cplx(x: &Mat<c64>, n: usize) -> (Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut re = vec![0.0f64; nn];
	let mut im = vec![0.0f64; nn];
	for i in 0..n {
		for j in 0..n {
			let z = x[(i, j)];
			re[i * n + j] = z.re;
			im[i * n + j] = z.im;
		}
	}
	(re, im)
}

/// Solve \(A X = B\) for complex \(A,B\) (\(n\times n\)), returning \(X\) row-major.
pub fn solve_cplx_multi(
	a_re: &[f64],
	a_im: &[f64],
	n: usize,
	b_re: &[f64],
	b_im: &[f64],
) -> (Vec<f64>, Vec<f64>) {
	ensure_seq();
	if n == 0 {
		return (Vec::new(), Vec::new());
	}
	let a = pack_cplx(a_re, a_im, n);
	let b = pack_cplx(b_re, b_im, n);
	let x = a.partial_piv_lu().solve(&b);
	unpack_cplx(&x, n)
}

/// Solve Hermitian \(A X = B\). Falls back to LU if \(LL^H\) fails.
pub fn solve_herm_multi(
	a_re: &[f64],
	a_im: &[f64],
	n: usize,
	b_re: &[f64],
	b_im: &[f64],
) -> (Vec<f64>, Vec<f64>) {
	ensure_seq();
	if n == 0 {
		return (Vec::new(), Vec::new());
	}
	let a = pack_cplx(a_re, a_im, n);
	let b = pack_cplx(b_re, b_im, n);
	if let Ok(llt) = a.llt(Side::Lower) {
		return unpack_cplx(&llt.solve(&b), n);
	}
	let a_reg = Mat::from_fn(n, n, |i, j| {
		let z = a[(i, j)];
		if i == j {
			z + c64::new(PIVOT_FLOOR, 0.0)
		} else {
			z
		}
	});
	if let Ok(llt) = a_reg.llt(Side::Lower) {
		return unpack_cplx(&llt.solve(&b), n);
	}
	solve_cplx_multi(a_re, a_im, n, b_re, b_im)
}

/// Diagonal of \(A^{-1}\) for real SPD `a`.
pub fn inverse_diag_spd(a: &[f64], n: usize) -> Vec<f64> {
	ensure_seq();
	if n == 0 {
		return Vec::new();
	}
	let mat = pack_real(a, n);
	let b = Mat::<f64>::identity(n, n);
	let x = if let Ok(llt) = mat.llt(Side::Lower) {
		llt.solve(&b)
	} else {
		let mat_reg = Mat::from_fn(n, n, |i, j| {
			if i == j {
				mat[(i, j)] + PIVOT_FLOOR
			} else {
				mat[(i, j)]
			}
		});
		mat_reg
			.llt(Side::Lower)
			.expect("SPD Cholesky failed after regularization")
			.solve(&b)
	};
	let mut ydiag = vec![0.0f64; n];
	for j in 0..n {
		ydiag[j] = x[(j, j)];
	}
	ydiag
}

/// Diagonal of \(A^{-1}\) for complex `a`.
pub fn inverse_diag_cplx(re: &[f64], im: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut b_re = vec![0.0f64; nn];
	let b_im = vec![0.0f64; nn];
	for j in 0..n {
		b_re[j * n + j] = 1.0;
	}
	let (x_re, x_im) = solve_cplx_multi(re, im, n, &b_re, &b_im);
	let mut y_re = vec![0.0f64; n];
	let mut y_im = vec![0.0f64; n];
	for j in 0..n {
		y_re[j] = x_re[j * n + j];
		y_im[j] = x_im[j * n + j];
	}
	(y_re, y_im)
}

/// Complex matrix product \(C = A B\), all \(n\times n\) row-major.
pub fn matmul_cplx(
	a_re: &[f64],
	a_im: &[f64],
	b_re: &[f64],
	b_im: &[f64],
	n: usize,
) -> (Vec<f64>, Vec<f64>) {
	ensure_seq();
	if n == 0 {
		return (Vec::new(), Vec::new());
	}
	let a = pack_cplx(a_re, a_im, n);
	let b = pack_cplx(b_re, b_im, n);
	let c = &a * &b;
	unpack_cplx(&c, n)
}
