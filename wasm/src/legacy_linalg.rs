//! Textbook dense LU / Cholesky used for parity tests and A/B benches.

const PIVOT_FLOOR: f64 = 1e-18;

#[inline]
fn cmul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
	(ar * br - ai * bi, ar * bi + ai * br)
}

#[inline]
fn cdiv(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
	let d = br * br + bi * bi;
	if d < 1e-30 {
		return (0.0, 0.0);
	}
	((ar * br + ai * bi) / d, (ai * br - ar * bi) / d)
}

/// In-place real Cholesky: lower triangle of `a` becomes \(L\) with \(A = L L^T\).
pub fn chol_factor_real(a: &mut [f64], n: usize) {
	for i in 0..n {
		for j in 0..=i {
			let mut sum = a[i * n + j];
			for k in 0..j {
				sum -= a[i * n + k] * a[j * n + k];
			}
			if i == j {
				a[i * n + i] = sum.max(PIVOT_FLOOR).sqrt();
			} else {
				a[i * n + j] = sum / a[j * n + j];
			}
		}
	}
}

pub fn chol_solve_real(l: &[f64], n: usize, b: &mut [f64]) {
	for i in 0..n {
		let mut s = b[i];
		for k in 0..i {
			s -= l[i * n + k] * b[k];
		}
		b[i] = s / l[i * n + i];
	}
	for i in (0..n).rev() {
		let mut s = b[i];
		for k in i + 1..n {
			s -= l[k * n + i] * b[k];
		}
		b[i] = s / l[i * n + i];
	}
}

pub fn chol_factor_herm(re: &mut [f64], im: &mut [f64], n: usize) {
	for i in 0..n {
		for j in 0..=i {
			let mut sum_re = re[i * n + j];
			let mut sum_im = im[i * n + j];
			for k in 0..j {
				let lir = re[i * n + k];
				let lii = im[i * n + k];
				let ljr = re[j * n + k];
				let lji = im[j * n + k];
				sum_re -= lir * ljr + lii * lji;
				sum_im -= lii * ljr - lir * lji;
			}
			if i == j {
				re[i * n + i] = sum_re.max(PIVOT_FLOOR).sqrt();
				im[i * n + i] = 0.0;
			} else {
				let d = re[j * n + j];
				re[i * n + j] = sum_re / d;
				im[i * n + j] = sum_im / d;
			}
		}
	}
}

pub fn chol_solve_herm(l_re: &[f64], l_im: &[f64], n: usize, b_re: &mut [f64], b_im: &mut [f64]) {
	for i in 0..n {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for k in 0..i {
			let lr = l_re[i * n + k];
			let li = l_im[i * n + k];
			sr -= lr * b_re[k] - li * b_im[k];
			si -= lr * b_im[k] + li * b_re[k];
		}
		let d = l_re[i * n + i];
		b_re[i] = sr / d;
		b_im[i] = si / d;
	}
	for i in (0..n).rev() {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for k in i + 1..n {
			let lr = l_re[k * n + i];
			let li = l_im[k * n + i];
			sr -= lr * b_re[k] + li * b_im[k];
			si -= lr * b_im[k] - li * b_re[k];
		}
		let d = l_re[i * n + i];
		b_re[i] = sr / d;
		b_im[i] = si / d;
	}
}

/// Complex LU with partial pivoting. `re`/`im` become \(L+U\) (unit \(L\) diagonal).
/// `piv[k]` is the row swapped with \(k\) at step \(k\).
pub fn lu_factor_cplx(re: &mut [f64], im: &mut [f64], n: usize, piv: &mut [usize]) {
	for i in 0..n {
		piv[i] = i;
	}
	for k in 0..n {
		let mut p = k;
		let mut max_abs = re[k * n + k].hypot(im[k * n + k]);
		for i in (k + 1)..n {
			let a = re[i * n + k].hypot(im[i * n + k]);
			if a > max_abs {
				max_abs = a;
				p = i;
			}
		}
		piv[k] = p;
		if p != k {
			for j in 0..n {
				re.swap(k * n + j, p * n + j);
				im.swap(k * n + j, p * n + j);
			}
		}
		let d_re = re[k * n + k];
		let d_im = im[k * n + k];
		if d_re.hypot(d_im) < PIVOT_FLOOR {
			re[k * n + k] = PIVOT_FLOOR;
			im[k * n + k] = 0.0;
		}
		let pk_re = re[k * n + k];
		let pk_im = im[k * n + k];
		for i in (k + 1)..n {
			let (mr, mi) = cdiv(re[i * n + k], im[i * n + k], pk_re, pk_im);
			re[i * n + k] = mr;
			im[i * n + k] = mi;
			for j in (k + 1)..n {
				let (pr, pi) = cmul(mr, mi, re[k * n + j], im[k * n + j]);
				re[i * n + j] -= pr;
				im[i * n + j] -= pi;
			}
		}
	}
}

pub fn lu_solve_cplx(
	lu_re: &[f64],
	lu_im: &[f64],
	n: usize,
	piv: &[usize],
	b_re: &mut [f64],
	b_im: &mut [f64],
) {
	for k in 0..n {
		let p = piv[k];
		if p != k {
			b_re.swap(k, p);
			b_im.swap(k, p);
		}
	}
	for i in 0..n {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for j in 0..i {
			let (pr, pi) = cmul(lu_re[i * n + j], lu_im[i * n + j], b_re[j], b_im[j]);
			sr -= pr;
			si -= pi;
		}
		b_re[i] = sr;
		b_im[i] = si;
	}
	for i in (0..n).rev() {
		let mut sr = b_re[i];
		let mut si = b_im[i];
		for j in (i + 1)..n {
			let (pr, pi) = cmul(lu_re[i * n + j], lu_im[i * n + j], b_re[j], b_im[j]);
			sr -= pr;
			si -= pi;
		}
		let (xr, xi) = cdiv(sr, si, lu_re[i * n + i], lu_im[i * n + i]);
		b_re[i] = xr;
		b_im[i] = xi;
	}
}

/// Solve \(A X = B\) for complex \(A,B\) (\(n\times n\)), returning \(X\) row-major.
pub fn solve_cplx_multi(
	a_re: &[f64],
	a_im: &[f64],
	n: usize,
	b_re: &[f64],
	b_im: &[f64],
) -> (Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut lu_re = a_re.to_vec();
	let mut lu_im = a_im.to_vec();
	let mut piv = vec![0usize; n];
	lu_factor_cplx(&mut lu_re, &mut lu_im, n, &mut piv);
	let mut x_re = vec![0.0f64; nn];
	let mut x_im = vec![0.0f64; nn];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		for i in 0..n {
			br[i] = b_re[i * n + j];
			bi[i] = b_im[i * n + j];
		}
		lu_solve_cplx(&lu_re, &lu_im, n, &piv, &mut br, &mut bi);
		for i in 0..n {
			x_re[i * n + j] = br[i];
			x_im[i * n + j] = bi[i];
		}
	}
	(x_re, x_im)
}

/// Solve Hermitian \(A X = B\).
pub fn solve_herm_multi(
	a_re: &[f64],
	a_im: &[f64],
	n: usize,
	b_re: &[f64],
	b_im: &[f64],
) -> (Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut l_re = a_re.to_vec();
	let mut l_im = a_im.to_vec();
	chol_factor_herm(&mut l_re, &mut l_im, n);
	let mut x_re = vec![0.0f64; nn];
	let mut x_im = vec![0.0f64; nn];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		for i in 0..n {
			br[i] = b_re[i * n + j];
			bi[i] = b_im[i * n + j];
		}
		chol_solve_herm(&l_re, &l_im, n, &mut br, &mut bi);
		for i in 0..n {
			x_re[i * n + j] = br[i];
			x_im[i * n + j] = bi[i];
		}
	}
	(x_re, x_im)
}

/// Diagonal of \(A^{-1}\) for real SPD `a`.
pub fn inverse_diag_spd(a: &[f64], n: usize) -> Vec<f64> {
	let mut work = a.to_vec();
	chol_factor_real(&mut work, n);
	let mut ydiag = vec![0.0f64; n];
	let mut rhs = vec![0.0f64; n];
	for j in 0..n {
		rhs.fill(0.0);
		rhs[j] = 1.0;
		chol_solve_real(&work, n, &mut rhs);
		ydiag[j] = rhs[j];
	}
	ydiag
}

/// Diagonal of \(A^{-1}\) for complex `a`.
pub fn inverse_diag_cplx(re: &[f64], im: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
	let mut work_re = re.to_vec();
	let mut work_im = im.to_vec();
	let mut piv = vec![0usize; n];
	lu_factor_cplx(&mut work_re, &mut work_im, n, &mut piv);
	let mut y_re = vec![0.0f64; n];
	let mut y_im = vec![0.0f64; n];
	let mut br = vec![0.0f64; n];
	let mut bi = vec![0.0f64; n];
	for j in 0..n {
		br.fill(0.0);
		bi.fill(0.0);
		br[j] = 1.0;
		lu_solve_cplx(&work_re, &work_im, n, &piv, &mut br, &mut bi);
		y_re[j] = br[j];
		y_im[j] = bi[j];
	}
	(y_re, y_im)
}

/// Complex matrix product \(C = A B\), all \(n\times n\) row-major.
#[allow(dead_code)]
pub fn matmul_cplx(
	a_re: &[f64],
	a_im: &[f64],
	b_re: &[f64],
	b_im: &[f64],
	n: usize,
) -> (Vec<f64>, Vec<f64>) {
	let nn = n * n;
	let mut c_re = vec![0.0f64; nn];
	let mut c_im = vec![0.0f64; nn];
	for i in 0..n {
		for j in 0..n {
			let mut re = 0.0;
			let mut im = 0.0;
			for k in 0..n {
				let ar = a_re[i * n + k];
				let ai = a_im[i * n + k];
				let br = b_re[k * n + j];
				let bi = b_im[k * n + j];
				re += ar * br - ai * bi;
				im += ar * bi + ai * br;
			}
			c_re[i * n + j] = re;
			c_im[i * n + j] = im;
		}
	}
	(c_re, c_im)
}

#[cfg(test)]
mod tests {
	use super::*;

	fn close(a: f64, b: f64, tol: f64, label: &str) {
		assert!(
			(a - b).abs() <= tol,
			"{label}: {a} vs {b} (tol {tol})"
		);
	}

	#[test]
	fn chol_solves_spd_identity() {
		let n = 3;
		#[rustfmt::skip]
		let a0 = [
			4.0, 1.0, 0.5,
			1.0, 3.0, 0.2,
			0.5, 0.2, 2.0,
		];
		let mut l = a0.to_vec();
		chol_factor_real(&mut l, n);
		for j in 0..n {
			let mut b = vec![0.0; n];
			b[j] = 1.0;
			chol_solve_real(&l, n, &mut b);
			let mut ax = vec![0.0; n];
			for i in 0..n {
				for k in 0..n {
					ax[i] += a0[i * n + k] * b[k];
				}
			}
			for i in 0..n {
				let expect = if i == j { 1.0 } else { 0.0 };
				close(ax[i], expect, 1e-12, "A x = e_j");
			}
		}
	}

	#[test]
	fn lu_solves_complex() {
		let n = 2;
		#[rustfmt::skip]
		let a_re = [
			1.0, 2.0,
			3.0, 4.0,
		];
		#[rustfmt::skip]
		let a_im = [
			1.0, 0.0,
			0.0, -1.0,
		];
		let mut lu_re = a_re.to_vec();
		let mut lu_im = a_im.to_vec();
		let mut piv = vec![0usize; n];
		lu_factor_cplx(&mut lu_re, &mut lu_im, n, &mut piv);
		let mut br = vec![1.0, 0.0];
		let mut bi = vec![0.0, 0.0];
		lu_solve_cplx(&lu_re, &lu_im, n, &piv, &mut br, &mut bi);
		let mut ax_re = [0.0; 2];
		let mut ax_im = [0.0; 2];
		for i in 0..n {
			for k in 0..n {
				let (pr, pi) = cmul(a_re[i * n + k], a_im[i * n + k], br[k], bi[k]);
				ax_re[i] += pr;
				ax_im[i] += pi;
			}
		}
		close(ax_re[0], 1.0, 1e-12, "lu Ax re0");
		close(ax_re[1], 0.0, 1e-12, "lu Ax re1");
		close(ax_im[0], 0.0, 1e-12, "lu Ax im0");
		close(ax_im[1], 0.0, 1e-12, "lu Ax im1");
	}

	#[test]
	fn herm_chol_agrees_with_real_on_real_matrix() {
		let n = 2;
		#[rustfmt::skip]
		let a = [
			4.0, 1.0,
			1.0, 3.0,
		];
		let mut l_re = a.to_vec();
		let mut l_im = vec![0.0; 4];
		chol_factor_herm(&mut l_re, &mut l_im, n);
		let mut br = vec![1.0, 0.0];
		let mut bi = vec![0.0, 0.0];
		chol_solve_herm(&l_re, &l_im, n, &mut br, &mut bi);
		let mut ax = [0.0; 2];
		for i in 0..2 {
			for k in 0..2 {
				ax[i] += a[i * 2 + k] * br[k];
			}
		}
		close(ax[0], 1.0, 1e-12, "herm Ax0");
		close(ax[1], 0.0, 1e-12, "herm Ax1");
		close(bi[0], 0.0, 1e-12, "herm im0");
		close(bi[1], 0.0, 1e-12, "herm im1");
	}
}
