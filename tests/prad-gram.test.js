/**
 * Radiated-power Gram: isolated fields on Gauss-μ × φ → Hermitian P_H.
 *
 * Run from the repo root:
 *   node --test tests/prad-gram.test.js
 */
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {after, before, describe, test} from 'node:test';
import {exponentFromPeakDbi, PATTERN_COS_N, PATTERN_ISOTROPIC} from '../js/phasedarray/element.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const P0 = 0.5;

async function loadKernel(kind){
	const dir = join(ROOT, 'js', 'wasm', kind);
	const glue = await import(pathToFileURL(join(dir, 'farfield_kernel.js')).href);
	const bytes = await readFile(join(dir, 'farfield_kernel_bg.wasm'));
	await glue.default({module_or_path: bytes});
	return new glue.RadiatedPowerKernel();
}

function at(re, im, n, p, q){
	const i = p * n + q;
	return {re: re[i], im: im[i]};
}

function mag(z){
	return Math.hypot(z.re, z.im);
}

function close(a, b, tol, label){
	assert.ok(Math.abs(a - b) <= tol, `${label}: ${a} vs ${b} (tol ${tol})`);
}

function hermitian(re, im, n, tol){
	for (let p = 0; p < n; p++){
		close(im[p * n + p], 0, 1e-6, `diag im[${p}]`);
		for (let q = 0; q < n; q++){
			const a = at(re, im, n, p, q);
			const b = at(re, im, n, q, p);
			close(a.re, b.re, tol, `re Hermitian ${p},${q}`);
			close(a.im, -b.im, tol, `im Hermitian ${p},${q}`);
		}
	}
}

describe('WASM radiated-power Gram', () => {
	/** @type {Record<string, import('../js/wasm/simd/farfield_kernel.js').RadiatedPowerKernel>} */
	const kernels = {};

	before(async () => {
		kernels.simd = await loadKernel('simd');
		kernels.scalar = await loadKernel('scalar');
	});

	after(() => {
		for (const k of Object.values(kernels)){
			if (k && typeof k.free === 'function') k.free();
		}
	});

	for (const kind of ['simd', 'scalar']){
		describe(kind, () => {
			function compute(x, y, freq, elemKind, elemN, nMu, nPhi){
				const k = kernels[kind];
				k.set_quadrature(nMu, nPhi);
				k.compute(x, y, freq, elemKind, elemN);
				return {
					re: k.take_re(),
					im: k.take_im(),
					n: k.n_elements(),
					m: k.n_samples(),
				};
			}

			test('N=1 isotropic radiates P0', () => {
				const {re, im, n, m} = compute(
					new Float32Array([0]),
					new Float32Array([0]),
					1,
					PATTERN_ISOTROPIC,
					0,
					8,
					16
				);
				assert.equal(n, 1);
				assert.equal(m, 8 * 16);
				close(re[0], P0, 2e-4, 'P11');
				close(im[0], 0, 1e-5, 'P11 im');
			});

			test('N=1 cos^n radiates P0', () => {
				const nExp = exponentFromPeakDbi(5);
				const {re, im} = compute(
					new Float32Array([0.1]),
					new Float32Array([-0.2]),
					1,
					PATTERN_COS_N,
					nExp,
					12,
					16
				);
				close(re[0], P0, 5e-4, 'P11 cos^n');
				close(im[0], 0, 1e-5, 'P11 im');
			});

			test('two coincident elements', () => {
				const {re, im, n} = compute(
					new Float32Array([0.3, 0.3]),
					new Float32Array([0.1, 0.1]),
					1,
					PATTERN_ISOTROPIC,
					0,
					8,
					16
				);
				assert.equal(n, 2);
				close(re[0], P0, 2e-4, 'P11');
				close(re[3], P0, 2e-4, 'P22');
				close(re[1], P0, 2e-4, 'P12');
				close(im[1], 0, 2e-4, 'P12 im');
				hermitian(re, im, 2, 1e-5);
			});

			test('large separation off-diagonal is small', () => {
				const {re, im} = compute(
					new Float32Array([0, 20]),
					new Float32Array([0, 0]),
					1,
					PATTERN_ISOTROPIC,
					0,
					48,
					96
				);
				close(re[0], P0, 1e-3, 'P11');
				assert.ok(mag(at(re, im, 2, 0, 1)) < 0.05, `|P12|=${mag(at(re, im, 2, 0, 1))}`);
			});

			test('Hermitian and nearly real for planar axisymmetric elements', () => {
				const {re, im, n} = compute(
					new Float32Array([0, 0.5, 1]),
					new Float32Array([0, 0.25, -0.25]),
					1.2,
					PATTERN_COS_N,
					1,
					16,
					32
				);
				assert.equal(n, 3);
				hermitian(re, im, 3, 1e-5);
				for (let i = 0; i < re.length; i++){
					close(im[i], 0, 3e-4, `im[${i}]`);
				}
			});

			test('quadrature convergence on a short baseline', () => {
				const x = new Float32Array([0, 0.5]);
				const y = new Float32Array([0, 0]);
				const coarse = compute(x, y, 1, PATTERN_ISOTROPIC, 0, 12, 24);
				const fine = compute(x, y, 1, PATTERN_ISOTROPIC, 0, 24, 48);
				close(fine.re[0], P0, 1e-4, 'fine P11');
				assert.ok(Math.abs(coarse.re[0] - fine.re[0]) < 5e-4, 'P11 Δ');
				assert.ok(Math.abs(coarse.re[1] - fine.re[1]) < 5e-3, 'P12 Δ');
			});
		});
	}

	test('SIMD and scalar Gram agree', () => {
		const x = new Float32Array([0, 0.5, 1.0, 0.25]);
		const y = new Float32Array([0, 0.25, -0.25, 0.4]);
		for (const k of Object.values(kernels)){
			k.set_quadrature(16, 24);
			k.compute(x, y, 1.1, PATTERN_COS_N, 0.5);
		}
		const reS = kernels.simd.take_re();
		const imS = kernels.simd.take_im();
		const reC = kernels.scalar.take_re();
		const imC = kernels.scalar.take_im();
		assert.equal(reS.length, reC.length);
		let max = 0;
		for (let i = 0; i < reS.length; i++){
			max = Math.max(max, Math.abs(reS[i] - reC[i]), Math.abs(imS[i] - imC[i]));
		}
		assert.ok(max < 2e-4, `simd vs scalar max|Δ|=${max}`);
	});
});

describe('sample-panel split', () => {
	test('sequential panels sum to the full Gram', async () => {
		const {splitRowRanges, mergeGrams} = await import('../js/wasm/farfield-pool.js');
		const k = await loadKernel('simd');
		try {
			const x = new Float32Array([0, 0.5, 1.0, 0.25]);
			const y = new Float32Array([0, 0.25, -0.25, 0.4]);
			k.set_quadrature(8, 12);
			k.compute(x, y, 1.05, PATTERN_ISOTROPIC, 0);
			const n = k.n_elements();
			const m = k.n_samples();
			const fullRe = k.take_re();
			const fullIm = k.take_im();
			const ranges = splitRowRanges(m, 3);
			const tiles = ranges.map((r) => {
				k.fill_isolated_range(x, y, 1.05, PATTERN_ISOTROPIC, 0, r.row0, r.rowCount);
				k.form_gram();
				return {re: k.take_re(), im: k.take_im()};
			});
			const merged = mergeGrams(n, tiles);
			let max = 0;
			for (let i = 0; i < fullRe.length; i++){
				max = Math.max(max, Math.abs(merged.re[i] - fullRe[i]), Math.abs(merged.im[i] - fullIm[i]));
			}
			assert.ok(max < 3e-5, `panel merge max|Δ|=${max}`);
		}
		finally {
			k.free();
		}
	});

	test('worker panels match the main-thread Gram', async () => {
		const {availableParallelism} = await import('node:os');
		const {Worker} = await import('node:worker_threads');
		const {
			startFarfieldPool,
			stopFarfieldPool,
			runPradJob,
			farfieldPoolSize,
		} = await import('../js/wasm/farfield-pool.js');
		const wasmPath = join(ROOT, 'js', 'wasm', 'simd', 'farfield_kernel_bg.wasm');
		stopFarfieldPool();
		const started = await startFarfieldPool({
			simd: true,
			Worker,
			workers: Math.max(2, availableParallelism()),
			wasmPath,
		});
		if ((started.workers || farfieldPoolSize()) < 2){
			stopFarfieldPool();
			assert.fail('worker pool did not start');
			return;
		}
		try {
			const k = await loadKernel('simd');
			const x = new Float32Array([0, 0.5, 1.0, 0.25, 0.8]);
			const y = new Float32Array([0, 0.25, -0.25, 0.4, -0.1]);
			k.set_quadrature(8, 12);
			k.compute(x, y, 1.05, PATTERN_COS_N, 0.7);
			const fullRe = k.take_re();
			const fullIm = k.take_im();
			k.free();
			const merged = await runPradJob({
				x, y,
				frequencyScale: 1.05,
				elementKind: PATTERN_COS_N,
				elementN: 0.7,
				nMu: 8,
				nPhi: 12,
			});
			let max = 0;
			for (let i = 0; i < fullRe.length; i++){
				max = Math.max(
					max,
					Math.abs(merged.re[i] - fullRe[i]),
					Math.abs(merged.im[i] - fullIm[i])
				);
			}
			assert.ok(max < 2e-4, `worker merge max|Δ|=${max}`);
		}
		finally {
			stopFarfieldPool();
		}
	});
});
