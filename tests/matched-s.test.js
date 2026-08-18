/**
 * Matched S-matrix from the J0 radiated-power Gram (§6–8).
 *
 * Run from the repo root:
 *   node --test tests/matched-s.test.js
 */
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {after, before, describe, test} from 'node:test';
import {PATTERN_ISOTROPIC} from '../js/phasedarray/element.js';
import {conjugatePhaseCycles, gemv, identityT, nMuFromGeometry} from '../js/phasedarray/matched.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const Z_REF = 50;
const TAU = 1e-3;

async function loadKernel(kind){
	const dir = join(ROOT, 'js', 'wasm', kind);
	const glue = await import(pathToFileURL(join(dir, 'farfield_kernel.js')).href);
	const bytes = await readFile(join(dir, 'farfield_kernel_bg.wasm'));
	await glue.default({module_or_path: bytes});
	return new glue.RadiatedPowerKernel();
}

function rectArray(nx, ny, dx, dy){
	const n = nx * ny;
	const x = new Float32Array(n);
	const y = new Float32Array(n);
	let k = 0;
	for (let ix = 0; ix < nx; ix++){
		for (let iy = 0; iy < ny; iy++, k++){
			x[k] = dx * ix;
			y[k] = dy * iy;
		}
	}
	return {x, y};
}

function mag(re, im){
	return Math.hypot(re, im);
}

function close(a, b, tol, label){
	assert.ok(Math.abs(a - b) <= tol, `${label}: ${a} vs ${b} (tol ${tol})`);
}

function db20(v){
	if (!(v > 0)) return -Infinity;
	return 20 * Math.log10(v);
}

function fmt(v, digits){
	if (!Number.isFinite(v)) return String(v);
	const a = Math.abs(v);
	if (a !== 0 && (a >= 1e3 || a < 1e-3)) return v.toExponential(digits);
	return v.toFixed(digits);
}

/**
 * S is N×N with N = nx*ny. Diagonal table is laid out as the array:
 * index = ix * ny + iy.
 */
function printSiiGrid(sRe, nx, ny){
	const n = nx * ny;
	const colW = 11;
	const header = [' iy\\ix', ...Array.from({length: nx}, (_, ix) => String(ix).padStart(colW))].join(' ');
	const lines = ['', `Sii (real) on ${nx}×${ny} lattice:`, header];
	for (let iy = 0; iy < ny; iy++){
		const row = [String(iy).padStart(6)];
		for (let ix = 0; ix < nx; ix++){
			const i = ix * ny + iy;
			row.push(fmt(sRe[i * n + i], 4).padStart(colW));
		}
		lines.push(row.join(' '));
	}
	return lines;
}

describe('matched S from J0 Prad', () => {
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
			test('N=1 isotropic: z0 = Zref, S = 0', () => {
				const k = kernels[kind];
				k.set_quadrature(8, 2);
				k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
				k.form_matched_s(Z_REF);
				const z0 = k.take_z0();
				const sRe = k.take_s_re();
				const sIm = k.take_s_im();
				assert.equal(k.n_elements(), 1);
				close(z0[0], Z_REF, 1e-5, 'z0');
				close(sRe[0], 0, 1e-8, 'S11 re');
				close(sIm[0], 0, 1e-12, 'S11 im');
				const tRe = k.take_t_re();
				const tIm = k.take_t_im();
				close(tRe[0], 1, 1e-6, 'T11');
				close(tIm[0], 0, 1e-12, 'T11 im');
				assert.ok(k.match_residual() < TAU, `residual ${k.match_residual()}`);
			});
		});
	}

	test('8×8 isotropic λ/2: Sii, z0, worst |Sij|', () => {
		const k = kernels.simd;
		const nx = 8;
		const ny = 8;
		const {x, y} = rectArray(nx, ny, 0.5, 0.5);
		const n = x.length;
		k.set_quadrature(32, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF);
		const z0 = k.take_z0();
		const sRe = k.take_s_re();
		const sIm = k.take_s_im();
		assert.equal(k.n_elements(), n);
		assert.equal(z0.length, n);
		assert.equal(sRe.length, n * n);
		assert.ok(k.match_residual() < TAU, `residual ${k.match_residual()}`);

		let maxSii = 0;
		let zMin = Infinity;
		let zMax = -Infinity;
		let zSum = 0;
		for (let i = 0; i < n; i++){
			zMin = Math.min(zMin, z0[i]);
			zMax = Math.max(zMax, z0[i]);
			zSum += z0[i];
			maxSii = Math.max(maxSii, mag(sRe[i * n + i], sIm[i * n + i]));
		}
		assert.ok(maxSii < 1e-3, `max |Sii|=${maxSii}`);

		let maxSij = 0;
		let worstI = 0;
		let worstJ = 1;
		let maxIm = 0;
		let maxAsym = 0;
		for (let i = 0; i < n; i++){
			for (let j = 0; j < n; j++){
				const a = i * n + j;
				const b = j * n + i;
				maxIm = Math.max(maxIm, Math.abs(sIm[a]));
				maxAsym = Math.max(maxAsym, Math.abs(sRe[a] - sRe[b]), Math.abs(sIm[a] - sIm[b]));
				if (i === j) continue;
				const m = mag(sRe[a], sIm[a]);
				if (m > maxSij){
					maxSij = m;
					worstI = i;
					worstJ = j;
				}
			}
		}
		assert.ok(maxIm < 1e-6, `max |Sim|=${maxIm}`);
		assert.ok(maxAsym < 1e-6, `S not symmetric ${maxAsym}`);

		const ixW = Math.floor(worstI / ny);
		const iyW = worstI % ny;
		const jxW = Math.floor(worstJ / ny);
		const jyW = worstJ % ny;
		const lines = printSiiGrid(sRe, nx, ny);
		lines.push('');
		lines.push(`match iterations ${k.match_iterations()}  residual ${k.match_residual().toExponential(3)}`);
		lines.push(`z0 Ω  min ${fmt(zMin, 4)}  max ${fmt(zMax, 4)}  mean ${fmt(zSum / n, 4)}`);
		lines.push(`max |Sii| ${maxSii.toExponential(3)}  (${db20(maxSii).toFixed(1)} dB)`);
		lines.push(
			`worst |Sij| ${maxSij.toExponential(4)}  (${db20(maxSij).toFixed(2)} dB)`
			+ `  i=${worstI} (ix=${ixW},iy=${iyW})  j=${worstJ} (ix=${jxW},iy=${jyW})`
		);
		lines.push('');
		console.log(lines.join('\n'));
	});

	test('N=1 T is identity and GEMV is a no-op', () => {
		const k = kernels.simd;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF);
		const tRe = k.take_t_re();
		const tIm = k.take_t_im();
		close(tRe[0], 1, 1e-6, 'T11');
		const w = gemv(tRe, tIm, [0.4], [-0.3]);
		close(w.re[0], 0.4, 1e-6, 'w re');
		close(w.im[0], -0.3, 1e-6, 'w im');
	});

	test('8×8 GEMV w = T a and conjugate phases differ from geometric', () => {
		const k = kernels.simd;
		const nx = 8;
		const ny = 8;
		const {x, y} = rectArray(nx, ny, 0.5, 0.5);
		const n = x.length;
		k.set_quadrature(32, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF);
		const tRe = k.take_t_re();
		const tIm = k.take_t_im();
		assert.equal(tRe.length, n * n);
		const aRe = new Float64Array(n);
		const aIm = new Float64Array(n);
		for (let i = 0; i < n; i++) aRe[i] = 1;
		const w = gemv(tRe, tIm, aRe, aIm);
		assert.equal(w.re.length, n);
		let maxOff = 0;
		for (let i = 0; i < n; i++){
			if (i !== 0) maxOff = Math.max(maxOff, Math.abs(tRe[i]));
		}
		assert.ok(maxOff > 1e-3, `T has off-diagonal ${maxOff}`);
		const ident = identityT(n);
		const wI = gemv(ident.re, ident.im, aRe, aIm);
		let maxDiff = 0;
		for (let i = 0; i < n; i++) maxDiff = Math.max(maxDiff, Math.abs(w.re[i] - wI.re[i]));
		assert.ok(maxDiff > 1e-3, `T a differs from a (${maxDiff})`);

		const theta = 30;
		const phi = 0;
		const cyclesEmb = conjugatePhaseCycles(x, y, theta, phi, 1, tRe, tIm, PATTERN_ISOTROPIC, 0);
		const cyclesIso = conjugatePhaseCycles(x, y, theta, phi, 1, null, null, PATTERN_ISOTROPIC, 0);
		const xf = Math.sin(theta * Math.PI / 180);
		const cx = 0.5 * 0.5 * (nx - 1);
		let maxGeoIso = 0;
		let maxEmbIso = 0;
		const d0 = cyclesIso[0] - (x[0] + cx) * xf;
		for (let i = 0; i < n; i++){
			const geo = (x[i] + cx) * xf;
			let d = cyclesIso[i] - geo - d0;
			d -= Math.round(d);
			maxGeoIso = Math.max(maxGeoIso, Math.abs(d));
			let e = cyclesEmb[i] - cyclesIso[i];
			e -= Math.round(e);
			maxEmbIso = Math.max(maxEmbIso, Math.abs(e));
		}
		assert.ok(maxGeoIso < 0.02, `isolated conjugate vs geometric ${maxGeoIso}`);
		assert.ok(maxEmbIso > 1e-4, `embedded conjugate differs from isolated ${maxEmbIso}`);
		assert.ok(nMuFromGeometry({r: [1.75]}, 1) >= 17);
	});
});
