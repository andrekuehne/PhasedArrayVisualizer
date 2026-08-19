/**
 * Matched S-matrix from the J0 radiated-power Gram (§6–8), including optional
 * \(jX(\Delta x,\Delta y)+jX_\mathrm{self}I\) at a real \(z_0\).
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
import {alignPhaseCycles, conjugatePhaseCycles, gemv, identityT, nMuFromGeometry} from '../js/phasedarray/matched.js';

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
				k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, 0, 0);
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
				const pRe = k.take_re();
				const zRe = k.take_z_re();
				const zIm = k.take_z_im();
				assert.equal(zRe.length, 1);
				assert.equal(zIm.length, 1);
				close(zRe[0], 2 * Z_REF * pRe[0], 1e-4, 'Z re = 2 Zref P_H');
				close(zIm[0], 0, 1e-8, 'Z im');
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
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, 0, 0);
		const z0 = k.take_z0();
		const sRe = k.take_s_re();
		const sIm = k.take_s_im();
		assert.equal(k.n_elements(), n);
		assert.equal(z0.length, n);
		assert.equal(sRe.length, n * n);
		const pRe = k.take_re();
		const zRe = k.take_z_re();
		const zIm = k.take_z_im();
		assert.equal(zRe.length, n * n);
		assert.equal(zIm.length, n * n);
		let maxZim = 0;
		let maxZerr = 0;
		for (let i = 0; i < n * n; i++){
			maxZim = Math.max(maxZim, Math.abs(zIm[i]));
			maxZerr = Math.max(maxZerr, Math.abs(zRe[i] - 2 * Z_REF * pRe[i]));
		}
		assert.ok(maxZim < 1e-8, `max |Zim|=${maxZim}`);
		assert.ok(maxZerr < 1e-3, `Z vs 2 Zref P_H max err ${maxZerr}`);
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
		k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, 0, 0);
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
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, 0, 0);
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

	test('alignPhaseCycles removes a global offset and anchors conjugate to geometric', () => {
		const n = 8;
		const target = new Float32Array(n);
		const source = new Float32Array(n);
		const offset = 0.37;
		for (let i = 0; i < n; i++){
			target[i] = 0.15 * i;
			source[i] = target[i] + offset;
		}
		const aligned = alignPhaseCycles(source, target);
		let maxErr = 0;
		for (let i = 0; i < n; i++) maxErr = Math.max(maxErr, Math.abs(aligned[i] - target[i]));
		assert.ok(maxErr < 1e-6, `constant offset residual ${maxErr}`);

		const k = kernels.simd;
		const nx = 8;
		const ny = 8;
		const {x, y} = rectArray(nx, ny, 0.5, 0.5);
		const nn = x.length;
		k.set_quadrature(32, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, 0, 0);
		const tRe = k.take_t_re();
		const tIm = k.take_t_im();
		const theta = 30;
		const xf = Math.sin(theta * Math.PI / 180);
		const cx = 0.5 * 0.5 * (nx - 1);
		const geo = new Float32Array(nn);
		for (let i = 0; i < nn; i++) geo[i] = (x[i] + cx) * xf;

		const iso = alignPhaseCycles(
			conjugatePhaseCycles(x, y, theta, 0, 1, null, null, PATTERN_ISOTROPIC, 0),
			geo
		);
		const emb = alignPhaseCycles(
			conjugatePhaseCycles(x, y, theta, 0, 1, tRe, tIm, PATTERN_ISOTROPIC, 0),
			geo
		);
		let maxIso = 0;
		let meanRe = 0;
		let meanIm = 0;
		const twoPi = 2 * Math.PI;
		for (let i = 0; i < nn; i++){
			maxIso = Math.max(maxIso, Math.abs(iso[i] - geo[i]));
			const d = twoPi * (emb[i] - geo[i]);
			meanRe += Math.cos(d);
			meanIm += Math.sin(d);
		}
		assert.ok(maxIso < 0.02, `aligned isolated vs geometric ${maxIso}`);
		assert.ok(Math.hypot(meanRe, meanIm) / nn > 0.99, 'aligned embedded offset is ~0');
		let maxEmbGeo = 0;
		for (let i = 0; i < nn; i++) maxEmbGeo = Math.max(maxEmbGeo, Math.abs(emb[i] - geo[i]));
		assert.ok(maxEmbGeo > 1e-4, `aligned embedded still differs from geometric ${maxEmbGeo}`);
	});

	test('three irregular points with Xnn: real z0, leftover Sii', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 0.5, 1.0]);
		const y = new Float32Array([0, 0.25, -0.25]);
		const n = x.length;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, 10, 2, 0, 0, 0, 0);
		const sRe = k.take_s_re();
		const sIm = k.take_s_im();
		const z0Im = k.take_z0_im();
		const zIm = k.take_z_im();
		assert.equal(sRe.length, n * n);
		assert.ok(k.match_residual() < TAU, `residual ${k.match_residual()}`);
		for (let i = 0; i < n; i++){
			close(z0Im[i], 0, 0, `z0_im[${i}]`);
			close(zIm[i * n + i], 0, 0, `X${i}${i}`);
		}
		let maxSii = 0;
		for (let i = 0; i < n; i++){
			maxSii = Math.max(maxSii, mag(sRe[i * n + i], sIm[i * n + i]));
		}
		assert.ok(maxSii > 1e-3, `|Sii| leftover reactance ${maxSii}`);
		let maxIm = 0;
		for (let i = 0; i < n * n; i++) maxIm = Math.max(maxIm, Math.abs(sIm[i]));
		assert.ok(maxIm > 1e-6, `S should be complex, max|Im|=${maxIm}`);
	});

	test('right-angle Xnn with Location A splits equal-distance pairs', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 1, 0]);
		const y = new Float32Array([0, 0, 1]);
		const n = 3;
		const xnn = 10;
		const aniso = 0.5;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, xnn, 2, 0, aniso, 0, 0);
		const zIm = k.take_z_im();
		close(zIm[1], xnn * (1 + aniso), 1e-6, 'X along +x');
		close(zIm[2], xnn * (1 - aniso), 1e-6, 'X along +y');
		close(zIm[1], zIm[n], 1e-12, 'X01 = X10');
		close(zIm[2], zIm[2 * n], 1e-12, 'X02 = X20');
		assert.ok(k.match_residual() < TAU, `residual ${k.match_residual()}`);
	});

	test('collinear Xnn with Oscillation β=π flips next-nearest sign', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 0.5, 1.0]);
		const y = new Float32Array([0, 0, 0]);
		const xnn = 10;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, xnn, 2, Math.PI, 0, 0, 0);
		const zIm = k.take_z_im();
		close(zIm[1], xnn, 1e-6, 'X01 nn');
		close(zIm[2], -xnn * 0.25, 1e-6, 'X02 opposite');
		assert.ok(zIm[1] * zIm[2] < 0, 'next-nearest opposite sign');
	});

	test('N=1 common Zref: S = 0, T = 1', () => {
		const k = kernels.simd;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, Z_REF, 0);
		close(k.take_z0()[0], Z_REF, 1e-5, 'z0');
		close(k.take_z0_im()[0], 0, 1e-12, 'z0 im');
		close(k.take_s_re()[0], 0, 1e-8, 'S11 re');
		close(k.take_t_re()[0], 1, 1e-6, 'T11');
		assert.equal(k.match_iterations(), 0);
	});

	test('N=1 common real zc: S11 = (R-zc)/(R+zc)', () => {
		const k = kernels.simd;
		const zc = 40;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, zc, 0);
		const r = k.take_z_re()[0];
		const s = (r - zc) / (r + zc);
		close(k.take_s_re()[0], s, 1e-8, 'S11');
		close(k.take_s_im()[0], 0, 1e-12, 'S11 im');
		close(k.take_z0()[0], zc, 0, 'z0');
	});

	test('N=1 common Self X on Z: S = (R+jX-zc)/(R+jX+zc)', () => {
		const k = kernels.simd;
		const zc = 50;
		const xSelf = 10;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, zc, xSelf);
		const zRe = k.take_z_re()[0];
		const zIm = k.take_z_im()[0];
		const numRe = zRe - zc;
		const numIm = zIm;
		const denRe = zRe + zc;
		const denIm = zIm;
		const d2 = denRe * denRe + denIm * denIm;
		const sRe = (numRe * denRe + numIm * denIm) / d2;
		const sIm = (numIm * denRe - numRe * denIm) / d2;
		close(zIm, xSelf, 0, 'X11');
		close(k.take_s_re()[0], sRe, 1e-8, 'S11 re');
		close(k.take_s_im()[0], sIm, 1e-8, 'S11 im');
		close(k.take_z0()[0], zc, 0, 'z0 re');
		close(k.take_z0_im()[0], 0, 0, 'z0 im');
	});

	test('N=1 per-port Self X leaves S11 = jX/(2R+jX)', () => {
		const k = kernels.simd;
		const xSelf = 10;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 0, 0, 0, xSelf);
		const r = k.take_z_re()[0];
		const denRe = 2 * r;
		const denIm = xSelf;
		const d2 = denRe * denRe + denIm * denIm;
		const sRe = (xSelf * denIm) / d2;
		const sIm = (xSelf * denRe) / d2;
		close(k.take_z0()[0], r, 1e-5, 'z0 = R');
		close(k.take_z0_im()[0], 0, 0, 'z0 im');
		close(k.take_z_im()[0], xSelf, 0, 'X11');
		close(k.take_s_re()[0], sRe, 1e-8, 'S11 re');
		close(k.take_s_im()[0], sIm, 1e-8, 'S11 im');
		assert.ok(mag(k.take_s_re()[0], k.take_s_im()[0]) > 1e-3, 'S11 from leftover X');
		assert.ok(k.match_residual() < TAU, `residual ${k.match_residual()}`);
	});

	test('8×8 common Z0 is flat and |Sii| exceeds per-port', () => {
		const k = kernels.simd;
		const nx = 8;
		const ny = 8;
		const {x, y} = rectArray(nx, ny, 0.5, 0.5);
		const n = x.length;
		k.set_quadrature(32, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, 0, 0);
		const sPer = k.take_s_re();
		const sImPer = k.take_s_im();
		let maxPer = 0;
		for (let i = 0; i < n; i++) maxPer = Math.max(maxPer, mag(sPer[i * n + i], sImPer[i * n + i]));
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, Z_REF, 0);
		const z0 = k.take_z0();
		const z0Im = k.take_z0_im();
		const sRe = k.take_s_re();
		const sIm = k.take_s_im();
		assert.equal(k.match_iterations(), 0);
		for (let i = 0; i < n; i++){
			close(z0[i], Z_REF, 0, `z0[${i}]`);
			close(z0Im[i], 0, 0, `z0_im[${i}]`);
		}
		let maxSii = 0;
		for (let i = 0; i < n; i++) maxSii = Math.max(maxSii, mag(sRe[i * n + i], sIm[i * n + i]));
		assert.ok(maxSii > 0.01, `common max |Sii|=${maxSii}`);
		assert.ok(maxSii > maxPer * 10, `common ${maxSii} vs per-port ${maxPer}`);
	});

	test('propagation: closest pair is Xnn, next-nearest flips, R unchanged', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 0.5, 1.0]);
		const y = new Float32Array([0, 0, 0]);
		const xnn = 10;
		const zc = 45;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s(Z_REF, x, y, 0, 0, 0, 0, zc, 0);
		const r01 = k.take_z_re()[1];
		k.form_matched_s_propagation(Z_REF, x, y, xnn, 0, 1, 1, 1, zc, 0);
		const zIm = k.take_z_im();
		const zRe = k.take_z_re();
		const z0 = k.take_z0();
		close(zIm[1], xnn, 1e-6, 'X01 nn');
		close(zIm[3], xnn, 1e-6, 'X10');
		close(zIm[2], -xnn, 1e-6, 'X02 flip');
		close(zRe[1], r01, 1e-6, 'R01 gram');
		close(z0[0], zc, 0, 'z0');
		close(z0[1], zc, 0, 'z0 1');
		close(k.take_z0_im()[0], 0, 0, 'z0 im');
	});

	test('propagation: εx≠εy splits equal-distance x/y pairs', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 0.5, 0]);
		const y = new Float32Array([0, 0, 0.5]);
		const n = 3;
		const xnn = 8;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s_propagation(Z_REF, x, y, xnn, 0, 1, 4, 1, Z_REF, 0);
		const zIm = k.take_z_im();
		close(zIm[1], xnn, 1e-6, 'X01 ref');
		assert.ok(Math.abs(zIm[2] - zIm[1]) > 1e-6, 'εy splits y-arm');
		close(zIm[2], zIm[2 * n], 1e-12, 'X02 = X20');
	});

	test('propagation: α_λ shrinks farther pairs', () => {
		const k = kernels.simd;
		const x = new Float32Array([0, 0.5, 1.0]);
		const y = new Float32Array([0, 0, 0]);
		const xnn = 8;
		const att = 2;
		k.set_quadrature(24, 2);
		k.compute_j0(x, y, 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s_propagation(Z_REF, x, y, xnn, att, 1, 1, 1, Z_REF, 0);
		const zIm = k.take_z_im();
		close(zIm[1], xnn, 1e-6, 'nn');
		const far = xnn * Math.exp(-att * (1.0 - 0.5)) * (-1);
		close(zIm[2], far, 1e-6, 'far');
		assert.ok(Math.abs(zIm[2]) < Math.abs(zIm[1]), 'farther weaker');
	});

	test('propagation: invalid zc clamps to Z_REF', () => {
		const k = kernels.simd;
		k.set_quadrature(8, 2);
		k.compute_j0(new Float32Array([0]), new Float32Array([0]), 1, PATTERN_ISOTROPIC, 0);
		k.form_matched_s_propagation(Z_REF, new Float32Array([0]), new Float32Array([0]), 0, 0, 1, 1, 1, 0, 0);
		close(k.take_z0()[0], Z_REF, 0, 'z0 clamp');
	});
});
