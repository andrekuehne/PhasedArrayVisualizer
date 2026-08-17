/**
 * Equivalence tests: frozen JS far-field loops vs SIMD and scalar WASM kernels.
 *
 * Run from the repo root:
 *   node --test tests/farfield-equivalence.test.js
 */
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {after, before, describe, test} from 'node:test';
import {linspace} from '../js/util.js';
import {jsLudwig3, jsSpherical, jsUV} from './farfield-js-reference.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const DOMAIN_SPHERICAL = 0;
const DOMAIN_UV = 1;
const DOMAIN_LUDWIG3 = 2;

// SIMD sincos vs JS Math.sin/cos: small ulp noise, not a different formula.
const RTOL = 2e-4;
const ATOL = 5e-5;

async function loadKernel(kind){
	const dir = join(ROOT, 'js', 'wasm', kind);
	const glueUrl = pathToFileURL(join(dir, 'farfield_kernel.js')).href;
	const glue = await import(glueUrl);
	const bytes = await readFile(join(dir, 'farfield_kernel_bg.wasm'));
	await glue.default({module_or_path: bytes});
	return new glue.FarfieldKernel();
}

function runWasm(kernel, domain, frequencyScale, x, y, mag, pha, ax1, ax2, tileRows){
	const n1 = ax1.length;
	const n2 = ax2.length;
	const step = tileRows === undefined ? n2 : tileRows;
	kernel.prepare(n1, n2);
	kernel.set_inputs(x, y, mag, pha, ax1, ax2);
	for (let row0 = 0; row0 < n2; row0 += step){
		kernel.accumulate_tile(domain, frequencyScale, row0, Math.min(step, n2 - row0));
	}
	const maxValue = kernel.finalize(x.length);
	return {maxValue, total: kernel.take_total()};
}

function maxAbsDiff(a, b){
	let m = 0;
	for (let i = 0; i < a.length; i++){
		const d = Math.abs(a[i] - b[i]);
		if (d > m) m = d;
	}
	return m;
}

function assertCloseMesh(label, wasm, js){
	assert.equal(wasm.total.length, js.total.length, `${label}: length`);
	const scale = Math.max(js.maxValue, wasm.maxValue, 1);
	const diff = maxAbsDiff(wasm.total, js.total);
	const maxDiff = Math.abs(wasm.maxValue - js.maxValue);
	assert.ok(
		diff <= ATOL + RTOL * scale,
		`${label}: mesh max|Δ|=${diff} (scale=${scale}, maxJS=${js.maxValue}, maxWASM=${wasm.maxValue})`
	);
	assert.ok(
		maxDiff <= ATOL + RTOL * scale,
		`${label}: maxValue JS=${js.maxValue} WASM=${wasm.maxValue}`
	);
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

function ones(n, value){
	const a = new Float32Array(n);
	a.fill(value === undefined ? 1 : value);
	return a;
}

/** Match PhasedArray.create_farfield_vectors for ideal steering, no illumination. */
function steeringPhase(x, y, thetaDeg, phiDeg){
	const xf = Math.sin(thetaDeg * Math.PI / 180) * Math.cos(phiDeg * Math.PI / 180);
	const yf = Math.sin(thetaDeg * Math.PI / 180) * Math.sin(phiDeg * Math.PI / 180);
	const pha = new Float32Array(x.length);
	const twoPi = 2 * Math.PI;
	for (let i = 0; i < x.length; i++){
		let cycles = (x[i] * xf + y[i] * yf) % 1.0;
		pha[i] = -twoPi * cycles;
	}
	return pha;
}

function taperMag(n){
	const mag = new Float32Array(n);
	for (let i = 0; i < n; i++) mag[i] = 0.4 + 0.6 * (i / Math.max(n - 1, 1));
	return mag;
}

describe('JS vs WASM far-field kernel', () => {
	let simd;
	let scalar;

	before(async () => {
		simd = await loadKernel('simd');
		scalar = await loadKernel('scalar');
	});

	after(() => {
		simd?.free();
		scalar?.free();
	});

	const cases = [
		{
			name: 'spherical remainder grid (n1=5, not multiple of 4)',
			domain: DOMAIN_SPHERICAL,
			freq: 1,
			setup(){
				const {x, y} = rectArray(2, 2, 0.5, 0.5);
				return {
					x, y,
					mag: new Float32Array([1, 0.8, 0.6, 0.4]),
					pha: new Float32Array([0, 0.3, -0.7, 1.1]),
					ax1: linspace(-Math.PI / 2, Math.PI / 2, 5),
					ax2: linspace(-Math.PI / 2, Math.PI / 2, 3),
					js: (a) => jsSpherical(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2, 1),
				};
			},
		},
		{
			name: 'spherical frequencyScale=1.37',
			domain: DOMAIN_SPHERICAL,
			freq: 1.37,
			setup(){
				const {x, y} = rectArray(3, 2, 0.4, 0.6);
				return {
					x, y,
					mag: ones(x.length),
					pha: Float32Array.from(x, (_, i) => 0.15 * i),
					ax1: linspace(-Math.PI / 2, Math.PI / 2, 7),
					ax2: linspace(-Math.PI / 2, Math.PI / 2, 9),
					js: (a) => jsSpherical(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2, 1.37),
				};
			},
		},
		{
			name: 'spherical 8x8 steered 30 deg (app-like)',
			domain: DOMAIN_SPHERICAL,
			freq: 1,
			setup(){
				const {x, y} = rectArray(8, 8, 0.5, 0.5);
				return {
					x, y,
					mag: ones(x.length),
					pha: steeringPhase(x, y, 30, 0),
					ax1: linspace(-Math.PI / 2, Math.PI / 2, 17),
					ax2: linspace(-Math.PI / 2, Math.PI / 2, 17),
					js: (a) => jsSpherical(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2, 1),
				};
			},
		},
		{
			name: 'UV 4-element, frequencyScale must not affect geometry',
			domain: DOMAIN_UV,
			freq: 2.5,
			setup(){
				const x = new Float32Array([0, 0.5, 0, 0.5]);
				const y = new Float32Array([0, 0, 0.5, 0.5]);
				return {
					x, y,
					mag: taperMag(4),
					pha: new Float32Array([0.2, -0.4, 0.8, -1.1]),
					ax1: linspace(-1, 1, 11),
					ax2: linspace(-1, 1, 5),
					js: (a) => jsUV(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2),
				};
			},
		},
		{
			name: 'Ludwig3 az/el with mixed amplitudes',
			domain: DOMAIN_LUDWIG3,
			freq: 0.8,
			setup(){
				const {x, y} = rectArray(3, 3, 0.5, 0.5);
				const sc = Math.PI / 180;
				return {
					x, y,
					mag: taperMag(x.length),
					pha: steeringPhase(x, y, 20, 15),
					ax1: linspace(-90 * sc, 90 * sc, 13),
					ax2: linspace(-90 * sc, 90 * sc, 8),
					js: (a) => jsLudwig3(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2),
				};
			},
		},
		{
			name: 'single isotropic element, 1x1 mesh',
			domain: DOMAIN_SPHERICAL,
			freq: 1,
			setup(){
				const x = new Float32Array([0]);
				const y = new Float32Array([0]);
				const mag = new Float32Array([1]);
				const pha = new Float32Array([0]);
				const ax1 = new Float32Array([0]);
				const ax2 = new Float32Array([0]);
				return {
					x, y, mag, pha, ax1, ax2,
					js: (a) => jsSpherical(a.x, a.y, a.mag, a.pha, a.ax1, a.ax2, 1),
				};
			},
		},
	];

	for (const tc of cases){
		test(tc.name, () => {
			const a = tc.setup();
			const js = a.js(a);
			const simdOut = runWasm(simd, tc.domain, tc.freq, a.x, a.y, a.mag, a.pha, a.ax1, a.ax2);
			const scalarOut = runWasm(scalar, tc.domain, tc.freq, a.x, a.y, a.mag, a.pha, a.ax1, a.ax2);
			assertCloseMesh(`${tc.name} SIMD`, simdOut, js);
			assertCloseMesh(`${tc.name} scalar`, scalarOut, js);
			assertCloseMesh(`${tc.name} SIMD vs scalar`, simdOut, {total: scalarOut.total, maxValue: scalarOut.maxValue});
		});
	}

	test('row tiles match a single full-grid accumulate', () => {
		const {x, y} = rectArray(4, 3, 0.5, 0.5);
		const mag = ones(x.length);
		const pha = steeringPhase(x, y, 10, 25);
		const ax1 = linspace(-Math.PI / 2, Math.PI / 2, 9);
		const ax2 = linspace(-Math.PI / 2, Math.PI / 2, 10);
		const full = runWasm(simd, DOMAIN_SPHERICAL, 1.05, x, y, mag, pha, ax1, ax2, 10);
		const tiled = runWasm(simd, DOMAIN_SPHERICAL, 1.05, x, y, mag, pha, ax1, ax2, 3);
		assert.equal(full.total.length, tiled.total.length);
		assert.equal(maxAbsDiff(full.total, tiled.total), 0);
		assert.equal(full.maxValue, tiled.maxValue);
	});

	test('UV ignores frequencyScale in the geometric kernel', () => {
		const x = new Float32Array([0, 0.7]);
		const y = new Float32Array([0.2, -0.1]);
		const mag = ones(2);
		const pha = new Float32Array([0.5, -0.25]);
		const u = linspace(-1, 1, 6);
		const v = linspace(-1, 1, 4);
		const a = runWasm(simd, DOMAIN_UV, 1.0, x, y, mag, pha, u, v);
		const b = runWasm(simd, DOMAIN_UV, 3.0, x, y, mag, pha, u, v);
		assert.equal(maxAbsDiff(a.total, b.total), 0);
	});
});
