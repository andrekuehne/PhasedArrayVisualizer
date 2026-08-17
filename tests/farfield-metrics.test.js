/**
 * Pattern metrics: SIMD vs scalar WASM, plus synthetic HPBW / SLL checks.
 *
 * Run from the repo root:
 *   node --test tests/farfield-metrics.test.js
 */
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {after, before, describe, test} from 'node:test';
import {linspace} from '../js/util.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const DOMAIN_SPHERICAL = 0;
const DOMAIN_UV = 1;
const DOMAIN_LUDWIG3 = 2;

const METRIC_KEYS = [
	'peak_i1', 'peak_i2', 'peak_ax1', 'peak_ax2',
	'hpbw_ax1', 'hpbw_ax2', 'hpbw_ax1_deg', 'hpbw_ax2_deg',
	'hpbw_ax1_clipped', 'hpbw_ax2_clipped',
	'hpbw_large', 'hpbw_small', 'hpbw_large_deg', 'hpbw_small_deg',
	'hpbw_large_clipped', 'hpbw_small_clipped',
	'hpbw_large_angle_deg', 'hpbw_small_angle_deg',
	'nearest_sll_db', 'largest_sll_db',
	'nearest_sll_ax1', 'nearest_sll_ax2',
	'largest_sll_ax1', 'largest_sll_ax2',
];

async function loadGlue(kind){
	const dir = join(ROOT, 'js', 'wasm', kind);
	const glueUrl = pathToFileURL(join(dir, 'farfield_kernel.js')).href;
	const glue = await import(glueUrl);
	const bytes = await readFile(join(dir, 'farfield_kernel_bg.wasm'));
	await glue.default({module_or_path: bytes});
	return glue;
}

function copyMetrics(m){
	const o = {};
	for (const k of METRIC_KEYS) o[k] = m[k];
	if (typeof m.free === 'function') m.free();
	return o;
}

function f32lin(a, b, n){
	return Float32Array.from(linspace(a, b, n));
}

function addGaussian(total, ax1, ax2, c1, c2, amp, sigma){
	const n1 = ax1.length;
	const n2 = ax2.length;
	const s2 = 2 * sigma * sigma;
	for (let i2 = 0; i2 < n2; i2++){
		for (let i1 = 0; i1 < n1; i1++){
			const d1 = ax1[i1] - c1;
			const d2 = ax2[i2] - c2;
			total[i2 * n1 + i1] += amp * Math.exp(-(d1 * d1 + d2 * d2) / s2);
		}
	}
}

function expectedHpbw(sigma){
	return 2 * sigma * Math.sqrt(2 * Math.LN2);
}

function assertMetricsEqual(label, a, b){
	for (const k of METRIC_KEYS){
		const va = a[k];
		const vb = b[k];
		if (typeof va === 'boolean'){
			assert.equal(va, vb, `${label}: ${k}`);
			continue;
		}
		if (Number.isNaN(va) && Number.isNaN(vb)) continue;
		if (va === vb) continue;
		const scale = Math.max(Math.abs(va), Math.abs(vb), 1);
		assert.ok(
			Math.abs(va - vb) <= 1e-5 * scale,
			`${label}: ${k} SIMD=${va} scalar=${vb}`
		);
	}
}

function runExtract(extract, domain, ax1, ax2, total){
	return copyMetrics(extract(domain, ax1, ax2, total));
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

function ones(n){
	const a = new Float32Array(n);
	a.fill(1);
	return a;
}

describe('WASM pattern metrics', () => {
	let simd;
	let scalar;
	let simdKernel;
	let scalarKernel;

	before(async () => {
		simd = await loadGlue('simd');
		scalar = await loadGlue('scalar');
		simdKernel = new simd.FarfieldKernel();
		scalarKernel = new scalar.FarfieldKernel();
	});

	after(() => {
		simdKernel?.free();
		scalarKernel?.free();
	});

	test('SIMD and scalar agree on a Ludwig3 Gaussian', () => {
		const n = 101;
		const ax1 = f32lin(-0.4, 0.4, n);
		const ax2 = f32lin(-0.4, 0.4, n);
		const total = new Float32Array(n * n);
		addGaussian(total, ax1, ax2, 0, 0, 1, 0.05);
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		assertMetricsEqual('ludwig3 gaussian', a, b);
		const want = expectedHpbw(0.05);
		assert.equal(a.hpbw_ax1_clipped, false);
		assert.equal(a.hpbw_ax2_clipped, false);
		assert.equal(a.hpbw_large_clipped, false);
		assert.equal(a.hpbw_small_clipped, false);
		assert.ok(Math.abs(a.hpbw_ax1 - want) / want < 0.05, `hpbw_ax1 ${a.hpbw_ax1} vs ${want}`);
		assert.ok(Math.abs(a.hpbw_ax2 - want) / want < 0.05, `hpbw_ax2 ${a.hpbw_ax2} vs ${want}`);
		const wantDeg = want * 180 / Math.PI;
		assert.ok(Math.abs(a.hpbw_large_deg - wantDeg) / wantDeg < 0.08);
		assert.ok(Math.abs(a.hpbw_small_deg - wantDeg) / wantDeg < 0.08);
	});

	test('SIMD and scalar agree on two-Gaussian SLL', () => {
		const n = 121;
		const ax1 = f32lin(-0.6, 0.6, n);
		const ax2 = f32lin(-0.6, 0.6, n);
		const total = new Float32Array(n * n);
		addGaussian(total, ax1, ax2, 0, 0, 1, 0.04);
		addGaussian(total, ax1, ax2, 0.18, 0, 0.1, 0.03);
		addGaussian(total, ax1, ax2, 0.35, 0.25, 0.25, 0.03);
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		assertMetricsEqual('two gaussian sll', a, b);
		assert.ok(Math.abs(a.nearest_sll_db - 10 * Math.log10(0.1)) < 0.5);
		assert.ok(Math.abs(a.largest_sll_db - 10 * Math.log10(0.25)) < 0.5);
		assert.ok(a.largest_sll_db > a.nearest_sll_db);
	});

	test('spherical pole uses orthogonal meridian, SIMD ≡ scalar', () => {
		const n = 101;
		const ax1 = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const ax2 = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const total = new Float32Array(n * n);
		const sigma = 0.08;
		for (let i2 = 0; i2 < n; i2++){
			for (let i1 = 0; i1 < n; i1++){
				const th = ax1[i1];
				total[i2 * n + i1] = Math.exp(-(th * th) / (2 * sigma * sigma));
			}
		}
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_SPHERICAL, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_SPHERICAL, ax1, ax2, total);
		assertMetricsEqual('spherical pole', a, b);
		assert.equal(a.hpbw_ax1_clipped, false);
		assert.equal(a.hpbw_ax2_clipped, false);
		const want = expectedHpbw(sigma);
		assert.ok(Math.abs(a.hpbw_ax1 - want) / want < 0.08);
		assert.ok(Math.abs(a.hpbw_ax2 - want) / want < 0.08);
	});

	test('clipped beam is NaN with flags, SIMD ≡ scalar', () => {
		const n = 11;
		const ax1 = f32lin(-0.1, 0.1, n);
		const ax2 = f32lin(-0.1, 0.1, n);
		const total = new Float32Array(n * n);
		total.fill(1);
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_LUDWIG3, ax1, ax2, total);
		assertMetricsEqual('clipped', a, b);
		assert.equal(a.hpbw_ax1_clipped, true);
		assert.equal(a.hpbw_ax2_clipped, true);
		assert.equal(a.hpbw_large_clipped, true);
		assert.equal(a.hpbw_small_clipped, true);
		assert.ok(Number.isNaN(a.hpbw_ax1));
		assert.ok(Number.isNaN(a.hpbw_large_deg));
		assert.ok(Number.isNaN(a.nearest_sll_db));
		assert.ok(Number.isNaN(a.largest_sll_db));
	});

	test('UV invisible samples are ignored for the peak', () => {
		const n = 5;
		const ax1 = f32lin(-1.2, 1.2, n);
		const ax2 = f32lin(-1.2, 1.2, n);
		const total = new Float32Array(n * n);
		for (let i2 = 0; i2 < n; i2++){
			for (let i1 = 0; i1 < n; i1++){
				const u = ax1[i1];
				const v = ax2[i2];
				if (u * u + v * v >= 1) total[i2 * n + i1] = 100;
			}
		}
		const mid = Math.floor(n / 2);
		total[mid * n + mid] = 1;
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_UV, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_UV, ax1, ax2, total);
		assertMetricsEqual('uv invisible', a, b);
		assert.equal(a.peak_i1, mid);
		assert.equal(a.peak_i2, mid);
	});

	test('spherical phi wrap at ±90 deg is not clipped', () => {
		const n = 101;
		const ax1 = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const ax2 = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const sigma = 0.08;
		const th0 = 0.3;
		const ph0 = Math.PI / 2;
		const st = Math.sin(th0);
		const ct = Math.cos(th0);
		const sp = Math.sin(ph0);
		const cp = Math.cos(ph0);
		const r0 = [st * cp, st * sp, ct];
		const total = new Float32Array(n * n);
		for (let i2 = 0; i2 < n; i2++){
			const s2 = Math.sin(ax2[i2]);
			const c2 = Math.cos(ax2[i2]);
			for (let i1 = 0; i1 < n; i1++){
				const s1 = Math.sin(ax1[i1]);
				const c1 = Math.cos(ax1[i1]);
				const r = [s1 * c2, s1 * s2, c1];
				const d = Math.min(1, Math.max(-1, r[0] * r0[0] + r[1] * r0[1] + r[2] * r0[2]));
				const ang = Math.acos(d);
				total[i2 * n + i1] = Math.exp(-(ang * ang) / (2 * sigma * sigma));
			}
		}
		const a = runExtract(simd.extract_pattern_metrics, DOMAIN_SPHERICAL, ax1, ax2, total);
		const b = runExtract(scalar.extract_pattern_metrics, DOMAIN_SPHERICAL, ax1, ax2, total);
		assertMetricsEqual('phi wrap', a, b);
		assert.equal(a.hpbw_ax2_clipped, false);
		assert.ok(Math.abs(a.peak_ax2) > 1.0);
	});

	test('boresight array HPBW agrees across domains', () => {
		const {x, y} = rectArray(8, 8, 0.5, 0.5);
		const mag = ones(x.length);
		const pha = ones(x.length);
		pha.fill(0);
		const n = 129;
		const theta = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const phi = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const u = f32lin(-1, 1, n);
		const v = f32lin(-1, 1, n);
		const az = f32lin(-Math.PI / 2, Math.PI / 2, n);
		const el = f32lin(-Math.PI / 2, Math.PI / 2, n);

		function runDomain(domain, ax1, ax2){
			simdKernel.prepare(ax1.length, ax2.length);
			simdKernel.set_inputs(x, y, mag, pha, ax1, ax2);
			simdKernel.accumulate_tile(domain, 1, 0, ax2.length);
			simdKernel.finalize(x.length);
			const total = simdKernel.take_total();
			return runExtract(simd.extract_pattern_metrics, domain, ax1, ax2, total);
		}

		const sph = runDomain(DOMAIN_SPHERICAL, theta, phi);
		const uv = runDomain(DOMAIN_UV, u, v);
		const l3 = runDomain(DOMAIN_LUDWIG3, az, el);
		assert.equal(sph.hpbw_ax1_clipped, false);
		assert.equal(uv.hpbw_ax1_clipped, false);
		assert.equal(l3.hpbw_ax1_clipped, false);
		const vals = [sph.hpbw_ax1_deg, uv.hpbw_ax1_deg, l3.hpbw_ax1_deg];
		const mean = vals.reduce((s, x) => s + x, 0) / vals.length;
		for (const w of vals){
			assert.ok(
				Math.abs(w - mean) / mean < 0.15,
				`domain HPBW ${vals.join(', ')} mean=${mean}`
			);
		}
	});
});
