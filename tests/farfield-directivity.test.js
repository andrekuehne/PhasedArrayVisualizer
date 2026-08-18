/**
 * Directivity quadrature: each domain integrates |AF|² dΩ over the front hemisphere.
 *
 * Run from the repo root:
 *   node --test tests/farfield-directivity.test.js
 */
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {before, describe, test} from 'node:test';
import {exponentFromPeakDbi, PATTERN_COS_N, PATTERN_ISOTROPIC} from '../js/phasedarray/element.js';
import {FarfieldLudwig3, FarfieldSpherical, FarfieldUV} from '../js/phasedarray/farfield.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const DOMAIN_SPHERICAL = 0;
const DOMAIN_UV = 1;
const DOMAIN_LUDWIG3 = 2;

function fill(ff, fn){
	const [p1, p2] = ff.meshPoints;
	let maxValue = -Infinity;
	for (let i2 = 0; i2 < p2; i2++){
		for (let i1 = 0; i1 < p1; i1++){
			const v = fn(i1, i2);
			ff.farfield_total[i2][i1] = v;
			if (v > maxValue) maxValue = v;
		}
	}
	ff.maxValue = maxValue;
}

function flatten(ff){
	const [p1, p2] = ff.meshPoints;
	const flat = new Float32Array(p1 * p2);
	for (let i2 = 0; i2 < p2; i2++){
		flat.set(ff.farfield_total[i2], i2 * p1);
	}
	return flat;
}

async function loadGlue(){
	const dir = join(ROOT, 'js', 'wasm', 'scalar');
	const glue = await import(pathToFileURL(join(dir, 'farfield_kernel.js')).href);
	const bytes = await readFile(join(dir, 'farfield_kernel_bg.wasm'));
	await glue.default({module_or_path: bytes});
	return glue;
}

describe('domain directivity integrals', () => {
	test('uniform pattern over the front hemisphere approaches D=2', () => {
		const n = 129;
		const sph = new FarfieldSpherical(n, n, 1);
		const uv = new FarfieldUV(n, n, 1, 1);
		const l3 = new FarfieldLudwig3(n, n, 1);
		fill(sph, () => 1);
		fill(uv, () => 1);
		fill(l3, () => 1);
		assert.ok(Math.abs(sph.compute_directivity() - 2) < 0.05);
		assert.ok(Math.abs(uv.compute_directivity() - 2) < 0.08);
		assert.ok(Math.abs(l3.compute_directivity() - 2) < 0.03);
	});

	test('cosine-element pattern approaches D=4', () => {
		const n = 129;
		const sph = new FarfieldSpherical(n, n, 1);
		const uv = new FarfieldUV(n, n, 1, 1);
		const l3 = new FarfieldLudwig3(n, n, 1);
		fill(sph, (it) => Math.abs(Math.cos(sph.theta[it])));
		fill(uv, (iu, iv) => {
			const w2 = 1 - uv.u[iu] ** 2 - uv.v[iv] ** 2;
			return w2 > 0 ? Math.sqrt(w2) : 0;
		});
		fill(l3, (ia, ie) => Math.abs(Math.cos(l3.az[ia]) * Math.cos(l3.el[ie])));
		assert.ok(Math.abs(sph.compute_directivity() - 4) < 0.08);
		assert.ok(Math.abs(uv.compute_directivity() - 4) < 0.08);
		assert.ok(Math.abs(l3.compute_directivity() - 4) < 0.05);
	});

	test('UV samples outside the unit circle do not contribute', () => {
		const inside = new FarfieldUV(65, 65, 1, 1);
		const outside = new FarfieldUV(65, 65, 1, 1.5);
		fill(inside, () => 1);
		fill(outside, () => 1);
		const dIn = inside.compute_directivity();
		const dOut = outside.compute_directivity();
		assert.ok(Math.abs(dIn - dOut) / dIn < 0.04);
	});
});

describe('element exponent', () => {
	test('matches 10^(G/10)/2 - 1', () => {
		assert.equal(exponentFromPeakDbi(3), 0);
		const n5 = exponentFromPeakDbi(5);
		assert.ok(Math.abs(n5 - (Math.pow(10, 0.5) / 2 - 1)) < 1e-12);
		assert.ok(Math.abs(n5 - 0.58113883) < 1e-6);
	});
});

describe('WASM element pattern apply', () => {
	/** @type {Awaited<ReturnType<typeof loadGlue>>} */
	let glue;

	before(async () => {
		glue = await loadGlue();
	});

	function applyAndDirectivity(ff, domain, ax1, ax2, kind, nExp){
		fill(ff, () => 1);
		const flat = flatten(ff);
		const peak = glue.apply_element_pattern(
			domain,
			Float32Array.from(ax1),
			Float32Array.from(ax2),
			flat,
			kind,
			nExp
		);
		ff.wrap_flat_total(flat);
		ff.maxValue = peak;
		return ff.compute_directivity();
	}

	test('WASM exponent matches JS', () => {
		assert.ok(Math.abs(glue.element_exponent_from_peak_dbi(5) - exponentFromPeakDbi(5)) < 1e-6);
		assert.equal(glue.element_exponent_from_peak_dbi(3), 0);
	});

	test('isotropic apply on a uniform pattern still gives D≈2', () => {
		const n = 129;
		const sph = new FarfieldSpherical(n, n, 1);
		const uv = new FarfieldUV(n, n, 1, 1);
		const l3 = new FarfieldLudwig3(n, n, 1);
		assert.ok(Math.abs(applyAndDirectivity(sph, DOMAIN_SPHERICAL, sph.theta, sph.phi, PATTERN_ISOTROPIC, 0.8) - 2) < 0.05);
		assert.ok(Math.abs(applyAndDirectivity(uv, DOMAIN_UV, uv.u, uv.v, PATTERN_ISOTROPIC, 0.8) - 2) < 0.08);
		assert.ok(Math.abs(applyAndDirectivity(l3, DOMAIN_LUDWIG3, l3.az, l3.el, PATTERN_ISOTROPIC, 0.8) - 2) < 0.03);
	});

	test('cos^n zeros UV samples outside the unit circle even at n=0', () => {
		const uv = new FarfieldUV(5, 5, 1, 1.5);
		fill(uv, () => 1);
		const flat = flatten(uv);
		glue.apply_element_pattern(
			DOMAIN_UV,
			Float32Array.from(uv.u),
			Float32Array.from(uv.v),
			flat,
			PATTERN_COS_N,
			0
		);
		const [p1] = uv.meshPoints;
		for (let iv = 0; iv < uv.vPoints; iv++){
			for (let iu = 0; iu < uv.uPoints; iu++){
				const r2 = uv.u[iu] ** 2 + uv.v[iv] ** 2;
				const val = flat[iv * p1 + iu];
				if (r2 >= 1) assert.equal(val, 0);
				else assert.ok(Math.abs(val - 1) < 1e-6);
			}
		}
	});

	test('cos^n at 5 dBi approaches D≈3.16 on all domains', () => {
		const nPts = 129;
		const nExp = exponentFromPeakDbi(5);
		const expected = Math.pow(10, 0.5);
		const sph = new FarfieldSpherical(nPts, nPts, 1);
		const uv = new FarfieldUV(nPts, nPts, 1, 1);
		const l3 = new FarfieldLudwig3(nPts, nPts, 1);
		assert.ok(Math.abs(applyAndDirectivity(sph, DOMAIN_SPHERICAL, sph.theta, sph.phi, PATTERN_COS_N, nExp) - expected) < 0.08);
		assert.ok(Math.abs(applyAndDirectivity(uv, DOMAIN_UV, uv.u, uv.v, PATTERN_COS_N, nExp) - expected) < 0.08);
		assert.ok(Math.abs(applyAndDirectivity(l3, DOMAIN_LUDWIG3, l3.az, l3.el, PATTERN_COS_N, nExp) - expected) < 0.08);
	});
});
