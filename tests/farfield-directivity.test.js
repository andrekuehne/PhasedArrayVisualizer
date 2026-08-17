/**
 * Directivity quadrature: each domain integrates |AF|² dΩ over the front hemisphere.
 *
 * Run from the repo root:
 *   node --test tests/farfield-directivity.test.js
 */
import assert from 'node:assert/strict';
import {describe, test} from 'node:test';
import {FarfieldLudwig3, FarfieldSpherical, FarfieldUV} from '../js/phasedarray/farfield.js';

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
