/**
 * Timing probe for Green PEC-dipole Z: unique-lag fill vs from_z LU.
 *
 * Run from the repo root:
 *   node --test tests/green-bench.test.js
 */
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {after, before, describe, test} from 'node:test';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const Z_REF = 50;
const H = 0.25;
const ELL = 0.1;
const A = 0.001;

const CASES = [
	{nx: 8, ny: 8},
	{nx: 16, ny: 16},
	{nx: 32, ny: 32},
];

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

function pad(s, w){
	s = String(s);
	return s.length >= w ? s : s + ' '.repeat(w - s.length);
}

function padL(s, w){
	s = String(s);
	return s.length >= w ? s : ' '.repeat(w - s.length) + s;
}

function printTable(title, rows){
	const header = [
		pad('build', 8),
		pad('array', 8),
		padL('N', 5),
		padL('fill ms', 9),
		padL('LU ms', 9),
		padL('total ms', 9),
	].join('  ');
	const lines = [`\n${title}`, header, '-'.repeat(header.length), ...rows, ''];
	console.log(lines.join('\n'));
}

describe('Green PEC unique-lag Z timings', () => {
	const kernels = {};

	before(async () => {
		kernels.simd = await loadKernel('simd');
	});

	after(() => {
		for (const k of Object.values(kernels)){
			if (k && typeof k.free === 'function') k.free();
		}
	});

	test('unique-lag fill vs from_z LU', {timeout: 120_000}, () => {
		const k = kernels.simd;
		const rows = [];
		for (const c of CASES){
			const {x, y} = rectArray(c.nx, c.ny, 0.5, 0.5);
			k.fill_green_pec_dipole_z(x, y, 1, H, ELL, A);
			k.form_from_z(Z_REF, Z_REF, 0);
			const t0 = performance.now();
			k.fill_green_pec_dipole_z(x, y, 1, H, ELL, A);
			const t1 = performance.now();
			k.form_from_z(Z_REF, Z_REF, 0);
			const t2 = performance.now();
			const fill = t1 - t0;
			const lu = t2 - t1;
			rows.push([
				pad('simd', 8),
				pad(`${c.nx}x${c.ny}`, 8),
				padL(k.n_elements(), 5),
				padL(fill.toFixed(1), 9),
				padL(lu.toFixed(1), 9),
				padL((fill + lu).toFixed(1), 9),
			].join('  '));
		}
		printTable('Green PEC (unique-lag Z + from_z)', rows);
	});
});
