/**
 * Timing probe for the radiated-power Gram: main thread vs sample-panel workers.
 *
 * Run from the repo root:
 *   node --test tests/prad-bench.test.js
 */
import {availableParallelism} from 'node:os';
import {readFile} from 'node:fs/promises';
import {dirname, join} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {Worker} from 'node:worker_threads';
import {after, before, describe, test} from 'node:test';
import {PATTERN_ISOTROPIC} from '../js/phasedarray/element.js';
import {
	farfieldPoolSize,
	pradWorkerCount,
	runPradJob,
	startFarfieldPool,
	stopFarfieldPool,
} from '../js/wasm/farfield-pool.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');

const CASES = [
	{nx: 8, ny: 8, nMu: 16, nPhi: 32},
	{nx: 8, ny: 8, nMu: 32, nPhi: 64},
	{nx: 16, ny: 16, nMu: 16, nPhi: 32},
	{nx: 16, ny: 16, nMu: 32, nPhi: 64},
	{nx: 32, ny: 32, nMu: 16, nPhi: 32},
	{nx: 32, ny: 32, nMu: 32, nPhi: 64},
	{nx: 64, ny: 64, nMu: 16, nPhi: 32},
	{nx: 64, ny: 64, nMu: 32, nPhi: 64},
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
		padL('nμ', 4),
		padL('nφ', 4),
		padL('M', 6),
		padL('W', 3),
		padL('fill ms', 9),
		padL('gram ms', 9),
		padL('total ms', 9),
	].join('  ');
	const lines = [`\n${title}`, header, '-'.repeat(header.length), ...rows, ''];
	console.log(lines.join('\n'));
}

describe('radiated-power Gram timings', () => {
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

	test('main thread fill-A vs Gram', () => {
		const rows = [];
		for (const build of ['simd', 'scalar']){
			const k = kernels[build];
			for (const c of CASES){
				if (build === 'scalar' && c.nx === 64 && c.nMu === 32) continue;
				const {x, y} = rectArray(c.nx, c.ny, 0.5, 0.5);
				k.set_quadrature(c.nMu, c.nPhi);
				k.fill_isolated(x, y, 1, PATTERN_ISOTROPIC, 0);
				k.form_gram();
				const t0 = performance.now();
				k.fill_isolated(x, y, 1, PATTERN_ISOTROPIC, 0);
				const t1 = performance.now();
				k.form_gram();
				const t2 = performance.now();
				const fill = t1 - t0;
				const gram = t2 - t1;
				rows.push([
					pad(build, 8),
					pad(`${c.nx}x${c.ny}`, 8),
					padL(k.n_elements(), 5),
					padL(c.nMu, 4),
					padL(c.nPhi, 4),
					padL(k.n_samples(), 6),
					padL(1, 3),
					padL(fill.toFixed(1), 9),
					padL(gram.toFixed(1), 9),
					padL((fill + gram).toFixed(1), 9),
				].join('  '));
			}
		}
		printTable('Main thread', rows);
	});
});

describe('radiated-power Gram worker timings', () => {
	const hw = Math.max(2, availableParallelism());

	after(() => {
		stopFarfieldPool();
	});

	async function benchBuild(build){
		stopFarfieldPool();
		const wasmPath = join(ROOT, 'js', 'wasm', build, 'farfield_kernel_bg.wasm');
		const started = await startFarfieldPool({
			simd: build === 'simd',
			Worker,
			workers: hw,
			wasmPath,
		});
		const wCount = started.workers || farfieldPoolSize();
		if (wCount < 2){
			console.log(`\nWorkers (${build}): pool did not start (got ${wCount})\n`);
			return;
		}
		const rows = [];
		for (const c of CASES){
			const {x, y} = rectArray(c.nx, c.ny, 0.5, 0.5);
			const n = x.length;
			const m = c.nMu * c.nPhi;
			const wUse = pradWorkerCount(n, m, wCount);
			const spec = {
				x, y,
				frequencyScale: 1,
				elementKind: PATTERN_ISOTROPIC,
				elementN: 0,
				nMu: c.nMu,
				nPhi: c.nPhi,
			};
			await runPradJob(spec);
			const t0 = performance.now();
			await runPradJob(spec);
			const total = performance.now() - t0;
			rows.push([
				pad(build, 8),
				pad(`${c.nx}x${c.ny}`, 8),
				padL(n, 5),
				padL(c.nMu, 4),
				padL(c.nPhi, 4),
				padL(m, 6),
				padL(wUse, 3),
				padL('—', 9),
				padL('—', 9),
				padL(total.toFixed(1), 9),
			].join('  '));
		}
		printTable(`Workers (${wCount} started, sample-panel split)`, rows);
		stopFarfieldPool();
	}

	test('worker sample-panel Gram', async () => {
		await benchBuild('simd');
		await benchBuild('scalar');
	});
});
