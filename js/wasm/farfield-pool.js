/**
 * Persistent worker pool for the far-field WASM kernel and radiated-power Gram.
 * Splits ax2 rows or quadrature samples across workers; no SharedArrayBuffer.
 */

const MAX_WORKERS = 8;
const INIT_TIMEOUT_MS = 30000;
/** Cap in-flight partial Gram copies (re+im f32) at about 512 MiB. */
const PRAD_GRAM_MEM_BUDGET = 512 * 1024 * 1024;

/** @typedef {{row0: number, rowCount: number}} RowRange */
/** @typedef {{total: Float32Array, maxValue: number}} TotalTile */
/** @typedef {{re: Float32Array, im: Float32Array}} GramTile */

/** @type {null | {
 *   workers: Worker[],
 *   readyCount: number,
 *   jobId: number,
 *   cancelCurrent: null | (() => void),
 *   handlers: Array<null | ((msg: object) => void)>,
 * }} */
let pool = null;

export function workerCount(override){
	if (override > 0) return Math.max(1, Math.min(override | 0, MAX_WORKERS));
	const hw = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
	return Math.max(1, Math.min(hw, MAX_WORKERS));
}

export function farfieldPoolSize(){
	return pool === null ? 0 : pool.readyCount;
}

/**
 * Contiguous ranges covering `0..n2`.
 * @param {number} n2
 * @param {number} nWorkers
 * @returns {RowRange[]}
 */
export function splitRowRanges(n2, nWorkers){
	const n = Math.max(1, Math.min(Math.max(0, nWorkers | 0), n2));
	if (n2 <= 0) return [];
	const base = Math.floor(n2 / n);
	const rem = n2 % n;
	const ranges = [];
	let row0 = 0;
	for (let i = 0; i < n; i++){
		const rowCount = base + (i < rem ? 1 : 0);
		if (rowCount <= 0) continue;
		ranges.push({row0, rowCount});
		row0 += rowCount;
	}
	return ranges;
}

/**
 * Concatenate per-worker intensity tiles into a full row-major buffer.
 * @param {number} n1
 * @param {RowRange[]} ranges
 * @param {TotalTile[]} tiles
 * @returns {TotalTile}
 */
export function mergeTotals(n1, ranges, tiles){
	let n2 = 0;
	for (let i = 0; i < ranges.length; i++) n2 += ranges[i].rowCount;
	const total = new Float32Array(n1 * n2);
	let maxValue = -Infinity;
	for (let i = 0; i < ranges.length; i++){
		const {row0, rowCount} = ranges[i];
		const tile = tiles[i];
		if (!tile || !tile.total) continue;
		const expect = rowCount * n1;
		if (tile.total.length !== expect){
			throw new Error(`mergeTotals: tile ${i} length ${tile.total.length} != ${expect}`);
		}
		total.set(tile.total, row0 * n1);
		if (tile.maxValue > maxValue) maxValue = tile.maxValue;
	}
	return {total, maxValue};
}

/**
 * Sum partial Hermitian Grams (sample-panel split).
 * @param {number} n
 * @param {GramTile[]} tiles
 * @returns {GramTile}
 */
export function mergeGrams(n, tiles){
	const nn = n * n;
	const re = new Float32Array(nn);
	const im = new Float32Array(nn);
	for (let t = 0; t < tiles.length; t++){
		const tile = tiles[t];
		if (!tile || !tile.re || !tile.im) continue;
		if (tile.re.length !== nn || tile.im.length !== nn){
			throw new Error(`mergeGrams: tile ${t} length ${tile.re.length} != ${nn}`);
		}
		for (let i = 0; i < nn; i++){
			re[i] += tile.re[i];
			im[i] += tile.im[i];
		}
	}
	return {re, im};
}

/**
 * @param {number} n
 * @param {Float32Array} re
 * @param {Float32Array} im
 * @param {GramTile} tile
 */
function addGramTile(n, re, im, tile){
	const nn = n * n;
	if (!tile || !tile.re || !tile.im) return;
	if (tile.re.length !== nn || tile.im.length !== nn){
		throw new Error(`addGramTile: length ${tile.re.length} != ${nn}`);
	}
	for (let i = 0; i < nn; i++){
		re[i] += tile.re[i];
		im[i] += tile.im[i];
	}
}

/**
 * How many workers to use for a Gram of size n so partial P copies stay in budget.
 * @param {number} n
 * @param {number} m
 * @param {number} available
 */
export function pradWorkerCount(n, m, available){
	const bytes = 8 * n * n;
	const maxByMem = bytes <= 0 ? available : Math.max(1, Math.floor(PRAD_GRAM_MEM_BUDGET / bytes));
	return Math.max(1, Math.min(available, m, maxByMem));
}

/**
 * @param {{
 *   simd?: boolean,
 *   workers?: number,
 *   Worker?: typeof Worker,
 *   workerUrl?: URL | string,
 *   wasmPath?: string,
 * }} [opts]
 * @returns {Promise<{workers: number}>}
 */
export async function startFarfieldPool(opts){
	if (pool !== null) return {workers: pool.readyCount};
	opts = opts || {};
	const WorkerCtor = opts.Worker || (typeof Worker !== 'undefined' ? Worker : null);
	if (!WorkerCtor) return {workers: 0};
	const n = workerCount(opts.workers);
	if (n < 2) return {workers: 0};

	const simd = Boolean(opts.simd);
	let url = opts.workerUrl;
	if (!url){
		try {
			url = new URL('./farfield-worker.js', import.meta.url);
		}
		catch {
			return {workers: 0};
		}
	}

	const started = await Promise.all(
		Array.from({length: n}, () => spawnWorker(WorkerCtor, url, simd, opts.wasmPath))
	);
	const workers = started.filter(Boolean);
	if (workers.length < 2){
		for (const w of workers) w.terminate();
		return {workers: 0};
	}

	pool = {
		workers,
		readyCount: workers.length,
		jobId: 0,
		cancelCurrent: null,
		handlers: new Array(workers.length).fill(null),
	};
	for (let i = 0; i < workers.length; i++){
		listen(workers[i], (data) => {
			const h = pool && pool.handlers[i];
			if (h) h(data);
		}, () => {});
	}
	return {workers: pool.readyCount};
}

export function stopFarfieldPool(){
	if (pool === null) return;
	if (typeof pool.cancelCurrent === 'function') pool.cancelCurrent();
	for (const w of pool.workers){
		try { w.terminate(); } catch { /* ignore */ }
	}
	pool = null;
}

/**
 * @param {Worker} w
 * @param {(data: object) => void} onData
 * @param {(err: object) => void} onError
 */
function listen(w, onData, onError){
	if (typeof w.on === 'function'){
		w.on('message', onData);
		w.on('error', onError);
		return;
	}
	w.onmessage = (ev) => onData(ev.data);
	w.onerror = onError;
}

/**
 * @param {typeof Worker} WorkerCtor
 * @param {URL | string} url
 * @param {boolean} simd
 * @param {string} [wasmPath]
 * @returns {Promise<Worker | null>}
 */
function spawnWorker(WorkerCtor, url, simd, wasmPath){
	return new Promise((resolve) => {
		let w;
		try {
			w = new WorkerCtor(url, {type: 'module'});
		}
		catch {
			resolve(null);
			return;
		}
		let settled = false;
		const finish = (worker) => {
			if (settled) return;
			settled = true;
			clearTimeout(timer);
			resolve(worker);
		};
		const timer = setTimeout(() => {
			try { w.terminate(); } catch { /* ignore */ }
			finish(null);
		}, INIT_TIMEOUT_MS);
		listen(w, (msg) => {
			if (msg && msg.type === 'ready'){
				finish(w);
				return;
			}
			if (msg && msg.type === 'error' && msg.id === undefined){
				try { w.terminate(); } catch { /* ignore */ }
				finish(null);
			}
		}, () => {
			try { w.terminate(); } catch { /* ignore */ }
			finish(null);
		});
		w.postMessage({type: 'init', simd, wasmPath: wasmPath || undefined});
	});
}

/**
 * @param {{
 *   domain: number,
 *   frequencyScale: number,
 *   nElements: number,
 *   x: Float32Array,
 *   y: Float32Array,
 *   mag: Float32Array,
 *   pha: Float32Array,
 *   ax1: Float32Array,
 *   ax2: Float32Array,
 *   tileRows: number,
 * }} spec
 * @param {(done: number, total: number) => void} [onProgress]
 * @returns {Promise<TotalTile>}
 */
export function runFarfieldJob(spec, onProgress){
	if (pool === null || pool.readyCount < 2){
		return Promise.reject(new Error('Farfield worker pool is not available.'));
	}
	if (typeof pool.cancelCurrent === 'function') pool.cancelCurrent();

	const jobId = ++pool.jobId;
	pool.handlers.fill(null);
	const ax2 = spec.ax2;
	const n1 = spec.ax1.length;
	const n2 = ax2.length;
	const nUse = Math.min(pool.readyCount, n2);
	const ranges = splitRowRanges(n2, nUse);

	return new Promise((resolve, reject) => {
		let settled = false;
		const settle = (fn, value) => {
			if (settled) return;
			settled = true;
			if (pool !== null && pool.jobId === jobId) pool.cancelCurrent = null;
			fn(value);
		};
		pool.cancelCurrent = () => settle(reject, new Error('cancelled'));

		/** @type {(TotalTile | null)[]} */
		const tiles = new Array(ranges.length).fill(null);
		const progressDone = new Array(ranges.length).fill(0);
		const progressMax = new Array(ranges.length).fill(1);
		let remaining = ranges.length;

		const reportProgress = () => {
			if (typeof onProgress !== 'function') return;
			let done = 0;
			let total = 0;
			for (let i = 0; i < ranges.length; i++){
				done += progressDone[i];
				total += progressMax[i];
			}
			onProgress(done, total);
		};

		for (let i = 0; i < ranges.length; i++){
			const worker = pool.workers[i];
			const range = ranges[i];
			pool.handlers[i] = (msg) => {
				if (!msg || msg.id !== jobId) return;
				if (msg.type === 'progress'){
					progressDone[i] = msg.done;
					progressMax[i] = msg.total;
					reportProgress();
					return;
				}
				if (msg.type === 'result'){
					tiles[i] = {total: msg.total, maxValue: msg.maxValue};
					remaining--;
					if (remaining === 0){
						try {
							settle(resolve, mergeTotals(n1, ranges, tiles));
						}
						catch (e){
							settle(reject, e);
						}
					}
					return;
				}
				if (msg.type === 'error'){
					settle(reject, new Error(msg.message || 'Farfield worker error'));
				}
			};
			const ax2Slice = ax2.slice(range.row0, range.row0 + range.rowCount);
			worker.postMessage({
				type: 'run',
				id: jobId,
				domain: spec.domain,
				frequencyScale: spec.frequencyScale,
				nElements: spec.nElements,
				x: spec.x,
				y: spec.y,
				mag: spec.mag,
				pha: spec.pha,
				ax1: spec.ax1,
				ax2: ax2Slice,
				tileRows: spec.tileRows,
			});
		}
	});
}

/**
 * Split the Gauss-μ × φ sample axis across workers and sum partial Grams.
 * @param {{
 *   x: Float32Array,
 *   y: Float32Array,
 *   frequencyScale: number,
 *   elementKind: number,
 *   elementN: number,
 *   nMu: number,
 *   nPhi: number,
 * }} spec
 * @param {(done: number, total: number) => void} [onProgress]
 * @returns {Promise<GramTile>}
 */
export function runPradJob(spec, onProgress){
	if (pool === null || pool.readyCount < 2){
		return Promise.reject(new Error('Farfield worker pool is not available.'));
	}
	if (typeof pool.cancelCurrent === 'function') pool.cancelCurrent();

	const jobId = ++pool.jobId;
	pool.handlers.fill(null);
	const n = spec.x.length;
	const m = Math.max(1, (spec.nMu | 0) * (spec.nPhi | 0));
	const nUse = pradWorkerCount(n, m, pool.readyCount);
	const ranges = splitRowRanges(m, nUse);

	return new Promise((resolve, reject) => {
		let settled = false;
		const settle = (fn, value) => {
			if (settled) return;
			settled = true;
			if (pool !== null && pool.jobId === jobId) pool.cancelCurrent = null;
			fn(value);
		};
		pool.cancelCurrent = () => settle(reject, new Error('cancelled'));

		const re = new Float32Array(n * n);
		const im = new Float32Array(n * n);
		const progressDone = new Array(ranges.length).fill(0);
		const progressMax = new Array(ranges.length).fill(1);
		let remaining = ranges.length;

		const reportProgress = () => {
			if (typeof onProgress !== 'function') return;
			let done = 0;
			let total = 0;
			for (let i = 0; i < ranges.length; i++){
				done += progressDone[i];
				total += progressMax[i];
			}
			onProgress(done, total);
		};

		for (let i = 0; i < ranges.length; i++){
			const worker = pool.workers[i];
			const range = ranges[i];
			pool.handlers[i] = (msg) => {
				if (!msg || msg.id !== jobId) return;
				if (msg.type === 'progress'){
					progressDone[i] = msg.done;
					progressMax[i] = msg.total;
					reportProgress();
					return;
				}
				if (msg.type === 'result'){
					try {
						addGramTile(n, re, im, msg);
					}
					catch (e){
						settle(reject, e);
						return;
					}
					remaining--;
					if (remaining === 0) settle(resolve, {re, im});
					return;
				}
				if (msg.type === 'error'){
					settle(reject, new Error(msg.message || 'Radiated-power worker error'));
				}
			};
			worker.postMessage({
				type: 'run_prad',
				id: jobId,
				x: spec.x,
				y: spec.y,
				frequencyScale: spec.frequencyScale,
				elementKind: spec.elementKind,
				elementN: spec.elementN,
				nMu: spec.nMu,
				nPhi: spec.nPhi,
				sample0: range.row0,
				sampleCount: range.rowCount,
			});
		}
	});
}
