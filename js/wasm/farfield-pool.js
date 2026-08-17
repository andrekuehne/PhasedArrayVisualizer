/**
 * Persistent Dedicated Worker pool for the far-field WASM kernel.
 * Splits ax2 rows across workers; no SharedArrayBuffer.
 */

const MAX_WORKERS = 8;
const INIT_TIMEOUT_MS = 30000;

/** @typedef {{row0: number, rowCount: number}} RowRange */
/** @typedef {{total: Float32Array, maxValue: number}} TotalTile */

/** @type {null | {
 *   workers: Worker[],
 *   readyCount: number,
 *   jobId: number,
 *   cancelCurrent: null | (() => void),
 *   handlers: Array<null | ((msg: object) => void)>,
 * }} */
let pool = null;

export function workerCount(){
	const hw = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
	return Math.max(1, Math.min(hw, MAX_WORKERS));
}

export function farfieldPoolSize(){
	return pool === null ? 0 : pool.readyCount;
}

/**
 * Contiguous ax2 row ranges covering `0..n2`.
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
 * @param {{simd: boolean}} opts
 * @returns {Promise<{workers: number}>}
 */
export async function startFarfieldPool(opts){
	if (pool !== null) return {workers: pool.readyCount};
	if (typeof Worker === 'undefined') return {workers: 0};
	const n = workerCount();
	if (n < 2) return {workers: 0};

	const simd = Boolean(opts && opts.simd);
	let url;
	try {
		url = new URL('./farfield-worker.js', import.meta.url);
	}
	catch {
		return {workers: 0};
	}

	const started = await Promise.all(Array.from({length: n}, () => spawnWorker(url, simd)));
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
		workers[i].onerror = null;
		workers[i].onmessage = (ev) => {
			const h = pool && pool.handlers[i];
			if (h) h(ev.data);
		};
	}
	return {workers: pool.readyCount};
}

/**
 * @param {URL} url
 * @param {boolean} simd
 * @returns {Promise<Worker | null>}
 */
function spawnWorker(url, simd){
	return new Promise((resolve) => {
		let w;
		try {
			w = new Worker(url, {type: 'module'});
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
		w.onerror = () => {
			try { w.terminate(); } catch { /* ignore */ }
			finish(null);
		};
		w.onmessage = (ev) => {
			const msg = ev.data;
			if (msg && msg.type === 'ready'){
				finish(w);
				return;
			}
			if (msg && msg.type === 'error' && msg.id === undefined){
				try { w.terminate(); } catch { /* ignore */ }
				finish(null);
			}
		};
		w.postMessage({type: 'init', simd});
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
