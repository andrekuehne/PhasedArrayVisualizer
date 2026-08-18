/**
 * Module worker: far-field row tiles and radiated-power sample panels.
 * Browser Dedicated Worker or Node worker_threads.
 */

/** @type {import('./simd/farfield_kernel.js').FarfieldKernel | null} */
let kernel = null;
/** @type {import('./simd/farfield_kernel.js').RadiatedPowerKernel | null} */
let prad = null;

/** @type {(msg: object, transfer?: Transferable[]) => void} */
let post = (msg, transfer) => {
	if (transfer && transfer.length) self.postMessage(msg, transfer);
	else self.postMessage(msg);
};

/**
 * @param {object} msg
 */
async function onData(msg){
	if (!msg || !msg.type) return;
	if (msg.type === 'init'){
		try {
			const mod = msg.simd
				? await import('./simd/farfield_kernel.js')
				: await import('./scalar/farfield_kernel.js');
			if (msg.wasmPath){
				const {readFile} = await import('node:fs/promises');
				const bytes = await readFile(msg.wasmPath);
				await mod.default({module_or_path: bytes});
			}
			else {
				await mod.default();
			}
			kernel = new mod.FarfieldKernel();
			prad = new mod.RadiatedPowerKernel();
			post({type: 'ready'});
		}
		catch (e){
			post({type: 'error', message: String(e && e.message ? e.message : e)});
		}
		return;
	}
	if (msg.type === 'run'){
		runFarfield(msg);
		return;
	}
	if (msg.type === 'run_prad'){
		runPrad(msg);
	}
}

/**
 * @param {object} msg
 */
function runFarfield(msg){
	if (kernel === null){
		post({type: 'error', id: msg.id, message: 'Farfield worker kernel is not initialized.'});
		return;
	}
	try {
		const n1 = msg.ax1.length;
		const n2 = msg.ax2.length;
		const tileRows = Math.max(1, msg.tileRows || 32);
		const tiles = Math.max(1, Math.ceil(n2 / tileRows));
		kernel.prepare(n1, n2);
		kernel.set_inputs(msg.x, msg.y, msg.mag, msg.pha, msg.ax1, msg.ax2);
		for (let t = 0; t < tiles; t++){
			const row0 = t * tileRows;
			const rowCount = Math.min(tileRows, n2 - row0);
			post({type: 'progress', id: msg.id, done: t, total: tiles});
			kernel.accumulate_tile(msg.domain, msg.frequencyScale, row0, rowCount);
		}
		post({type: 'progress', id: msg.id, done: tiles, total: tiles});
		const maxValue = kernel.finalize(msg.nElements);
		const total = kernel.take_total();
		post({type: 'result', id: msg.id, total, maxValue}, [total.buffer]);
	}
	catch (e){
		post({type: 'error', id: msg.id, message: String(e && e.message ? e.message : e)});
	}
}

/**
 * @param {object} msg
 */
function runPrad(msg){
	if (prad === null){
		post({type: 'error', id: msg.id, message: 'Radiated-power worker kernel is not initialized.'});
		return;
	}
	try {
		post({type: 'progress', id: msg.id, done: 0, total: 2});
		prad.set_quadrature(msg.nMu, msg.nPhi);
		prad.fill_isolated_range(
			msg.x,
			msg.y,
			msg.frequencyScale,
			msg.elementKind,
			msg.elementN,
			msg.sample0,
			msg.sampleCount
		);
		post({type: 'progress', id: msg.id, done: 1, total: 2});
		prad.form_gram();
		const re = prad.take_re();
		const im = prad.take_im();
		post({type: 'progress', id: msg.id, done: 2, total: 2});
		post({type: 'result', id: msg.id, re, im, n: prad.n_elements()}, [re.buffer, im.buffer]);
	}
	catch (e){
		post({type: 'error', id: msg.id, message: String(e && e.message ? e.message : e)});
	}
}

if (typeof WorkerGlobalScope !== 'undefined'){
	self.onmessage = (ev) => {
		onData(ev.data);
	};
}
else {
	const {parentPort} = await import('node:worker_threads');
	post = (msg, transfer) => {
		parentPort.postMessage(msg, transfer || []);
	};
	parentPort.on('message', onData);
}
