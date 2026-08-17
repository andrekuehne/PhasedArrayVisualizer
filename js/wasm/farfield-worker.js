/**
 * Module worker: one FarfieldKernel instance, jobs are ax2-row slices.
 */

/** @type {import('./simd/farfield_kernel.js').FarfieldKernel | null} */
let kernel = null;

self.onmessage = async (ev) => {
	const msg = ev.data;
	if (msg.type === 'init'){
		try {
			const mod = msg.simd
				? await import('./simd/farfield_kernel.js')
				: await import('./scalar/farfield_kernel.js');
			await mod.default();
			kernel = new mod.FarfieldKernel();
			self.postMessage({type: 'ready'});
		}
		catch (e){
			self.postMessage({type: 'error', message: String(e && e.message ? e.message : e)});
		}
		return;
	}
	if (msg.type !== 'run') return;
	if (kernel === null){
		self.postMessage({type: 'error', id: msg.id, message: 'Farfield worker kernel is not initialized.'});
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
			self.postMessage({type: 'progress', id: msg.id, done: t, total: tiles});
			kernel.accumulate_tile(msg.domain, msg.frequencyScale, row0, rowCount);
		}
		self.postMessage({type: 'progress', id: msg.id, done: tiles, total: tiles});
		const maxValue = kernel.finalize(msg.nElements);
		const total = kernel.take_total();
		self.postMessage({type: 'result', id: msg.id, total, maxValue}, [total.buffer]);
	}
	catch (e){
		self.postMessage({type: 'error', id: msg.id, message: String(e && e.message ? e.message : e)});
	}
};
