/**
 * Load the far-field WASM kernel (SIMD build when the engine supports it).
 */

import {startFarfieldPool} from "./farfield-pool.js";

const SIMD_TEST = new Uint8Array([
	0, 97, 115, 109, 1, 0, 0, 0, 1, 5, 1, 96, 0, 1, 123, 3, 2, 1, 0, 10, 10, 1, 8, 0, 65, 0, 253, 15, 253, 98, 11
]);

/** @type {import('./simd/farfield_kernel.js').FarfieldKernel | null} */
let kernel = null;

export function wasmSupportsSimd(){
	try {
		return WebAssembly.validate(SIMD_TEST);
	}
	catch {
		return false;
	}
}

/**
 * @returns {import('./simd/farfield_kernel.js').FarfieldKernel}
 */
export function getFarfieldKernel(){
	if (kernel === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return kernel;
}

export async function initFarfieldWasm(){
	const useSimd = wasmSupportsSimd();
	const mod = useSimd
		? await import('./simd/farfield_kernel.js')
		: await import('./scalar/farfield_kernel.js');
	await mod.default();
	kernel = new mod.FarfieldKernel();
	try {
		await startFarfieldPool({simd: useSimd});
	}
	catch {
		// Main-thread kernel remains available.
	}
	return {simd: useSimd};
}
