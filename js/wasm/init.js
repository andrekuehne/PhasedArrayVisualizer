/**
 * Load the far-field WASM kernel (SIMD build when the engine supports it).
 */

import {startFarfieldPool} from "./farfield-pool.js";

const SIMD_TEST = new Uint8Array([
	0, 97, 115, 109, 1, 0, 0, 0, 1, 5, 1, 96, 0, 1, 123, 3, 2, 1, 0, 10, 10, 1, 8, 0, 65, 0, 253, 15, 253, 98, 11
]);

/** @type {import('./simd/farfield_kernel.js').FarfieldKernel | null} */
let kernel = null;
/** @type {import('./simd/farfield_kernel.js').RadiatedPowerKernel | null} */
let pradKernel = null;
/** @type {typeof import('./simd/farfield_kernel.js').extract_pattern_metrics | null} */
let extractFn = null;
/** @type {typeof import('./simd/farfield_kernel.js').apply_element_pattern | null} */
let applyElementFn = null;
/** @type {typeof import('./simd/farfield_kernel.js').apply_green_pec_pattern | null} */
let applyGreenPecFn = null;
/** @type {typeof import('./simd/farfield_kernel.js').element_exponent_from_peak_dbi | null} */
let exponentFn = null;
/** @type {typeof import('./simd/farfield_kernel.js').z_self_pec_dipole | null} */
let zSelfPecFn = null;

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

/**
 * @returns {import('./simd/farfield_kernel.js').RadiatedPowerKernel}
 */
export function getRadiatedPowerKernel(){
	if (pradKernel === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return pradKernel;
}

/**
 * Extract HPBW, sidelobe, and beam-pointing metrics from a computed intensity map.
 * @param {number} domain
 * @param {Float32Array} ax1
 * @param {Float32Array} ax2
 * @param {Float32Array} total
 * @param {number} reqThetaRad requested spherical theta (rad)
 * @param {number} reqPhiRad requested spherical phi (rad)
 * @returns {import('./simd/farfield_kernel.js').PatternMetrics}
 */
export function extractPatternMetrics(domain, ax1, ax2, total, reqThetaRad, reqPhiRad){
	if (extractFn === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return extractFn(domain, ax1, ax2, total, reqThetaRad, reqPhiRad);
}

/**
 * Multiply far-field intensity by an element power pattern. Mutates `total`.
 * @param {number} domain
 * @param {Float32Array} ax1
 * @param {Float32Array} ax2
 * @param {Float32Array} total
 * @param {number} kind
 * @param {number} n
 * @returns {number} peak intensity after apply
 */
export function applyElementPattern(domain, ax1, ax2, total, kind, n){
	if (applyElementFn === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return applyElementFn(domain, ax1, ax2, total, kind, n);
}

/**
 * Multiply far-field intensity by PEC-dipole |F^iso|^2. Mutates `total`.
 * @param {number} domain
 * @param {Float32Array} ax1
 * @param {Float32Array} ax2
 * @param {Float32Array} total
 * @param {number} h height in wavelengths at f0
 * @param {number} ell dipole length in wavelengths at f0
 * @param {number} freqScale
 * @returns {number} peak intensity after apply
 */
export function applyGreenPecPattern(domain, ax1, ax2, total, h, ell, freqScale){
	if (applyGreenPecFn === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return applyGreenPecFn(domain, ax1, ax2, total, h, ell, freqScale);
}

/**
 * Power-conserving cos^n exponent from peak element gain in dBi.
 * @param {number} gainDbi
 * @returns {number}
 */
export function elementExponentFromPeakDbi(gainDbi){
	if (exponentFn === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return exponentFn(gainDbi);
}

/**
 * Isolated PEC-dipole self impedance Z11 in ohms.
 * @param {number} h
 * @param {number} ell
 * @param {number} a
 * @param {number} freqScale
 * @returns {Float64Array} `[re, im]`
 */
export function zSelfPecDipole(h, ell, a, freqScale){
	if (zSelfPecFn === null){
		throw new Error('Farfield WASM kernel is not initialized.');
	}
	return zSelfPecFn(h, ell, a, freqScale);
}

export async function initFarfieldWasm(){
	const useSimd = wasmSupportsSimd();
	const mod = useSimd
		? await import('./simd/farfield_kernel.js')
		: await import('./scalar/farfield_kernel.js');
	await mod.default();
	kernel = new mod.FarfieldKernel();
	pradKernel = new mod.RadiatedPowerKernel();
	extractFn = mod.extract_pattern_metrics;
	applyElementFn = mod.apply_element_pattern;
	applyGreenPecFn = mod.apply_green_pec_pattern;
	exponentFn = mod.element_exponent_from_peak_dbi;
	zSelfPecFn = mod.z_self_pec_dipole;
	try {
		await startFarfieldPool({simd: useSimd});
	}
	catch {
		// Main-thread kernel remains available.
	}
	return {simd: useSimd};
}
