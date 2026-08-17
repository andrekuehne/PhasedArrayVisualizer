/* tslint:disable */
/* eslint-disable */

export class FarfieldKernel {
    free(): void;
    [Symbol.dispose](): void;
    accumulate_tile(domain: number, frequency_scale: number, row0: number, row_count: number): void;
    finalize(n_elements: number): number;
    constructor();
    prepare(n1: number, n2: number): void;
    set_inputs(x: Float32Array, y: Float32Array, mag: Float32Array, pha: Float32Array, ax1: Float32Array, ax2: Float32Array): void;
    take_total(): Float32Array;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_farfieldkernel_free: (a: number, b: number) => void;
    readonly farfieldkernel_accumulate_tile: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly farfieldkernel_finalize: (a: number, b: number) => number;
    readonly farfieldkernel_new: () => number;
    readonly farfieldkernel_prepare: (a: number, b: number, c: number) => void;
    readonly farfieldkernel_set_inputs: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => void;
    readonly farfieldkernel_take_total: (a: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
