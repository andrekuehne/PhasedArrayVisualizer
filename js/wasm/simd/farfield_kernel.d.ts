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

/**
 * Pattern-feature metrics extracted from a computed intensity map.
 */
export class PatternMetrics {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    hpbw_ax1_clipped: boolean;
    hpbw_ax1_deg: number;
    hpbw_ax1: number;
    hpbw_ax2_clipped: boolean;
    hpbw_ax2_deg: number;
    hpbw_ax2: number;
    hpbw_large_angle_deg: number;
    hpbw_large_clipped: boolean;
    hpbw_large_deg: number;
    hpbw_large: number;
    hpbw_small_angle_deg: number;
    hpbw_small_clipped: boolean;
    hpbw_small_deg: number;
    hpbw_small: number;
    largest_sll_ax1: number;
    largest_sll_ax2: number;
    largest_sll_db: number;
    nearest_sll_ax1: number;
    nearest_sll_ax2: number;
    nearest_sll_db: number;
    peak_ax1: number;
    peak_ax2: number;
    peak_i1: number;
    peak_i2: number;
    peak_phi_deg: number;
    peak_theta_deg: number;
    requested_phi_deg: number;
    requested_theta_deg: number;
    squint_ax1_deg: number;
    squint_ax2_deg: number;
    squint_deg: number;
}

export class RadiatedPowerKernel {
    free(): void;
    [Symbol.dispose](): void;
    compute(x: Float32Array, y: Float32Array, frequency_scale: number, element_kind: number, element_n: number): void;
    compute_j0(x: Float32Array, y: Float32Array, frequency_scale: number, element_kind: number, element_n: number): void;
    fill_isolated(x: Float32Array, y: Float32Array, frequency_scale: number, element_kind: number, element_n: number): void;
    fill_isolated_range(x: Float32Array, y: Float32Array, frequency_scale: number, element_kind: number, element_n: number, sample0: number, sample_count: number): void;
    form_gram(): void;
    form_matched_s(z_ref: number): void;
    match_iterations(): number;
    match_residual(): number;
    n_elements(): number;
    n_samples(): number;
    constructor();
    set_quadrature(n_mu: number, n_phi: number): void;
    take_im(): Float32Array;
    take_re(): Float32Array;
    take_s_im(): Float64Array;
    take_s_re(): Float64Array;
    take_t_im(): Float64Array;
    take_t_re(): Float64Array;
    take_z0(): Float64Array;
}

/**
 * Multiply AF intensity `total` (row-major `i2 * n1 + i1`) by the element
 * power pattern. Returns the new peak. Cos^n uses `[max(w,0)]^n` and is 0
 * for invisible/back directions even when `n == 0` (`0^0` would otherwise be 1).
 */
export function apply_element_pattern(domain: number, ax1: Float32Array, ax2: Float32Array, total: Float32Array, kind: number, n: number): number;

/**
 * Power-conserving cos^n exponent from peak element gain in dBi:
 * `n = 10^(element_gain/10)/2 - 1`, clamped at 0.
 */
export function element_exponent_from_peak_dbi(gain_dbi: number): number;

export function extract_pattern_metrics(domain: number, ax1: Float32Array, ax2: Float32Array, total: Float32Array, req_theta_rad: number, req_phi_rad: number): PatternMetrics;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_farfieldkernel_free: (a: number, b: number) => void;
    readonly __wbg_get_patternmetrics_hpbw_ax1: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_ax1_clipped: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_ax1_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_ax2: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_ax2_clipped: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_ax2_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_large: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_large_angle_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_large_clipped: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_large_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_small: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_small_angle_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_small_clipped: (a: number) => number;
    readonly __wbg_get_patternmetrics_hpbw_small_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_largest_sll_ax1: (a: number) => number;
    readonly __wbg_get_patternmetrics_largest_sll_ax2: (a: number) => number;
    readonly __wbg_get_patternmetrics_largest_sll_db: (a: number) => number;
    readonly __wbg_get_patternmetrics_nearest_sll_ax1: (a: number) => number;
    readonly __wbg_get_patternmetrics_nearest_sll_ax2: (a: number) => number;
    readonly __wbg_get_patternmetrics_nearest_sll_db: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_ax1: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_ax2: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_i1: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_i2: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_phi_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_peak_theta_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_requested_phi_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_requested_theta_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_squint_ax1_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_squint_ax2_deg: (a: number) => number;
    readonly __wbg_get_patternmetrics_squint_deg: (a: number) => number;
    readonly __wbg_patternmetrics_free: (a: number, b: number) => void;
    readonly __wbg_radiatedpowerkernel_free: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax1: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax1_clipped: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax1_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax2: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax2_clipped: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_ax2_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_large: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_large_angle_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_large_clipped: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_large_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_small: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_small_angle_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_small_clipped: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_hpbw_small_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_largest_sll_ax1: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_largest_sll_ax2: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_largest_sll_db: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_nearest_sll_ax1: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_nearest_sll_ax2: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_nearest_sll_db: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_ax1: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_ax2: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_i1: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_i2: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_phi_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_peak_theta_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_requested_phi_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_requested_theta_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_squint_ax1_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_squint_ax2_deg: (a: number, b: number) => void;
    readonly __wbg_set_patternmetrics_squint_deg: (a: number, b: number) => void;
    readonly apply_element_pattern: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: any, i: number, j: number) => number;
    readonly element_exponent_from_peak_dbi: (a: number) => number;
    readonly extract_pattern_metrics: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => number;
    readonly farfieldkernel_accumulate_tile: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly farfieldkernel_finalize: (a: number, b: number) => number;
    readonly farfieldkernel_new: () => number;
    readonly farfieldkernel_prepare: (a: number, b: number, c: number) => void;
    readonly farfieldkernel_set_inputs: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => void;
    readonly farfieldkernel_take_total: (a: number) => [number, number];
    readonly radiatedpowerkernel_compute: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => void;
    readonly radiatedpowerkernel_compute_j0: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => void;
    readonly radiatedpowerkernel_fill_isolated: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => void;
    readonly radiatedpowerkernel_fill_isolated_range: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => void;
    readonly radiatedpowerkernel_form_gram: (a: number) => void;
    readonly radiatedpowerkernel_form_matched_s: (a: number, b: number) => void;
    readonly radiatedpowerkernel_match_iterations: (a: number) => number;
    readonly radiatedpowerkernel_match_residual: (a: number) => number;
    readonly radiatedpowerkernel_n_elements: (a: number) => number;
    readonly radiatedpowerkernel_n_samples: (a: number) => number;
    readonly radiatedpowerkernel_new: () => number;
    readonly radiatedpowerkernel_set_quadrature: (a: number, b: number, c: number) => void;
    readonly radiatedpowerkernel_take_im: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_re: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_s_im: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_s_re: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_t_im: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_t_re: (a: number) => [number, number];
    readonly radiatedpowerkernel_take_z0: (a: number) => [number, number];
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
