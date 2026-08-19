/* @ts-self-types="./farfield_kernel.d.ts" */

export class FarfieldKernel {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        FarfieldKernelFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_farfieldkernel_free(ptr, 0);
    }
    /**
     * @param {number} domain
     * @param {number} frequency_scale
     * @param {number} row0
     * @param {number} row_count
     */
    accumulate_tile(domain, frequency_scale, row0, row_count) {
        wasm.farfieldkernel_accumulate_tile(this.__wbg_ptr, domain, frequency_scale, row0, row_count);
    }
    /**
     * @param {number} n_elements
     * @returns {number}
     */
    finalize(n_elements) {
        const ret = wasm.farfieldkernel_finalize(this.__wbg_ptr, n_elements);
        return ret;
    }
    constructor() {
        const ret = wasm.farfieldkernel_new();
        this.__wbg_ptr = ret;
        FarfieldKernelFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * @param {number} n1
     * @param {number} n2
     */
    prepare(n1, n2) {
        wasm.farfieldkernel_prepare(this.__wbg_ptr, n1, n2);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {Float32Array} mag
     * @param {Float32Array} pha
     * @param {Float32Array} ax1
     * @param {Float32Array} ax2
     */
    set_inputs(x, y, mag, pha, ax1, ax2) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passArrayF32ToWasm0(mag, wasm.__wbindgen_malloc);
        const len2 = WASM_VECTOR_LEN;
        const ptr3 = passArrayF32ToWasm0(pha, wasm.__wbindgen_malloc);
        const len3 = WASM_VECTOR_LEN;
        const ptr4 = passArrayF32ToWasm0(ax1, wasm.__wbindgen_malloc);
        const len4 = WASM_VECTOR_LEN;
        const ptr5 = passArrayF32ToWasm0(ax2, wasm.__wbindgen_malloc);
        const len5 = WASM_VECTOR_LEN;
        wasm.farfieldkernel_set_inputs(this.__wbg_ptr, ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, ptr5, len5);
    }
    /**
     * @returns {Float32Array}
     */
    take_total() {
        const ret = wasm.farfieldkernel_take_total(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
}
if (Symbol.dispose) FarfieldKernel.prototype[Symbol.dispose] = FarfieldKernel.prototype.free;

/**
 * Pattern-feature metrics extracted from a computed intensity map.
 */
export class PatternMetrics {
    static __wrap(ptr) {
        const obj = Object.create(PatternMetrics.prototype);
        obj.__wbg_ptr = ptr;
        PatternMetricsFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        PatternMetricsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_patternmetrics_free(ptr, 0);
    }
    /**
     * @returns {boolean}
     */
    get hpbw_ax1_clipped() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax1_clipped(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    get hpbw_ax1_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax1_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_ax1() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax1(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {boolean}
     */
    get hpbw_ax2_clipped() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax2_clipped(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    get hpbw_ax2_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax2_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_ax2() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_ax2(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_large_angle_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_large_angle_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {boolean}
     */
    get hpbw_large_clipped() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_large_clipped(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    get hpbw_large_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_large_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_large() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_large(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_small_angle_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_small_angle_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {boolean}
     */
    get hpbw_small_clipped() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_small_clipped(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    get hpbw_small_deg() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_small_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get hpbw_small() {
        const ret = wasm.__wbg_get_patternmetrics_hpbw_small(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get largest_sll_ax1() {
        const ret = wasm.__wbg_get_patternmetrics_largest_sll_ax1(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get largest_sll_ax2() {
        const ret = wasm.__wbg_get_patternmetrics_largest_sll_ax2(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get largest_sll_db() {
        const ret = wasm.__wbg_get_patternmetrics_largest_sll_db(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get nearest_sll_ax1() {
        const ret = wasm.__wbg_get_patternmetrics_nearest_sll_ax1(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get nearest_sll_ax2() {
        const ret = wasm.__wbg_get_patternmetrics_nearest_sll_ax2(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get nearest_sll_db() {
        const ret = wasm.__wbg_get_patternmetrics_nearest_sll_db(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get peak_ax1() {
        const ret = wasm.__wbg_get_patternmetrics_peak_ax1(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get peak_ax2() {
        const ret = wasm.__wbg_get_patternmetrics_peak_ax2(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get peak_i1() {
        const ret = wasm.__wbg_get_patternmetrics_peak_i1(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    get peak_i2() {
        const ret = wasm.__wbg_get_patternmetrics_peak_i2(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    get peak_phi_deg() {
        const ret = wasm.__wbg_get_patternmetrics_peak_phi_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get peak_theta_deg() {
        const ret = wasm.__wbg_get_patternmetrics_peak_theta_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get requested_phi_deg() {
        const ret = wasm.__wbg_get_patternmetrics_requested_phi_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get requested_theta_deg() {
        const ret = wasm.__wbg_get_patternmetrics_requested_theta_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get squint_ax1_deg() {
        const ret = wasm.__wbg_get_patternmetrics_squint_ax1_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get squint_ax2_deg() {
        const ret = wasm.__wbg_get_patternmetrics_squint_ax2_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    get squint_deg() {
        const ret = wasm.__wbg_get_patternmetrics_squint_deg(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {boolean} arg0
     */
    set hpbw_ax1_clipped(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax1_clipped(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_ax1_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax1_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_ax1(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax1(this.__wbg_ptr, arg0);
    }
    /**
     * @param {boolean} arg0
     */
    set hpbw_ax2_clipped(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax2_clipped(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_ax2_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax2_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_ax2(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_ax2(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_large_angle_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_large_angle_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {boolean} arg0
     */
    set hpbw_large_clipped(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_large_clipped(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_large_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_large_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_large(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_large(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_small_angle_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_small_angle_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {boolean} arg0
     */
    set hpbw_small_clipped(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_small_clipped(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_small_deg(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_small_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set hpbw_small(arg0) {
        wasm.__wbg_set_patternmetrics_hpbw_small(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set largest_sll_ax1(arg0) {
        wasm.__wbg_set_patternmetrics_largest_sll_ax1(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set largest_sll_ax2(arg0) {
        wasm.__wbg_set_patternmetrics_largest_sll_ax2(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set largest_sll_db(arg0) {
        wasm.__wbg_set_patternmetrics_largest_sll_db(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set nearest_sll_ax1(arg0) {
        wasm.__wbg_set_patternmetrics_nearest_sll_ax1(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set nearest_sll_ax2(arg0) {
        wasm.__wbg_set_patternmetrics_nearest_sll_ax2(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set nearest_sll_db(arg0) {
        wasm.__wbg_set_patternmetrics_nearest_sll_db(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_ax1(arg0) {
        wasm.__wbg_set_patternmetrics_peak_ax1(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_ax2(arg0) {
        wasm.__wbg_set_patternmetrics_peak_ax2(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_i1(arg0) {
        wasm.__wbg_set_patternmetrics_peak_i1(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_i2(arg0) {
        wasm.__wbg_set_patternmetrics_peak_i2(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_phi_deg(arg0) {
        wasm.__wbg_set_patternmetrics_peak_phi_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set peak_theta_deg(arg0) {
        wasm.__wbg_set_patternmetrics_peak_theta_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set requested_phi_deg(arg0) {
        wasm.__wbg_set_patternmetrics_requested_phi_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set requested_theta_deg(arg0) {
        wasm.__wbg_set_patternmetrics_requested_theta_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set squint_ax1_deg(arg0) {
        wasm.__wbg_set_patternmetrics_squint_ax1_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set squint_ax2_deg(arg0) {
        wasm.__wbg_set_patternmetrics_squint_ax2_deg(this.__wbg_ptr, arg0);
    }
    /**
     * @param {number} arg0
     */
    set squint_deg(arg0) {
        wasm.__wbg_set_patternmetrics_squint_deg(this.__wbg_ptr, arg0);
    }
}
if (Symbol.dispose) PatternMetrics.prototype[Symbol.dispose] = PatternMetrics.prototype.free;

export class RadiatedPowerKernel {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RadiatedPowerKernelFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_radiatedpowerkernel_free(ptr, 0);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} element_kind
     * @param {number} element_n
     */
    compute(x, y, frequency_scale, element_kind, element_n) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_compute(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, element_kind, element_n);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} element_kind
     * @param {number} element_n
     */
    compute_j0(x, y, frequency_scale, element_kind, element_n) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_compute_j0(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, element_kind, element_n);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} h
     * @param {number} ell
     * @param {number} a
     */
    fill_green_pec_dipole_z(x, y, frequency_scale, h, ell, a) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_fill_green_pec_dipole_z(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, h, ell, a);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} element_kind
     * @param {number} element_n
     */
    fill_isolated(x, y, frequency_scale, element_kind, element_n) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_fill_isolated(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, element_kind, element_n);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} element_kind
     * @param {number} element_n
     * @param {number} sample0
     * @param {number} sample_count
     */
    fill_isolated_range(x, y, frequency_scale, element_kind, element_n, sample0, sample_count) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_fill_isolated_range(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, element_kind, element_n, sample0, sample_count);
    }
    /**
     * @param {number} z_ref
     * @param {number} z_common_re
     */
    form_from_z(z_ref, z_common_re) {
        wasm.radiatedpowerkernel_form_from_z(this.__wbg_ptr, z_ref, z_common_re);
    }
    form_gram() {
        wasm.radiatedpowerkernel_form_gram(this.__wbg_ptr);
    }
    /**
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} frequency_scale
     * @param {number} h
     * @param {number} ell
     * @param {number} a
     * @param {number} z_ref
     * @param {number} z_common_re
     */
    form_green_pec_dipole(x, y, frequency_scale, h, ell, a, z_ref, z_common_re) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_form_green_pec_dipole(this.__wbg_ptr, ptr0, len0, ptr1, len1, frequency_scale, h, ell, a, z_ref, z_common_re);
    }
    /**
     * @param {number} z_ref
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} x_nn
     * @param {number} alpha
     * @param {number} beta
     * @param {number} aniso
     * @param {number} z_common_re
     * @param {number} x_self
     */
    form_matched_s(z_ref, x, y, x_nn, alpha, beta, aniso, z_common_re, x_self) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_form_matched_s(this.__wbg_ptr, z_ref, ptr0, len0, ptr1, len1, x_nn, alpha, beta, aniso, z_common_re, x_self);
    }
    /**
     * @param {number} z_ref
     * @param {Float32Array} x
     * @param {Float32Array} y
     * @param {number} x_nn
     * @param {number} att
     * @param {number} eps_x
     * @param {number} eps_y
     * @param {number} freq
     * @param {number} z_common_re
     * @param {number} x_self
     */
    form_matched_s_propagation(z_ref, x, y, x_nn, att, eps_x, eps_y, freq, z_common_re, x_self) {
        const ptr0 = passArrayF32ToWasm0(x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF32ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        wasm.radiatedpowerkernel_form_matched_s_propagation(this.__wbg_ptr, z_ref, ptr0, len0, ptr1, len1, x_nn, att, eps_x, eps_y, freq, z_common_re, x_self);
    }
    /**
     * @returns {number}
     */
    match_iterations() {
        const ret = wasm.radiatedpowerkernel_match_iterations(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    match_residual() {
        const ret = wasm.radiatedpowerkernel_match_residual(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    n_elements() {
        const ret = wasm.radiatedpowerkernel_n_elements(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    n_samples() {
        const ret = wasm.radiatedpowerkernel_n_samples(this.__wbg_ptr);
        return ret >>> 0;
    }
    constructor() {
        const ret = wasm.radiatedpowerkernel_new();
        this.__wbg_ptr = ret;
        RadiatedPowerKernelFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * @param {number} n_mu
     * @param {number} n_phi
     */
    set_quadrature(n_mu, n_phi) {
        wasm.radiatedpowerkernel_set_quadrature(this.__wbg_ptr, n_mu, n_phi);
    }
    /**
     * @returns {Float32Array}
     */
    take_im() {
        const ret = wasm.radiatedpowerkernel_take_im(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @returns {Float32Array}
     */
    take_re() {
        const ret = wasm.radiatedpowerkernel_take_re(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_s_im() {
        const ret = wasm.radiatedpowerkernel_take_s_im(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_s_re() {
        const ret = wasm.radiatedpowerkernel_take_s_re(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_t_im() {
        const ret = wasm.radiatedpowerkernel_take_t_im(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_t_re() {
        const ret = wasm.radiatedpowerkernel_take_t_re(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_z0() {
        const ret = wasm.radiatedpowerkernel_take_z0(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_z0_im() {
        const ret = wasm.radiatedpowerkernel_take_z0_im(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_z_im() {
        const ret = wasm.radiatedpowerkernel_take_z_im(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    take_z_re() {
        const ret = wasm.radiatedpowerkernel_take_z_re(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
}
if (Symbol.dispose) RadiatedPowerKernel.prototype[Symbol.dispose] = RadiatedPowerKernel.prototype.free;

/**
 * Multiply AF intensity `total` (row-major `i2 * n1 + i1`) by the element
 * power pattern. Returns the new peak. Cos^n uses `[max(w,0)]^n` and is 0
 * for invisible/back directions even when `n == 0` (`0^0` would otherwise be 1).
 * @param {number} domain
 * @param {Float32Array} ax1
 * @param {Float32Array} ax2
 * @param {Float32Array} total
 * @param {number} kind
 * @param {number} n
 * @returns {number}
 */
export function apply_element_pattern(domain, ax1, ax2, total, kind, n) {
    const ptr0 = passArrayF32ToWasm0(ax1, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF32ToWasm0(ax2, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    var ptr2 = passArrayF32ToWasm0(total, wasm.__wbindgen_malloc);
    var len2 = WASM_VECTOR_LEN;
    const ret = wasm.apply_element_pattern(domain, ptr0, len0, ptr1, len1, ptr2, len2, total, kind, n);
    return ret;
}

/**
 * Power-conserving cos^n exponent from peak element gain in dBi:
 * `n = 10^(element_gain/10)/2 - 1`, clamped at 0.
 * @param {number} gain_dbi
 * @returns {number}
 */
export function element_exponent_from_peak_dbi(gain_dbi) {
    const ret = wasm.element_exponent_from_peak_dbi(gain_dbi);
    return ret;
}

/**
 * @param {number} domain
 * @param {Float32Array} ax1
 * @param {Float32Array} ax2
 * @param {Float32Array} total
 * @param {number} req_theta_rad
 * @param {number} req_phi_rad
 * @returns {PatternMetrics}
 */
export function extract_pattern_metrics(domain, ax1, ax2, total, req_theta_rad, req_phi_rad) {
    const ptr0 = passArrayF32ToWasm0(ax1, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF32ToWasm0(ax2, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF32ToWasm0(total, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.extract_pattern_metrics(domain, ptr0, len0, ptr1, len1, ptr2, len2, req_theta_rad, req_phi_rad);
    return PatternMetrics.__wrap(ret);
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_copy_to_typed_array_c7f28e53671b41e8: function(arg0, arg1, arg2) {
            new Uint8Array(arg2.buffer, arg2.byteOffset, arg2.byteLength).set(getArrayU8FromWasm0(arg0, arg1));
        },
        __wbg___wbindgen_throw_bb96b2010945f0bc: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./farfield_kernel_bg.js": import0,
    };
}

const FarfieldKernelFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_farfieldkernel_free(ptr, 1));
const PatternMetricsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_patternmetrics_free(ptr, 1));
const RadiatedPowerKernelFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_radiatedpowerkernel_free(ptr, 1));

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

let cachedFloat32ArrayMemory0 = null;
function getFloat32ArrayMemory0() {
    if (cachedFloat32ArrayMemory0 === null || cachedFloat32ArrayMemory0.byteLength === 0) {
        cachedFloat32ArrayMemory0 = new Float32Array(wasm.memory.buffer);
    }
    return cachedFloat32ArrayMemory0;
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function passArrayF32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getFloat32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedFloat32ArrayMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (!module.ok) {
            throw new Error(`failed to fetch Wasm: ${module.status} ${module.statusText} fetching '${module.url}'`);
        }

        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('farfield_kernel_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
