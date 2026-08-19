/**
 * Matched-port helpers: T a, S a, conjugate embedded phases, J0 quadrature order.
 */
import {PATTERN_COS_N, PATTERN_ISOTROPIC} from "./element.js";

export const Z_REF = 50;

/** Auto-switch Coupling to Isolated when a new geometry exceeds this many elements. */
export const MATCHED_AUTO_ISOLATE_N = 512;
/** Skip Green PEC LU (T = I) when N exceeds this. 32×32 stays matched. */
export const GREEN_PEC_AUTO_ISOLATE_N = 1024;

/**
 * Hemispheric power-conserving D(μ) used by the radiated-power Gram.
 * @param {number} mu cosθ
 * @param {number} kind PATTERN_ISOTROPIC | PATTERN_COS_N
 * @param {number} n
 * @returns {number}
 */
export function elementDirectivity(mu, kind, n){
	if (mu < 0) return 0;
	if (kind === PATTERN_COS_N){
		const nn = Number.isFinite(n) ? Math.max(0, n) : 0;
		if (nn === 0) return 2;
		return 2 * (nn + 1) * Math.pow(mu, nn);
	}
	return 2;
}

/**
 * Gauss-μ order from electrical diameter: ceil(π D_λ) + 12.
 * @param {{r?: ArrayLike<number>}} geometry
 * @param {number} frequencyScale
 * @returns {number}
 */
export function nMuFromGeometry(geometry, frequencyScale){
	let rMax = 0;
	const r = geometry && geometry.r;
	if (r){
		for (let i = 0; i < r.length; i++){
			const v = r[i];
			if (v > rMax) rMax = v;
		}
	}
	const f = Number(frequencyScale);
	const dLambda = 2 * rMax * (Number.isFinite(f) ? Math.abs(f) : 1);
	return Math.max(8, Math.ceil(Math.PI * dLambda) + 12);
}

/**
 * Identity T (isolated conjugate / missing cache).
 * @param {number} n
 * @returns {{re: Float64Array, im: Float64Array}}
 */
export function identityT(n){
	const nn = n * n;
	const re = new Float64Array(nn);
	const im = new Float64Array(nn);
	for (let i = 0; i < n; i++) re[i * n + i] = 1;
	return {re, im};
}

/**
 * Complex GEMV w = M a. M is row-major N×N.
 * @param {ArrayLike<number>} mRe
 * @param {ArrayLike<number>} mIm
 * @param {ArrayLike<number>} aRe
 * @param {ArrayLike<number>} aIm
 * @returns {{re: Float64Array, im: Float64Array}}
 */
export function gemv(mRe, mIm, aRe, aIm){
	const n = aRe.length;
	const re = new Float64Array(n);
	const im = new Float64Array(n);
	const hasIm = mIm && mIm.length === n * n;
	for (let i = 0; i < n; i++){
		let wr = 0;
		let wi = 0;
		const row = i * n;
		for (let j = 0; j < n; j++){
			const mr = mRe[row + j];
			const mi = hasIm ? mIm[row + j] : 0;
			const ar = aRe[j];
			const ai = aIm[j];
			wr += mr * ar - mi * ai;
			wi += mr * ai + mi * ar;
		}
		re[i] = wr;
		im[i] = wi;
	}
	return {re, im};
}

/**
 * Reflected-power fraction Γ = ||S a||² / ||a||².
 * @param {ArrayLike<number>} sRe
 * @param {ArrayLike<number>} sIm
 * @param {ArrayLike<number>} aRe
 * @param {ArrayLike<number>} aIm
 * @returns {number}
 */
export function reflectionRatio(sRe, sIm, aRe, aIm){
	let a2 = 0;
	for (let i = 0; i < aRe.length; i++) a2 += aRe[i] * aRe[i] + aIm[i] * aIm[i];
	if (!(a2 > 0)) return 0;
	const sa = gemv(sRe, sIm, aRe, aIm);
	let s2 = 0;
	for (let i = 0; i < sa.re.length; i++) s2 += sa.re[i] * sa.re[i] + sa.im[i] * sa.im[i];
	return s2 / a2;
}

/**
 * Conjugate-match phase in cycles: arg(F_emb) / 2π so kernel pha = -arg(F_emb).
 * F_emb = T^T F_iso with F_iso_p = sqrt(D) exp(j k r_p · rhat).
 *
 * @param {ArrayLike<number>} x
 * @param {ArrayLike<number>} y
 * @param {number} thetaDeg
 * @param {number} phiDeg
 * @param {number} frequencyScale
 * @param {ArrayLike<number>|null} tRe
 * @param {ArrayLike<number>|null} tIm
 * @param {number} elementKind
 * @param {number} elementN
 * @returns {Float32Array}
 */
export function conjugatePhaseCycles(
	x, y, thetaDeg, phiDeg, frequencyScale, tRe, tIm, elementKind, elementN
){
	const n = x.length;
	const th = Number(thetaDeg) * Math.PI / 180;
	const ph = Number(phiDeg) * Math.PI / 180;
	const st = Math.sin(th);
	const u = st * Math.cos(ph);
	const v = st * Math.sin(ph);
	const mu = Math.cos(th);
	const k = 2 * Math.PI * Number(frequencyScale);
	const amp = Math.sqrt(Math.max(0, elementDirectivity(mu, elementKind, elementN)));
	const fRe = new Float64Array(n);
	const fIm = new Float64Array(n);
	for (let p = 0; p < n; p++){
		const phase = k * (x[p] * u + y[p] * v);
		fRe[p] = amp * Math.cos(phase);
		fIm[p] = amp * Math.sin(phase);
	}
	const T = (tRe && tRe.length === n * n) ? {re: tRe, im: tIm} : identityT(n);
	const cycles = new Float32Array(n);
	const twoPi = 2 * Math.PI;
	const hasIm = T.im && T.im.length === n * n;
	for (let kk = 0; kk < n; kk++){
		let er = 0;
		let ei = 0;
		for (let p = 0; p < n; p++){
			const tr = T.re[p * n + kk];
			const ti = hasIm ? T.im[p * n + kk] : 0;
			er += tr * fRe[p] - ti * fIm[p];
			ei += tr * fIm[p] + ti * fRe[p];
		}
		cycles[kk] = Math.atan2(ei, er) / twoPi;
	}
	return cycles;
}

/**
 * Shift `source` by one global cycle offset so it sits as close as possible
 * to `reference` as unit-magnitude complex vectors:
 * α = exp(j arg Σ exp(j 2π (ref − src))), then unwrap onto the reference branch.
 *
 * @param {ArrayLike<number>} source
 * @param {ArrayLike<number>} reference
 * @returns {Float32Array}
 */
export function alignPhaseCycles(source, reference){
	const n = source.length;
	const twoPi = 2 * Math.PI;
	let re = 0;
	let im = 0;
	for (let i = 0; i < n; i++){
		const d = twoPi * (reference[i] - source[i]);
		re += Math.cos(d);
		im += Math.sin(d);
	}
	const offset = (re === 0 && im === 0) ? 0 : Math.atan2(im, re) / twoPi;
	const out = new Float32Array(n);
	for (let i = 0; i < n; i++){
		let d = source[i] + offset - reference[i];
		d -= Math.round(d);
		out[i] = reference[i] + d;
	}
	return out;
}
