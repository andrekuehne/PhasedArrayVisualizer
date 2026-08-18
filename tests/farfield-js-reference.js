/**
 * Independent JS far-field loops used as the equivalence reference for the WASM kernel.
 * Keep this as a line-by-line reference of the intended math. Do not "optimize"
 * it to match the Rust kernel — tests exist to catch those drifts.
 */

/**
 * @param {Float32Array[]} re
 * @param {Float32Array[]} im
 * @param {number} nElements
 */
export function jsCalculateTotal(re, im, nElements){
	const p2 = re.length;
	const p1 = re[0].length;
	const total = new Float32Array(p2 * p1);
	let maxValue = -Infinity;
	for (let i2 = 0; i2 < p2; i2++){
		for (let i1 = 0; i1 < p1; i1++){
			const c = Math.abs(re[i2][i1] ** 2 + im[i2][i1] ** 2) / nElements;
			total[i2 * p1 + i1] = c;
			if (c > maxValue) maxValue = c;
		}
	}
	return {total, maxValue};
}

function allocMesh(n1, n2){
	const re = new Array(n2);
	const im = new Array(n2);
	for (let i = 0; i < n2; i++){
		re[i] = new Float32Array(n1);
		im[i] = new Float32Array(n1);
	}
	return {re, im};
}

/**
 * Original FarfieldSpherical.calculator_loop inner accumulation.
 */
export function jsSpherical(x, y, mag, pha, theta, phi, frequencyScale){
	const {re, im} = allocMesh(theta.length, phi.length);
	const sc = 2 * Math.PI * frequencyScale;
	const sinThetaPi = Float32Array.from(theta, (t) => sc * Math.sin(t));
	for (let i = 0; i < x.length; i++){
		for (let ip = 0; ip < phi.length; ip++){
			const xxv = x[i] * Math.cos(phi[ip]);
			const yyv = y[i] * Math.sin(phi[ip]);
			for (let it = 0; it < theta.length; it++){
				const jk = sinThetaPi[it];
				const v = xxv * jk + yyv * jk + pha[i];
				re[ip][it] += mag[i] * Math.cos(v);
				im[ip][it] += mag[i] * Math.sin(v);
			}
		}
	}
	return jsCalculateTotal(re, im, x.length);
}

/**
 * UV array-factor accumulation. Geometric k is 2π·frequencyScale, matching spherical.
 */
export function jsUV(x, y, mag, pha, u, v, frequencyScale){
	const {re, im} = allocMesh(u.length, v.length);
	const sc = 2 * Math.PI * frequencyScale;
	for (let i = 0; i < x.length; i++){
		for (let iv = 0; iv < v.length; iv++){
			const xxv = x[i];
			const yyv = y[i] * v[iv];
			for (let iu = 0; iu < u.length; iu++){
				const phase = (xxv * u[iu] + yyv) * sc + pha[i];
				re[iv][iu] += mag[i] * Math.cos(phase);
				im[iv][iu] += mag[i] * Math.sin(phase);
			}
		}
	}
	return jsCalculateTotal(re, im, x.length);
}

/**
 * Ludwig-3 array-factor accumulation. Geometric k is 2π·frequencyScale, matching spherical.
 */
export function jsLudwig3(x, y, mag, pha, az, el, frequencyScale){
	const {re, im} = allocMesh(az.length, el.length);
	const sc = 2 * Math.PI * frequencyScale;
	for (let i = 0; i < x.length; i++){
		for (let iv = 0; iv < el.length; iv++){
			const xxv = x[i] * Math.cos(el[iv]);
			const yyv = y[i] * Math.sin(el[iv]);
			for (let iu = 0; iu < az.length; iu++){
				const w = (xxv * Math.sin(az[iu]) + yyv) * sc + pha[i];
				re[iv][iu] += mag[i] * Math.cos(w);
				im[iv][iu] += mag[i] * Math.sin(w);
			}
		}
	}
	return jsCalculateTotal(re, im, x.length);
}
