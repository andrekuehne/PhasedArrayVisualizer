export const PATTERN_ISOTROPIC = 0;
export const PATTERN_COS_N = 1;
export const PATTERN_GREEN_PEC = 2;
export const PATTERN_GREEN_SLAB = 3;
export const GREEN_PEC_DEFAULT_H = 0.25;
export const GREEN_PEC_DEFAULT_ELL = 0.1;
export const GREEN_PEC_DEFAULT_A = 0.001;
export const GREEN_PEC_DEFAULT_Z0 = 50;
export const GREEN_PEC_DEFAULT_XSELF = 0;
export const GREEN_SLAB_DEFAULT_EPS_R = 10;
export const GREEN_SLAB_DEFAULT_H_SUB = 0.05;
export const GREEN_SLAB_DEFAULT_TAN_DELTA = 0;

/** Hemispherically isotropic peak gain: 10·log10(2) ≈ 3.01 dBi. */
export const MIN_ELEMENT_GAIN_DBI = 10 * Math.log10(2);

/**
 * Power-conserving cos^n exponent: `n = 10^(element_gain/10)/2 - 1`.
 * @param {number} gainDbi
 * @returns {number}
 */
export function exponentFromPeakDbi(gainDbi){
	const g = Math.pow(10, Number(gainDbi) / 10);
	if (!Number.isFinite(g)) return 0;
	return Math.max(0, g / 2 - 1);
}

export class ElementIsotropic{
	static title = "Isotropic";
	static args = [];
	static controls = {};
	constructor(){
		this.kind = PATTERN_ISOTROPIC;
		this.n = 0;
	}
}

export class ElementCosN{
	static title = "cos^n";
	static args = ['element-gain'];
	static controls = {
		'element-gain': {
			title: "Gain (dBi)",
			type: "float",
			default: 5,
			step: 0.1,
			min: MIN_ELEMENT_GAIN_DBI,
			max: 20,
		},
	};
	constructor(gainDbi){
		this.kind = PATTERN_COS_N;
		this.gainDbi = Number(gainDbi);
		this.n = exponentFromPeakDbi(this.gainDbi);
	}
}

const GREEN_SHARED_CONTROLS = {
	'element-h': {
		title: "h (λ)",
		type: "float",
		default: GREEN_PEC_DEFAULT_H,
		step: 0.01,
		min: 0.001,
		max: 2,
		help: "Height of the dipole center above z=0 (PEC ground, or the dielectric surface for a slab), in wavelengths at f0. Default 0.25.",
	},
	'element-ell': {
		title: "ℓ (λ)",
		type: "float",
		default: GREEN_PEC_DEFAULT_ELL,
		step: 0.01,
		min: 0.01,
		max: 0.5,
		help: "Dipole length in wavelengths at f0. Short-dipole model; keep well below element spacing. Must be greater than 2a.",
	},
	'element-a': {
		title: "a (λ)",
		type: "float",
		default: GREEN_PEC_DEFAULT_A,
		step: 0.0005,
		min: 0.0001,
		max: 0.05,
		help: "Equivalent cylindrical radius in wavelengths at f0. Sets free-space self reactance through ln(ℓ/2a).",
	},
	'element-z0': {
		title: "Z0 (Ω)",
		type: "float",
		default: GREEN_PEC_DEFAULT_Z0,
		step: 1,
		min: 1,
		max: 500,
		help: "Common real source impedance on every port. Kurokawa S and T use this real zc (not a complex matching network).",
	},
	'element-xself': {
		title: "Self X (Ω)",
		type: "float",
		default: GREEN_PEC_DEFAULT_XSELF,
		step: 1,
		help: "Series reactance added to every diagonal of Z (like a series inductor). Does not replace the physical Xii. Default 0.",
	},
};

/** Horizontal short dipole over infinite PEC. Z and F^iso from the same Green function. */
export class ElementGreenPec{
	static title = "PEC dipole";
	static help = "Horizontal short dipole over infinite PEC. Mutual Z and the isolated pattern come from the same Green function.";
	static args = ['element-h', 'element-ell', 'element-a', 'element-z0', 'element-xself'];
	static controls = GREEN_SHARED_CONTROLS;
	constructor(h, ell, a, z0, xself){
		this.kind = PATTERN_GREEN_PEC;
		let hh = Number(h);
		let ee = Number(ell);
		let aa = Number(a);
		let zz = Number(z0);
		let xx = Number(xself);
		if (!Number.isFinite(hh) || hh <= 0) hh = GREEN_PEC_DEFAULT_H;
		if (!Number.isFinite(aa) || aa <= 0) aa = GREEN_PEC_DEFAULT_A;
		if (!Number.isFinite(ee) || ee <= 0) ee = GREEN_PEC_DEFAULT_ELL;
		if (ee <= 2 * aa) ee = Math.max(GREEN_PEC_DEFAULT_ELL, 2 * aa * 1.001);
		if (!Number.isFinite(zz) || zz <= 0) zz = GREEN_PEC_DEFAULT_Z0;
		if (!Number.isFinite(xx)) xx = GREEN_PEC_DEFAULT_XSELF;
		this.h = hh;
		this.ell = ee;
		this.a = aa;
		this.z0 = zz;
		this.xself = xx;
		this.n = 0;
	}
}

/** Horizontal short dipole over a grounded dielectric slab. Z and F^iso from the same Green function. */
export class ElementGreenSlab{
	static title = "Slab dipole";
	static help = "Horizontal short dipole over a grounded dielectric slab. Mutual Z and the isolated pattern come from the same Green function. Surface-wave power is not in the array-factor pattern.";
	static args = ['element-h', 'element-ell', 'element-a', 'element-eps-r', 'element-h-sub', 'element-tan-delta', 'element-z0', 'element-xself'];
	static controls = {
		...GREEN_SHARED_CONTROLS,
		'element-eps-r': {
			title: "ε_r",
			type: "float",
			default: GREEN_SLAB_DEFAULT_EPS_R,
			step: 0.1,
			min: 1,
			max: 100,
			help: "Relative permittivity of the grounded slab.",
		},
		'element-h-sub': {
			title: "h_sub (λ)",
			type: "float",
			default: GREEN_SLAB_DEFAULT_H_SUB,
			step: 0.005,
			min: 0.0001,
			max: 1,
			help: "Substrate thickness in wavelengths at f0. PEC is at z = −h_sub.",
		},
		'element-tan-delta': {
			title: "tan δ",
			type: "float",
			default: GREEN_SLAB_DEFAULT_TAN_DELTA,
			step: 0.001,
			min: 0,
			max: 1,
			help: "Dielectric loss tangent. 0 is lossless.",
		},
	};
	constructor(h, ell, a, epsR, hSub, tanDelta, z0, xself){
		this.kind = PATTERN_GREEN_SLAB;
		let hh = Number(h);
		let ee = Number(ell);
		let aa = Number(a);
		let er = Number(epsR);
		let hs = Number(hSub);
		let td = Number(tanDelta);
		let zz = Number(z0);
		let xx = Number(xself);
		if (!Number.isFinite(hh) || hh <= 0) hh = GREEN_PEC_DEFAULT_H;
		if (!Number.isFinite(aa) || aa <= 0) aa = GREEN_PEC_DEFAULT_A;
		if (!Number.isFinite(ee) || ee <= 0) ee = GREEN_PEC_DEFAULT_ELL;
		if (ee <= 2 * aa) ee = Math.max(GREEN_PEC_DEFAULT_ELL, 2 * aa * 1.001);
		if (!Number.isFinite(er) || er < 1) er = GREEN_SLAB_DEFAULT_EPS_R;
		if (!Number.isFinite(hs) || hs <= 0) hs = GREEN_SLAB_DEFAULT_H_SUB;
		if (!Number.isFinite(td) || td < 0) td = GREEN_SLAB_DEFAULT_TAN_DELTA;
		if (!Number.isFinite(zz) || zz <= 0) zz = GREEN_PEC_DEFAULT_Z0;
		if (!Number.isFinite(xx)) xx = GREEN_PEC_DEFAULT_XSELF;
		this.h = hh;
		this.ell = ee;
		this.a = aa;
		this.epsR = er;
		this.hSub = hs;
		this.tanDelta = td;
		this.z0 = zz;
		this.xself = xx;
		this.n = 0;
	}
}

export const ElementTypes = [
	ElementIsotropic,
	ElementCosN,
	ElementGreenPec,
	ElementGreenSlab,
];
