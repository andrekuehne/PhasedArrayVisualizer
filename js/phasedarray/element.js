export const PATTERN_ISOTROPIC = 0;
export const PATTERN_COS_N = 1;

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

export const ElementTypes = [
	ElementIsotropic,
	ElementCosN,
];
