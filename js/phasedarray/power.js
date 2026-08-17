export const POWER_UNITS = ['mW', 'W', 'dBm', 'dBW'];
export const POWER_SCOPES = ['element', 'array'];

const DEFAULT_DBM = 20;

/**
 * Convert a power value in `unit` to watts.
 *
 * @param {Number} value
 * @param {String} unit One of W, mW, dBm, dBW
 *
 * @return {Number}
 */
export function wattsFrom(value, unit){
	if (value === '' || value === null || value === undefined) return NaN;
	const v = Number(value);
	if (!Number.isFinite(v)) return NaN;
	if (unit === 'W'){
		if (v < 0) return NaN;
		return v;
	}
	if (unit === 'mW'){
		if (v < 0) return NaN;
		return v / 1000;
	}
	if (unit === 'dBW') return 10 ** (v / 10);
	if (unit === 'dBm') return 10 ** ((v - 30) / 10);
	throw Error(`Unknown power unit ${unit}`);
}

/**
 * Convert watts to a power value in `unit`.
 *
 * @param {Number} watts
 * @param {String} unit One of W, mW, dBm, dBW
 *
 * @return {Number}
 */
export function wattsTo(watts, unit){
	const w = Number(watts);
	if (!Number.isFinite(w) || w < 0) return NaN;
	if (unit === 'W') return w;
	if (unit === 'mW') return w * 1000;
	if (w === 0) return -Infinity;
	if (unit === 'dBW') return 10 * Math.log10(w);
	if (unit === 'dBm') return 10 * Math.log10(w) + 30;
	throw Error(`Unknown power unit ${unit}`);
}

/**
 * Default full-scale per-antenna power in watts (20 dBm).
 *
 * @return {Number}
 */
export function defaultWatts(){
	return wattsFrom(DEFAULT_DBM, 'dBm');
}

function formatLinear(v){
	if (!Number.isFinite(v)) return '';
	if (v === 0) return '0';
	const av = Math.abs(v);
	if (av >= 1e4 || av < 1e-4) return v.toExponential(4);
	return String(Number(v.toPrecision(6)));
}

function formatDb(v, decimals){
	if (v === -Infinity) return '-Inf';
	if (!Number.isFinite(v)) return '';
	return v.toFixed(decimals);
}

/**
 * Format watts as an input-field value in `unit` (no unit suffix).
 *
 * @param {Number} watts
 * @param {String} unit
 *
 * @return {String}
 */
export function formatPowerValue(watts, unit){
	const v = wattsTo(watts, unit);
	if (unit === 'dBm' || unit === 'dBW'){
		if (v === -Infinity) return '-Inf';
		if (!Number.isFinite(v)) return '';
		return String(Number(v.toFixed(4)));
	}
	return formatLinear(v);
}

/**
 * Format watts for display, including the unit suffix.
 *
 * @param {Number} watts
 * @param {String} unit
 * @param {Number} [decimals=2] Decimal places for dBm/dBW
 *
 * @return {String}
 */
export function formatPower(watts, unit, decimals){
	if (decimals === undefined) decimals = 2;
	const v = wattsTo(watts, unit);
	if (unit === 'dBm' || unit === 'dBW') return `${formatDb(v, decimals)} ${unit}`;
	const t = formatLinear(v);
	if (t === '') return '';
	return `${t} ${unit}`;
}

/**
 * True if `unit` is a known power unit.
 *
 * @param {String} unit
 *
 * @return {Boolean}
 */
export function isPowerUnit(unit){
	return POWER_UNITS.includes(unit);
}

/**
 * True if `scope` is per-element or full-array.
 *
 * @param {String} scope
 *
 * @return {Boolean}
 */
export function isPowerScope(scope){
	return POWER_SCOPES.includes(scope);
}
