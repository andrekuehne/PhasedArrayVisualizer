import {ones, normalize, zeros} from "../util.js";
import { IlluminationTypicalPlaneWave } from "./illumination.js"
/** @import { GeometryHint } from "./geometry.js" */
/** @import { IlluminationHint } from "./illumination.js" */

export class PhasedArray{
	/**
	* Create a Phased Array object.
	*
	* @param {GeometryHint} geometry
	* */
	constructor(geometry){
		this.geometry = geometry;
		this.set_theta_phi(0, 0);
		this.size = geometry.length;

		// ideal vectors that the element should be assigned to before quantization.
		this.vIdealPhaseFactor = new Float32Array(this.size);

		// requested vectors. this is the ideal vector with manual override applied.
		this.vRequestedPhaseFactor = new Float32Array(this.size);
		this.vRequestedMag = new Float32Array(this.size);

		// tasked vectors. These are assigned after steering and attenuation.
		this.vSteerPhaseFactor = new Float32Array(this.size);
		this.vTaperMag = new Float32Array(this.size);

		// Illumination vector that is illumination the array.
		this.vIllumMag = new Float32Array(this.size);
		this.vIllumPhaseFactor = new Float32Array(this.size);

		// qunatized vector.
		this.vQuantizePhaseFactor = new Float32Array(this.size);
		this.vQuantizeMag = new Float32Array(this.size);

		// vector ready for farfield calculation.
		this.vFarfieldMag = new Float32Array(this.size);

		// vector overrides.
		this.vectorPhaseManual = zeros(this.size);
		this.vectorPhaseIsManual = Array.from({length: this.size}, () => false);
		this.vectorMagManual = zeros(this.size);
		this.vectorMagIsManual = Array.from({length: this.size}, () => false);
		this.elementDisabled = Array.from({length: this.size}, () => false);

		this.requestUpdate = true;
		this.illum = null;
		this.powerWeightSum = 0;
	}
	set_theta_phi(theta, phi){
		this.theta = Number(theta);
		this.phi = Number(phi);
		this.requestUpdate = true;
	}
	create_taper_vectors(){
		const vmag = new Float32Array(this.size);
		const vpha = new Float32Array(this.size);
		const vector = this.vTaperMag;
		for (let i = 0; i < this.geometry.length; i++){
			let m = vector[i];
			vmag[i] = Math.abs(m);
			if (m < 0) vpha[i] = 0.5;
			else vpha[i] = 0.0;
		}
		return [vpha, vmag];
	}
	create_final_vectors(){
		return [this.vQuantizePhaseFactor, this.vFarfieldMag];
	}
	create_steer_vectors(){
		return [this.vSteerPhaseFactor, ones(this.size)];
	}
	create_ideal_vectors(){
		const vmag = new Float32Array(this.size);
		const vpha = new Float32Array(this.size);
		const omag = this.vTaperMag;
		const opha = this.vIdealPhaseFactor;

		for (let i = 0; i < this.geometry.length; i++){
			let m = omag[i];
			let p = opha[i];
			vmag[i] = Math.abs(m);
			if (m < 0) p += 0.5;
			vpha[i] = p % 1.0;
		}
		return [vpha, vmag];
	}
	/**
	* Set illumination.
	*
	* @param {IlluminationHint} illum
	* */
	set_illumination_type(illum){
		this.illum = illum;
	}
	compute_illumination(){
		if (this.illum === null) this.illum = new IlluminationTypicalPlaneWave();
		this.illum.compute_illumination(this);
	}
	compute_phase(){
		const xf = Math.sin(this.theta*Math.PI/180)*Math.cos(this.phi*Math.PI/180);
		const yf = Math.sin(this.theta*Math.PI/180)*Math.sin(this.phi*Math.PI/180);

		const x = this.geometry.x;
		const y = this.geometry.y;
		const cx = this.geometry.x_center;
		const cy = this.geometry.y_center;
		for (let i = 0; i < this.geometry.length; i++){
			let s = (x[i] + cx) * xf + (y[i] + cy) * yf;
			this.vSteerPhaseFactor[i] = s;
			this.vIdealPhaseFactor[i] = s - this.vIllumPhaseFactor[i];
		}
	}
	set_manual_phase(index, override, phaseRad){
		let ov = this.vectorPhaseIsManual[index];
		if (ov === false && override === false) return;
		let cp = this.vectorPhaseManual[index];
		if (ov == override && cp == phaseRad) return;
		this.vectorPhaseIsManual[index] = Boolean(override);
		this.vectorPhaseManual[index] = phaseRad;
		this.requestUpdate = true;
	}
	clear_all_manual_phase(){
		for (let i = 0; i < this.geometry.length; i++) this.vectorPhaseIsManual[i] = false;
		this.requestUpdate = true;
	}
	set_manual_magnitude(index, override, mag, disable){
		let ov = this.vectorMagIsManual[index];
		if (ov === false && override === false) return;
		let cp = this.vectorMagManual[index];
		let cd = this.elementDisabled[index];
		if (ov == override && cp == mag && cd == disable) return;
		this.vectorMagIsManual[index] = Boolean(override);
		this.vectorMagManual[index] = mag;
		this.elementDisabled[index] = disable;
		this.requestUpdate = true;
	}
	clear_all_manual_magnitude(){
		for (let i = 0; i < this.geometry.length; i++) {
			this.vectorMagIsManual[i] = false;
			this.elementDisabled[i] = false;
		}
		this.requestUpdate = true;
	}
	set_magnitude_weight(vector){
		this.vTaperMag = vector;
		this.requestUpdate = true;
	}
	calculate_requested_vector(){
		const pf = 2 * Math.PI;
		for (let i = 0; i < this.geometry.length; i++){
			let p, m;
			if (this.vectorPhaseIsManual[i]) p = this.vectorPhaseManual[i] / pf;
			else p = this.vIdealPhaseFactor[i];
			if (this.vectorMagIsManual[i]) {
				if (this.elementDisabled[i]) m = 0;
				else m = this.vectorMagManual[i];
			}
			else m = this.vTaperMag[i];
			if (m < 0) {
				m = Math.abs(m)
				p += 0.5;
			}
			this.vRequestedPhaseFactor[i] = p;
			this.vRequestedMag[i] = m;
		}
		this.requestUpdate = false;
	}
	quantize_phase(bits, dither){
		let lsb = 0;
		if (dither === undefined) dither = false;
		if (bits > 0) lsb = 1/2**bits;

		for (let i = 0; i < this.geometry.length; i++){
			let p = this.vRequestedPhaseFactor[i] % 1.0;
			if (bits <= 0) this.vQuantizePhaseFactor[i] = p;
			else{
				if (dither) p += lsb / 2 * Math.round(Math.random())
				while (p < 0) p += 1;
				this.vQuantizePhaseFactor[i] = lsb * Math.round(p / lsb);
			}
		}
	}
	quantize_attenuation(bits, lsb){
		const maxQ = lsb*(2**bits - 1);
		const maxV = Math.max(...this.vRequestedMag);
		let powerWeightSum = 0;
		for (let i = 0; i < this.geometry.length; i++){
			let m = this.vRequestedMag[i]/maxV;
			let a = -20*Math.log10(Math.abs(m));
			if (bits <= 0 || lsb <= 0) this.vQuantizeMag[i] = m;
			else{
				a = lsb*Math.round(a/lsb);
				if (a > maxQ) this.vQuantizeMag[i] = 0;
				else this.vQuantizeMag[i] = 10**(-a/20.0);
			}
			this.vFarfieldMag[i] = this.vQuantizeMag[i] * this.vIllumMag[i];
			const q = this.vQuantizeMag[i];
			powerWeightSum += q * q;
		}
		this.powerWeightSum = powerWeightSum;
	}
	/**
	* Power delivered to element `i` in watts, given full-scale per-antenna power.
	*
	* @param {Number} pAnt Full-scale per-antenna power (W)
	* @param {Number} i Element index
	*
	* @return {Number}
	*/
	elementPowerWatts(pAnt, i){
		const m = this.vQuantizeMag[i];
		return Number(pAnt) * m * m;
	}
	/**
	* Total power sent into the array in watts, given full-scale per-antenna power.
	*
	* @param {Number} pAnt Full-scale per-antenna power (W)
	*
	* @return {Number}
	*/
	totalPowerWatts(pAnt){
		return Number(pAnt) * this.powerWeightSum;
	}
	create_farfield_vectors(freq_scale){
		const pha = new Float32Array(this.size);
		const mag = new Float32Array(this.size);
		const f = 2 * Math.PI;
		for (let i = 0; i < this.size; i++){
			pha[i] = -f * ((this.vQuantizePhaseFactor[i] + this.vIllumPhaseFactor[i] * freq_scale) % 1.0);
		}
		return [pha, this.vFarfieldMag];
	}
}
