import {ones, normalize, zeros} from "../util.js";
import { FeedUniform } from "./feed.js"
/** @import { GeometryHint } from "./geometry.js" */
/** @import { FeedHint } from "./feed.js" */

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

		this.vectorPhaseRawFactor = new Float32Array(this.size);
		this.vectorPhaseManual = zeros(this.size);
		this.vectorPhaseIsManual = Array.from({length: this.size}, () => false);
		this.vectorPhaseSubtract = new Float32Array(this.size);
		this.vectorPhaseQuantizeFactor = new Float32Array(this.size);
		this.vectorPhasePureFactor = new Float32Array(this.size);
		this.vectorPhaseFarfieldCalc = new Float32Array(this.size);

		this.vectorPhaseIllumFactor = new Float32Array(this.size);
		this.vectorMagIllum = new Float32Array(this.size);
		this.vectorIllumAtten = new Float32Array(this.size);

		this.vectorMagRaw = new Float32Array(this.size);
		this.vectorMagManual = zeros(this.size);
		this.vectorMagIsManual = Array.from({length: this.size}, () => false);
		this.elementDisabled = Array.from({length: this.size}, () => false);
		this.vectorMag = new Float32Array(this.size);
		this.vectorMagFarfieldCalc = new Float32Array(this.size);
		this.vectorMagPure = new Float32Array(this.size);
		this.vectorAtten = new Float32Array(this.size);
		this.requestUpdate = true;
		this.feed = null;
	}
	set_theta_phi(theta, phi){
		this.theta = Number(theta);
		this.phi = Number(phi);
		this.requestUpdate = true;
	}
	/**
	* Create a Phased Array object.
	*
	* @param {FeedHint} feed
	* */
	set_feed_type(feed){
		this.feed = feed;
	}
	compute_illumination(){
		if (this.feed === null) this.feed = new FeedUniform();
		this.feed.compute_illumination(this);
		this.vectorIllumAtten = Float32Array.from(this.vectorMagIllum, v => 20 * Math.log10(Math.abs(v)))
	}
	compute_phase(){
		const xf = Math.sin(this.theta*Math.PI/180)*Math.cos(this.phi*Math.PI/180);
		const yf = Math.sin(this.theta*Math.PI/180)*Math.sin(this.phi*Math.PI/180);

		const x = this.geometry.x;
		const y = this.geometry.y;
		const cx = this.geometry.x_center;
		const cy = this.geometry.y_center;
		const p = -2 * Math.PI;
		for (let i = 0; i < this.geometry.length; i++){
			this.vectorPhaseRawFactor[i] = (((x[i] + cx) * xf + (y[i] + cy) * yf) - this.vectorPhaseIllumFactor[i]);
		}
	}
	calculate_r_taper(){
		const x = normalize(this.geometry.x);
		const y = normalize(this.geometry.y);
		const r = Float32Array.from(x, (ix, i) => Math.sqrt(ix**2 + y[i]**2));
		const maxR = Math.max(...r);
		this.vectorMagX = this.taperX.calculate_weights(Float32Array.from(r, (v) => v/maxR*0.5));
		this.vectorMagY = ones(this.vectorMagX.length);
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
		this.vectorMagRaw = vector;
		this.requestUpdate = true;
	}
	calculate_final_vector(){
		const p = 2 * Math.PI;
		for (let i = 0; i < this.geometry.length; i++){
			let p, m;
			if (this.vectorPhaseIsManual[i]) p = this.vectorPhaseManual[i] / p;
			else p = this.vectorPhaseRawFactor[i];
			if (this.vectorMagIsManual[i]) {
				if (this.elementDisabled[i]) m = 0;
				else m = this.vectorMagManual[i];
			}
			else m = this.vectorMagRaw[i];
			if (m < 0) {
				m = Math.abs(m)
				p += 0.5;
			}
			this.vectorPhasePureFactor[i] = p;
			this.vectorMagPure[i] = m;
		}
		this.requestUpdate = false;
	}
	quantize_phase(bits, dither){
		let lsb = 0;
		const f = 2 * Math.PI;
		if (dither === undefined) dither = false;
		if (bits > 0) lsb = 1/2**bits;

		for (let i = 0; i < this.geometry.length; i++){
			let p = this.vectorPhasePureFactor[i] % 1.0;
			if (bits <= 0) this.vectorPhaseQuantizeFactor[i] = p;
			else{
				if (dither) p += lsb / 2 * Math.round(Math.random())
				while (p < 0) p += 1;
				this.vectorPhaseQuantizeFactor[i] = lsb * Math.round(p / lsb);
			}
			this.vectorPhaseFarfieldCalc[i] = -f * ((this.vectorPhaseQuantizeFactor[i] + this.vectorPhaseIllumFactor[i]) % 1.0);
		}
	}
	quantize_attenuation(bits, lsb){
		const maxQ = lsb*(2**bits - 1);
		const maxV = Math.max(...this.vectorMagPure);
		for (let i = 0; i < this.geometry.length; i++){
			let m = this.vectorMagPure[i]/maxV;
			let a = -20*Math.log10(Math.abs(m));
			if (bits <= 0 || lsb <= 0){
				this.vectorMag[i] = m;
				this.vectorAtten[i] = -a;
			}
			else{
				a = lsb*Math.round(a/lsb);
				if (a > maxQ){
					this.vectorAtten[i] = -Infinity;
					this.vectorMag[i] = 0;
				}
				else{
					this.vectorAtten[i] = -a;
					this.vectorMag[i] = 10**(-a/20.0);
				}
			}
			this.vectorMagFarfieldCalc[i] = this.vectorMag[i] * this.vectorMagIllum[i];
		}
	}
}
