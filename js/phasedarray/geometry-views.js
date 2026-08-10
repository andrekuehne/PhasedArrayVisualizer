
/** @import { PhasedArray } from "../../phasedarray/phasedarray.js" */

export class GeometryViewABC{
	static title = "";
	static args = [];
	static controls = {};
	static allow_manual = false;
	static desc = "";
}
export class GeometryViewConversion{
	static user_title = "";
	static title = "";
	static args = [];
	static controls = {};
	static show_scale = false;
	static popup_type = null;
	default_cmap(){return "viridis";}
	convert_from_view(view, pa, scale){
		let v = view.retrieve_vectors(pa);
		return this.convert(v[0], v[1], scale)
	}
	string_from_view(view, pa, index){
		let v = view.retrieve_vectors(pa);
		return this.value_to_string(v[0][index], v[1][index])
	}
	phase_factor_to_string(phase_factor){return `${((phase_factor % 1.0) * 360).toFixed(2)} deg`}
	mag_to_string(mag){return `${(20 * Math.log10(mag)).toFixed(2)} dB`}

	scale_phase(v){
		const min = -180;
		const max = 180;
		const pd = max - min;
		const res = new Float32Array(v.length);
		for (let i = 0; i < v.length; i++){
			let pha = v[i] * 360;
			while (pha > 180) pha -= 360;
			while (pha < -180) pha += 360;
			res[i] = (pha - min) / pd;
		}
		return res;
	}
	scale_atten(v, scale){
		const ma = Math.max(...v);
		scale = Math.abs(scale);
		return Float32Array.from(v, e => {
			let o = -20 * Math.log10(e / ma) / scale;
			if (o > 1) o = NaN;
			return o;
		});
	}
}
export class Atten extends GeometryViewConversion{
	static user_title = "Magnitude";
	static title = "dB";
	static show_scale = true;
	static popup_type = "atten";
	default_cmap(){return "inferno_r";}
	convert(phaseFactor, magnitude, scale){ return this.scale_atten(magnitude, scale) };
	value_to_string(phaseFactor, magnitude){return this.mag_to_string(magnitude);}
}

export class Phase extends GeometryViewConversion{
	static user_title = "Phase";
	static title = "deg";
	static popup_type = "phase";
	default_cmap(){return "hsv";}
	convert(phaseFactor, magnitude, scale){ return this.scale_phase(phaseFactor) };
	value_to_string(phaseFactor, magnitude){return this.phase_factor_to_string(phaseFactor);}
}

export class Element extends GeometryViewABC{
	static title = "Element";
	static desc = "Phase and amplitude of each element after quantization."
	static allow_manual = true;
	retrieve_vectors(pa){return pa.create_final_vectors()}
}

export class Illumination extends GeometryViewABC{
	static title = "Illumination";
	static desc = "Phase and amplitude of each element's illumination/incoming wave."
	retrieve_vectors(pa){return [pa.vIllumPhaseFactor, pa.vIllumMag]}
}

export class Steer extends GeometryViewABC{
	static title = "Steer";
	static desc = "Phase and amplitude of each element's steering vector."
	retrieve_vectors(pa){return pa.create_steer_vectors()}
}

export class Taper extends GeometryViewABC{
	static title = "Taper";
	static desc = "Phase and amplitude of each element's tapering vector."
	retrieve_vectors(pa){return pa.create_taper_vectors()}
}

export class Ideal extends GeometryViewABC{
	static title = "Ideal";
	static desc = "Phase and amplitude of each element after steering/tapering and illumination compensation but before manual override and quantization."
	retrieve_vectors(pa){return pa.create_ideal_vectors()}
}

export class Requested extends GeometryViewABC{
	static title = "Requested";
	static desc = "Phase and amplitude of each element after calculations and manual overrides but before quantization."
	retrieve_vectors(pa){return [pa.vRequestedPhaseFactor, pa.vRequestedMag]}
}

export class Quantized extends GeometryViewABC{
	static title = "Quantized";
	static desc = "Quantized phase and amplitude of each element after calculations and manual overrides."
	retrieve_vectors(pa){return [pa.vRequestedPhaseFactor, pa.vRequestedMag]}
}

export const GeometryViews = [
	Element,
	Steer,
	Taper,
	Ideal,
	Requested,
	Quantized,
	Illumination,
]

export const GeometryViewConversions = [
	Atten,
	Phase,
]
