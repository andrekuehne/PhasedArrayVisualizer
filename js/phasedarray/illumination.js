/** @import { PhasedArray } from "./phasedarray.js" */
/**
 * @typedef {IlluminationTypicalPlaneWave} IlluminationHint
 */

export class IlluminationTypicalPlaneWave{
	static title = "Typical (Plane)";
	static args = [];
	static controls = {};
	constructor(){};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		pa.vIllumPhaseFactor.fill(0.0)
		pa.vIllumMag.fill(1.0)
	};
}
export class IlluminationArbitraryPlaneWave{
	static title = "Arbitrary Plane";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3'];
	static controls = {
		'illum-par-1': {'title': "Theta (deg)", 'type': "float", 'default': 0, 'step': 0.1},
		'illum-par-2': {'title': "Phi (deg)", 'type': "float", 'default': 0, 'step': 0.1},
		'illum-par-3': {'title': "Dielectric Constant", 'type': "float", 'default': 1, 'step': 0.1, 'min': 1.0}
	};
	constructor(theta, phi, dk){
		this.theta = theta * Math.PI / 180;
		this.phi = phi * Math.PI / 180;
		this.dk = dk;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const xf = Math.sin(this.theta)*Math.cos(this.phi);
		const yf = Math.sin(this.theta)*Math.sin(this.phi);
		const ne = Math.sqrt(Math.max(1, this.dk));

		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const cx = pa.geometry.x_center;
		const cy = pa.geometry.y_center;
		pa.vIllumMag.fill(1.0)
		for (let i = 0; i < pa.geometry.length; i++){
			let s = (x[i] + cx) * xf + (y[i] + cy) * yf;
			pa.vIllumPhaseFactor[i] = s * ne;
		}
	};
}


export class IlluminationSphericalWave{
	static title = "Spherical";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3', 'illum-par-4'];
	static controls = {
		'illum-par-1': {'title': "Origin X (λ)", 'type': "float", 'default': 0, 'step': 0.1},
		'illum-par-2': {'title': "Origin Y (λ)", 'type': "float", 'default': 0, 'step': 0.1},
		'illum-par-3': {'title': "Origin Z (λ)", 'type': "float", 'default': 10, 'step': 0.1},
		'illum-par-4': {'title': "Dielectric Constant", 'type': "float", 'default': 1, 'step': 0.1, 'min': 1.0}
	};
	constructor(x, y, z, dk){
		this.x = x;
		this.y = y;
		this.z = z;
		this.dk = dk;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const ox = this.x;
		const oy = this.y;
		const oz = this.z;
		const ne = Math.sqrt(Math.max(1, this.dk));

		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const cx = pa.geometry.x_center;
		const cy = pa.geometry.y_center;
		const rs = Float32Array.from({length: pa.geometry.length}, (_, i) => Math.sqrt((ox - x[i] + cx)**2 + (oy - y[i] + cy)**2 + oz**2));
		const min_r = Math.min(...rs);
		for (let i = 0; i < pa.geometry.length; i++){
			pa.vIllumMag[i] = min_r / rs[i];
			pa.vIllumPhaseFactor[i] = rs[i] * ne;
		}
	};
}

export class IlluminationGuidedWaveEntering{
	static title = "Guided (Entering)";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3', 'illum-par-4'];
	static controls = {
		'illum-par-1': {'title': "Relative Origin X", 'type': "float", 'default': -1, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-2': {'title': "Relative Origin Y", 'type': "float", 'default': 0, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-3': {'title': "Angle (deg)", 'type': "float", 'default': 0, 'min': 0, 'max': 360.0, 'step': 1.0,},
		'illum-par-4': {'title': "Dielectric Constant", 'type': "float", 'default': 1, 'step': 0.1, 'min': 1.0}
	};
	constructor(x, y, ang, dk){
		this.x = x;
		this.y = y;
		this.ang = ang;
		this.dk = dk;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const ang = this.ang * Math.PI / 180;
		const ne = Math.sqrt(Math.max(1, this.dk));
		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const min_x = Math.min(...x);
		const max_x = Math.max(...x);
		const min_y = Math.min(...y);
		const max_y = Math.max(...y);
		const cx = (max_x + min_x) / 2 + (max_x - min_x) * this.x / 2;
		const cy = (max_y + min_y) / 2 + (max_y - min_y) * this.y / 2;

		const fcos = Math.cos(ang);
		const fsin = Math.sin(ang);
		pa.vIllumMag.fill(1.0)
		for (let i = 0; i < pa.geometry.length; i++){
			let nx = Math.abs((x[i] - cx) * fcos - (y[i] - cy) * fsin);
			pa.vIllumPhaseFactor[i] = nx * ne;
		}
	};
}
export class IlluminationGuidedWaveLeaving{
	static title = "Guided (Leaving)";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3', 'illum-par-4'];
	static controls = {
		'illum-par-1': {'title': "Relative Origin X", 'type': "float", 'default': -1, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-2': {'title': "Relative Origin Y", 'type': "float", 'default': 0, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-3': {'title': "Angle (deg)", 'type': "float", 'default': 0, 'min': 0, 'max': 360.0, 'step': 1.0,},
		'illum-par-4': {'title': "Dielectric Constant", 'type': "float", 'default': 1, 'step': 0.1, 'min': 1.0}
	};
	constructor(x, y, ang, dk){
		this.x = x;
		this.y = y;
		this.ang = ang;
		this.dk = dk;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const ang = this.ang * Math.PI / 180;
		const ne = Math.sqrt(Math.max(1, this.dk));
		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const min_x = Math.min(...x);
		const max_x = Math.max(...x);
		const min_y = Math.min(...y);
		const max_y = Math.max(...y);
		const cx = (max_x + min_x) / 2 + (max_x - min_x) * this.x / 2;
		const cy = (max_y + min_y) / 2 + (max_y - min_y) * this.y / 2;

		const fcos = Math.cos(ang);
		const fsin = Math.sin(ang);
		pa.vIllumMag.fill(1.0)
		for (let i = 0; i < pa.geometry.length; i++){
			let nx = Math.abs((x[i] - cx) * fcos - (y[i] - cy) * fsin);
			pa.vIllumPhaseFactor[i] = -nx * ne;
		}
	};
}

export class IlluminationOutwardCylindricalWave{
	static title = "Cylindrical (Outward)";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3'];
	static controls = {
		'illum-par-1': {'title': "Relative X", 'type': "float", 'default': 0, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-2': {'title': "Relative Y", 'type': "float", 'default': 0, 'min': -1, 'max': 1.0, 'step': 0.1,},
		'illum-par-3': {'title': "Dielectric Constant", 'type': "float", 'default': 1, 'step': 0.1, 'min': 1.0}
	};
	constructor(x, y, dk){
		this.x = x;
		this.y = y;
		this.dk = dk;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const ne = Math.sqrt(Math.max(1, this.dk));

		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const min_x = Math.min(...x);
		const max_x = Math.max(...x);
		const cx = (max_x + min_x) / 2 + (max_x - min_x) * this.x / 2;
		const min_y = Math.min(...y);
		const max_y = Math.max(...y);
		const cy = (max_y + min_y) / 2 + (max_y - min_y) * this.y / 2;
		pa.vIllumMag.fill(1.0)
		for (let i = 0; i < pa.geometry.length; i++){
			pa.vIllumPhaseFactor[i] = Math.sqrt((cx - x[i])**2 + (cy - y[i])**2) * ne;
		}
	};
}
export class IlluminationInwardCylindricalWave extends IlluminationOutwardCylindricalWave{
	static title = "Cylindrical (Inward)";
	static args = ['illum-par-1', 'illum-par-2', 'illum-par-3'];
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const ne = Math.sqrt(Math.max(1, this.dk));

		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const min_x = Math.min(...x);
		const max_x = Math.max(...x);
		const cx = (max_x + min_x) / 2 + (max_x - min_x) * this.x / 2;
		const min_y = Math.min(...y);
		const max_y = Math.max(...y);
		const cy = (max_y + min_y) / 2 + (max_y - min_y) * this.y / 2;
		const rs = Float32Array.from(x, (ix, i) => Math.sqrt((ix - cx)**2 + (y[i] - cy)**2))
		const max_r = Math.max(...rs);
		pa.vIllumMag.fill(1.0)
		for (let i = 0; i < pa.geometry.length; i++){
			pa.vIllumPhaseFactor[i] = (max_r - rs[i]) * ne;
		}
	};
}
export const Illuminations = [
	IlluminationTypicalPlaneWave,
	IlluminationArbitraryPlaneWave,
	IlluminationSphericalWave,
	IlluminationGuidedWaveLeaving,
	IlluminationGuidedWaveEntering,
	IlluminationOutwardCylindricalWave,
	IlluminationInwardCylindricalWave,
]
