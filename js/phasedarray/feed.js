
/** @import { PhasedArray } from "./phasedarray.js" */
/**
 * @typedef {FeedUniform} FeedHint
 */

export class FeedUniform{
	static title = "Phased Array";
	static args = [];
	static controls = {};
	constructor(){};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		pa.vectorPhaseIllumFactor.fill(0.0)
		pa.vectorMagIllum.fill(1.0)
	};
}

export class FeedReflectArray{
	static title = "Reflect Array";
	static args = ['feed-par-1', 'feed-par-2', 'feed-par-3'];
	static controls = {
		'feed-par-1': {'title': "Feed X (λ)", 'type': "float", 'default': 0, 'step': 0.1},
		'feed-par-2': {'title': "Feed Y (λ)", 'type': "float", 'default': 0, 'step': 0.1},
		'feed-par-3': {'title': "Feed Z (λ)", 'type': "float", 'default': 10, 'step': 0.1}
	};
	constructor(x, y, z){
		this.x = x;
		this.y = y;
		this.z = z;
	};
	/**
	* Compute element phase from geometry.
	*
	* @param {PhasedArray} pa
	* */
	compute_illumination(pa){
		const feed_x = this.x;
		const feed_y = this.y;
		const feed_z = this.z;

		const x = pa.geometry.x;
		const y = pa.geometry.y;
		const cx = pa.geometry.x_center;
		const cy = pa.geometry.y_center;
		const rs = Float32Array.from({length: pa.geometry.length}, (_, i) => Math.sqrt((feed_x - x[i] + cx)**2 + (feed_y - y[i] + cy)**2 + feed_z**2));
		const min_r = Math.min(...rs);
		for (let i = 0; i < pa.geometry.length; i++){
			pa.vectorMagIllum[i] = min_r / rs[i];
			pa.vectorPhaseIllumFactor[i] = rs[i];
		}
	};
}

export const PhasedArrayFeeds = [
	FeedUniform,
	FeedReflectArray,
]
