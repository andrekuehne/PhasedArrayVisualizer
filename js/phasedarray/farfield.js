import {linspace} from "../util.js";
import {PATTERN_COS_N} from "./element.js";
import {applyElementPattern, extractPatternMetrics, getFarfieldKernel} from "../wasm/init.js";
import {farfieldPoolSize, runFarfieldJob} from "../wasm/farfield-pool.js";

/**
 * @typedef {FarfieldSpherical | FarfieldUV | FarfieldLudwig3} FarfieldHint
 */

const CON_FREQ = {"title": "freq. scale", "type": "float", "default": 1.0, "min": 0.0, "max": 200, "step": 0.02};
const CON_POINTS = {'title': "Points", 'type': "int", 'default': 257, 'min': 11, 'max': 2049};
const DOMAIN_SPHERICAL = 0;
const DOMAIN_UV = 1;
const DOMAIN_LUDWIG3 = 2;
const TILE_ROWS = 32;

/**
 * @param {ArrayLike<number>} a
 * @returns {Float32Array}
 */
function as_f32(a){
	return a instanceof Float32Array ? a : Float32Array.from(a);
}

export class FarfieldABC{
	// One Points control is applied to both mesh axes.
	static args = ['farfield-points', 'farfield-points', 'farfield-frequency'];
	static controls = {
		'farfield-points': CON_POINTS,
		'farfield-frequency': CON_FREQ,
	};
	constructor(ax1Points, ax2Points, frequencyScale){
		ax1Points = Number(ax1Points)
		ax2Points = Number(ax2Points)
		// ensure samples are odd
		if (ax1Points % 2 == 0) ax1Points++;
		if (ax2Points % 2 == 0) ax2Points++;

		this.farfield_total = new Array(ax2Points);
		this.farfield_log = new Array(ax2Points);
		this.maxValue = -Infinity;
		this.dirMax = null;
		this.idealDirectivity = null;
		this.patternMetrics = null;

		this.frequency_scale = 1.0;
		for (let i = 0; i < ax2Points; i++){
			this.farfield_total[i] = new Float32Array(ax1Points);
			this.farfield_log[i] = new Float32Array(ax1Points);
		}
		this.meshPoints = [ax1Points, ax2Points];
		this.frequencyScale = frequencyScale;
	}
	get domain(){ return this.constructor.domain; };
	_yield(text){
		this.ac++;
		return {
			text: text,
			progress: this.ac,
			max: this.maxProgress
		};
	}
	reset_parameters(){
		this.maxValue = -Infinity;
		const [p1, p2] = this.meshPoints;

		this.farfield_im = new Array(p2);
		this.farfield_re = new Array(p2);

		for (let i = 0; i < p2; i++){
			this.farfield_im[i] = new Float32Array(p1);
			this.farfield_re[i] = new Float32Array(p1);
		}
	}
	clear_parameters(){
		const [p1, p2] = this.meshPoints;
		for (let i2 = 0; i2 < p2; i2++){
			for (let i1 = 0; i1 < p1; i1++){
				this.farfield_im[i2][i1] = 0;
				this.farfield_re[i2][i1] = 0;
			}
		}
	}
	calculate_total(total){
		const [p1, p2] = this.meshPoints;
		for (let i2 = 0; i2 < p2; i2++){
			for (let i1 = 0; i1 < p1; i1++){
				const c = Math.abs(this.farfield_re[i2][i1]**2 + this.farfield_im[i2][i1]**2)/total;
				this.farfield_total[i2][i1] = c;
			}
			this.maxValue = Math.max(this.maxValue, ...this.farfield_total[i2]);
		}
		delete this.farfield_im;
		delete this.farfield_re;
	}
	calculate_log(){
		const [p1, p2] = this.meshPoints;
		for (let i2 = 0; i2 < p2; i2++){
			for (let i1 = 0; i1 < p1; i1++){
				this.farfield_log[i2][i1] = 10*Math.log10(this.farfield_total[i2][i1]/this.maxValue);
			}
		}
	}
	create_parameters(pa){
		let ac = 0;
		const maxProgress = pa.geometry.x.length + 4;
		let [pha, mag] = pa.create_farfield_vectors(this.frequencyScale)
		return {
			yield: (text) => {
				ac++;
				return {
					text: text,
					progress: ac,
					max: maxProgress
				}
			},
			x: pa.geometry.x,
			y: pa.geometry.y,
			pha: pha,
			mag: mag,
			elementKind: pa.elementPattern?.kind ?? 0,
			elementN: pa.elementPattern?.n ?? 0,
		}
	}
	/**
	 * Wrap a flat row-major intensity buffer as `farfield_total[i2][i1]`.
	 * @param {Float32Array} flat
	 */
	wrap_flat_total(flat){
		const [p1, p2] = this.meshPoints;
		this._totalFlat = flat;
		this.farfield_total = new Array(p2);
		for (let i2 = 0; i2 < p2; i2++){
			this.farfield_total[i2] = flat.subarray(i2 * p1, (i2 + 1) * p1);
		}
	}
	/**
	 * Multiply AF intensity by the element power pattern and refresh `maxValue`.
	 * @param {number} domain
	 * @param {ArrayLike<number>} ax1
	 * @param {ArrayLike<number>} ax2
	 * @param {ReturnType<FarfieldABC['create_parameters']>} pars
	 */
	apply_element_pattern(domain, ax1, ax2, pars){
		if (this._totalFlat == null) return;
		const kind = Number(pars.elementKind) || 0;
		const n = Number(pars.elementN) || 0;
		if (kind !== PATTERN_COS_N) return;
		this.maxValue = applyElementPattern(
			domain,
			as_f32(ax1),
			as_f32(ax2),
			this._totalFlat,
			kind,
			n
		);
	}
	/**
	 * Copy WASM `PatternMetrics` into a plain object and free the wasm handle.
	 * @param {import('../wasm/simd/farfield_kernel.js').PatternMetrics} m
	 */
	static copy_pattern_metrics(m){
		const o = {
			peak_i1: m.peak_i1,
			peak_i2: m.peak_i2,
			peak_ax1: m.peak_ax1,
			peak_ax2: m.peak_ax2,
			hpbw_ax1: m.hpbw_ax1,
			hpbw_ax2: m.hpbw_ax2,
			hpbw_ax1_deg: m.hpbw_ax1_deg,
			hpbw_ax2_deg: m.hpbw_ax2_deg,
			hpbw_ax1_clipped: m.hpbw_ax1_clipped,
			hpbw_ax2_clipped: m.hpbw_ax2_clipped,
			hpbw_large: m.hpbw_large,
			hpbw_small: m.hpbw_small,
			hpbw_large_deg: m.hpbw_large_deg,
			hpbw_small_deg: m.hpbw_small_deg,
			hpbw_large_clipped: m.hpbw_large_clipped,
			hpbw_small_clipped: m.hpbw_small_clipped,
			hpbw_large_angle_deg: m.hpbw_large_angle_deg,
			hpbw_small_angle_deg: m.hpbw_small_angle_deg,
			nearest_sll_db: m.nearest_sll_db,
			largest_sll_db: m.largest_sll_db,
			nearest_sll_ax1: m.nearest_sll_ax1,
			nearest_sll_ax2: m.nearest_sll_ax2,
			largest_sll_ax1: m.largest_sll_ax1,
			largest_sll_ax2: m.largest_sll_ax2,
			peak_theta_deg: m.peak_theta_deg,
			peak_phi_deg: m.peak_phi_deg,
			requested_theta_deg: m.requested_theta_deg,
			requested_phi_deg: m.requested_phi_deg,
			squint_deg: m.squint_deg,
			squint_ax1_deg: m.squint_ax1_deg,
			squint_ax2_deg: m.squint_ax2_deg,
		};
		if (typeof m.free === 'function') m.free();
		return o;
	}
	/**
	 * Extract HPBW / SLL / pointing from the current intensity grid.
	 * @param {number} domain
	 * @param {ArrayLike<number>} ax1
	 * @param {ArrayLike<number>} ax2
	 * @param {{theta?: number, phi?: number}} [pa]
	 */
	compute_pattern_metrics(domain, ax1, ax2, pa){
		this.patternMetrics = null;
		if (this._totalFlat == null) return;
		const deg2rad = Math.PI / 180;
		const reqTheta = (pa && Number.isFinite(pa.theta)) ? pa.theta * deg2rad : 0;
		const reqPhi = (pa && Number.isFinite(pa.phi)) ? pa.phi * deg2rad : 0;
		try {
			this.patternMetrics = FarfieldABC.copy_pattern_metrics(
				extractPatternMetrics(domain, as_f32(ax1), as_f32(ax2), this._totalFlat, reqTheta, reqPhi)
			);
		}
		catch (err){
			console.warn('Pattern metrics failed.', err);
		}
	}
	/**
	 * Tile the WASM array-factor kernel over ax2 rows (main thread).
	 * @param {ReturnType<FarfieldABC['create_parameters']>} pars
	 * @param {number} domain
	 * @param {ArrayLike<number>} ax1
	 * @param {ArrayLike<number>} ax2
	 */
	*accumulate_wasm_sync(pars, domain, ax1, ax2){
		const kernel = getFarfieldKernel();
		const [n1, n2] = this.meshPoints;
		const tiles = Math.max(1, Math.ceil(n2 / TILE_ROWS));
		yield {text: 'Preparing farfield...', progress: 0, max: tiles + 2};
		kernel.prepare(n1, n2);
		kernel.set_inputs(
			as_f32(pars.x),
			as_f32(pars.y),
			as_f32(pars.mag),
			as_f32(pars.pha),
			as_f32(ax1),
			as_f32(ax2)
		);
		for (let t = 0; t < tiles; t++){
			const row0 = t * TILE_ROWS;
			const rowCount = Math.min(TILE_ROWS, n2 - row0);
			yield {text: 'Calculating farfield...', progress: t + 1, max: tiles + 2};
			kernel.accumulate_tile(domain, Number(this.frequencyScale), row0, rowCount);
		}
		yield {text: 'Calculating total...', progress: tiles + 1, max: tiles + 2};
		this.maxValue = kernel.finalize(pars.x.length);
		this.wrap_flat_total(kernel.take_total());
		this.apply_element_pattern(domain, ax1, ax2, pars);
	}
	/**
	 * Tile the WASM array-factor kernel over ax2 rows, using workers when available.
	 * @param {ReturnType<FarfieldABC['create_parameters']>} pars
	 * @param {number} domain
	 * @param {ArrayLike<number>} ax1
	 * @param {ArrayLike<number>} ax2
	 */
	async *accumulate_wasm(pars, domain, ax1, ax2){
		const [, n2] = this.meshPoints;
		if (farfieldPoolSize() >= 2 && n2 >= 2){
			try {
				yield* this.accumulate_wasm_workers(pars, domain, ax1, ax2);
				return;
			}
			catch (err){
				if (err && err.message === 'cancelled') return;
				console.warn('Farfield workers failed; using main thread.', err);
			}
		}
		yield* this.accumulate_wasm_sync(pars, domain, ax1, ax2);
	}
	/**
	 * @param {ReturnType<FarfieldABC['create_parameters']>} pars
	 * @param {number} domain
	 * @param {ArrayLike<number>} ax1
	 * @param {ArrayLike<number>} ax2
	 */
	async *accumulate_wasm_workers(pars, domain, ax1, ax2){
		const events = [];
		let notify = null;
		const push = (ev) => {
			events.push(ev);
			if (notify){
				const n = notify;
				notify = null;
				n();
			}
		};
		const wait = () => {
			if (events.length) return Promise.resolve();
			return new Promise((resolve) => { notify = resolve; });
		};

		runFarfieldJob({
			domain,
			frequencyScale: Number(this.frequencyScale),
			nElements: pars.x.length,
			x: as_f32(pars.x),
			y: as_f32(pars.y),
			mag: as_f32(pars.mag),
			pha: as_f32(pars.pha),
			ax1: as_f32(ax1),
			ax2: as_f32(ax2),
			tileRows: TILE_ROWS,
		}, (done, total) => {
			push({kind: 'progress', done, total});
		}).then((result) => {
			push({kind: 'done', result});
		}, (err) => {
			push({kind: 'error', err});
		});

		yield {text: 'Preparing farfield...', progress: 0, max: 1};
		while (true){
			await wait();
			while (events.length){
				const ev = events.shift();
				if (ev.kind === 'progress'){
					yield {
						text: 'Calculating farfield...',
						progress: ev.done,
						max: Math.max(ev.total, 1) + 1,
					};
				}
				else if (ev.kind === 'done'){
					this.maxValue = ev.result.maxValue;
					this.wrap_flat_total(ev.result.total);
					this.apply_element_pattern(domain, ax1, ax2, pars);
					yield {text: 'Calculating total...', progress: 1, max: 1};
					return;
				}
				else {
					throw ev.err || new Error('Farfield worker job failed');
				}
			}
		}
	}
	cut(xc, xs, ys, axis){
		xc = Number(xc);
		const mp = Float32Array.from(xs, (x) => Math.abs(x - xc));
		let mv = Infinity;
		let mi = -1;

		for (let i = 0; i < mp.length; i++){
			if (mp[i] < mv){
				mv = mp[i];
				mi = i;
			}
		}
		if (mi < 0) return null;
		if (axis == 0) return ys[mi]
		return Float32Array.from(xs, (_, i) => ys[i][mi])
	}
}

export class FarfieldSpherical extends FarfieldABC{
	static title = 'Spherical';
	static domain = 'spherical';
	constructor(thetaPoints, phiPoints, frequencyScale){
		super(thetaPoints, phiPoints, frequencyScale);
		[thetaPoints, phiPoints] = this.meshPoints;
		this.thetaPoints = thetaPoints;
		this.phiPoints = phiPoints;

		this.theta = linspace(-Math.PI/2, Math.PI/2, this.thetaPoints);
		this.phi = linspace(-Math.PI/2, Math.PI/2, this.phiPoints);
	}
	async *calculator_loop(pa){
		const pars = this.create_parameters(pa);
		this.idealDirectivity = 4 * Math.PI * pa.geometry.area * this.frequencyScale ** 2;
		yield* this.accumulate_wasm(pars, DOMAIN_SPHERICAL, this.theta, this.phi);
		yield pars.yield('Calculating spherical directivity...');
		this.dirMax = this.compute_directivity();
		yield pars.yield('Calculating spherical pattern metrics...');
		this.compute_pattern_metrics(DOMAIN_SPHERICAL, this.theta, this.phi, pa);
		yield pars.yield('Calculating spherical log...');
		this.calculate_log();
	}
	compute_directivity(){
		let bsa = 0;
		const step = Math.PI/(this.thetaPoints - 1)*Math.PI/(this.phiPoints - 1);
		for (let it = 0; it < this.thetaPoints; it++){
			let st = Math.abs(Math.sin(this.theta[it]))*step;
			for (let ip = 0; ip < this.phiPoints; ip++){
				bsa += this.farfield_total[ip][it]*st;
			}
		}
		return 4 * Math.PI * this.maxValue / bsa;
	}
	constant_phi(phi){
		const y = this.cut(Number(phi)*Math.PI/180, this.phi, this.farfield_log, 0);
		if (y === null) return [null, null];
		return [this.theta, y]
	}
	constant_theta(theta){
		const y = this.cut(Number(theta)*Math.PI/180, this.theta, this.farfield_log, 1);
		if (y === null) return [null, null];
		return [this.phi, y]
	}
}

export class FarfieldUV extends FarfieldABC{
	static title = 'UV';
	static domain = 'uv';
	static args = [...FarfieldABC.args, 'farfield-uv-bound'];
	static controls = {
		...FarfieldABC.controls,
		'farfield-domain': {'title': null},
		'farfield-uv-bound': {'title': "U/V Bound", 'type': "float", 'default': 1, 'min': 0.1, 'step': 0.1},
	};
	constructor(uPoints, vPoints, frequencyScale, uMax, vMax){
		super(uPoints, vPoints, frequencyScale);
		[uPoints, vPoints] = this.meshPoints;
		if (uMax === undefined) uMax = 1;
		if (vMax === undefined) vMax = uMax;
		this.uPoints = uPoints;
		this.vPoints = vPoints;
		this.u = linspace(-uMax, uMax, this.uPoints);
		this.v = linspace(-vMax, vMax, this.vPoints);
	}
	async *calculator_loop(pa){
		const pars = this.create_parameters(pa);
		this.idealDirectivity = 4 * Math.PI * pa.geometry.area * this.frequencyScale ** 2;
		yield* this.accumulate_wasm(pars, DOMAIN_UV, this.u, this.v);
		yield pars.yield('Calculating UV directivity...');
		this.dirMax = this.compute_directivity();
		yield pars.yield('Calculating UV pattern metrics...');
		this.compute_pattern_metrics(DOMAIN_UV, this.u, this.v, pa);
		yield pars.yield('Calculating UV log...');
		this.calculate_log();
	}
	compute_directivity(){
		// Front-hemisphere solid angle: dΩ = du dv / √(1-u²-v²) for u²+v² < 1.
		const du = this.u[1] - this.u[0];
		const dv = this.v[1] - this.v[0];
		let bsa = 0;
		for (let iv = 0; iv < this.vPoints; iv++){
			const v = this.v[iv];
			for (let iu = 0; iu < this.uPoints; iu++){
				const r2 = this.u[iu] * this.u[iu] + v * v;
				if (r2 >= 1) continue;
				bsa += this.farfield_total[iv][iu] / Math.sqrt(1 - r2);
			}
		}
		bsa *= du * dv;
		return 4 * Math.PI * this.maxValue / bsa;
	}
	constant_u(u){
		const y = this.cut(u, this.u, this.farfield_log, 1);
		if (y === null) return [null, null];
		return [this.v, y]
	}
	constant_v(v){
		const y = this.cut(v, this.v, this.farfield_log, 0);
		if (y === null) return [null, null];
		return [this.u, y]
	}
}

export class FarfieldLudwig3 extends FarfieldABC{
	static domain = 'ludwig3';
	static title = 'Ludwig3';
	static controls = {
		...FarfieldABC.controls,
		'farfield-domain': {'title': null},
	};
	constructor(azPoints, elPoints, frequencyScale, azMax, elMax){
		super(azPoints, elPoints, frequencyScale);
		[azPoints, elPoints] = this.meshPoints;
		if (azMax === undefined) azMax = 90;
		if (elMax === undefined) elMax = 90;
		this.azPoints = azPoints;
		this.elPoints = elPoints;
		const sc = Math.PI/180
		this.az = linspace(-azMax*sc, azMax*sc, this.azPoints);
		this.el = linspace(-elMax*sc, elMax*sc, this.elPoints);
	}
	async *calculator_loop(pa){
		const pars = this.create_parameters(pa);
		this.idealDirectivity = 4 * Math.PI * pa.geometry.area * this.frequencyScale ** 2;
		yield* this.accumulate_wasm(pars, DOMAIN_LUDWIG3, this.az, this.el);
		yield pars.yield('Calculating Ludwig3 directivity...');
		this.dirMax = this.compute_directivity();
		yield pars.yield('Calculating Ludwig3 pattern metrics...');
		this.compute_pattern_metrics(DOMAIN_LUDWIG3, this.az, this.el, pa);
		yield pars.yield('Calculating Ludwig3 log...');
		this.calculate_log();
	}
	compute_directivity(){
		// Front-hemisphere solid angle: dΩ = |cos(el)| daz del.
		const daz = this.az[1] - this.az[0];
		const del = this.el[1] - this.el[0];
		const step = daz * del;
		let bsa = 0;
		for (let ie = 0; ie < this.elPoints; ie++){
			const w = Math.abs(Math.cos(this.el[ie])) * step;
			for (let ia = 0; ia < this.azPoints; ia++){
				bsa += this.farfield_total[ie][ia] * w;
			}
		}
		return 4 * Math.PI * this.maxValue / bsa;
	}
	constant_az(az){
		const y = this.cut(az, this.az, this.farfield_log, 1);
		if (y === null) return [null, null];
		return [this.el, y]
	}
	constant_el(el){
		const y = this.cut(el, this.el, this.farfield_log, 0);
		if (y === null) return [null, null];
		return [this.az, y]
	}
}

/** Stub Domain entries: plot Z/S matrices, not a farfield mesh. */
export class FarfieldMatrixZ{
	static title = 'Z';
	static domain = 'z';
	static args = ['farfield-frequency'];
	static controls = {
		'farfield-frequency': CON_FREQ,
	};
}

export class FarfieldMatrixS{
	static title = 'S';
	static domain = 's';
	static args = ['farfield-frequency'];
	static controls = {
		'farfield-frequency': CON_FREQ,
	};
}

export const FarfieldDomains = [
	FarfieldSpherical,
	FarfieldUV,
	FarfieldLudwig3,
	FarfieldMatrixZ,
	FarfieldMatrixS,
]
