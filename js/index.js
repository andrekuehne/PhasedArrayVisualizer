import {SceneControlPhasedArray, SceneControlFarfieldDomain, SceneTaperCuts} from "./index-scenes.js";
import {ScenePlotFarfieldCuts} from "./scene/plot-1d/scene-plot-farfield-cuts.js";
import {ScenePlotFarfield2DPlotly} from "./scene/plot-2d/scene-plot-2d-farfield-plotly.js";
import {ScenePlot2DGeometryGeneric} from "./scene/plot-2d/scene-plot-2d-geometry.js";
import {SceneParent} from "./scene/scene-abc.js"
import {SceneTheme} from "./scene/scene-util.js";
import {initFarfieldWasm} from "./wasm/init.js";

document.addEventListener('DOMContentLoaded', async () => {
	new SceneTheme();
	await initFarfieldWasm();
	const scene = new PhasedArrayScene('pa');
	scene.build_queue();
});

/**
 * @param {Record<string, number|boolean>|null} m
 * @param {'ax1'|'ax2'|'large'|'small'} axis
 * @param {boolean} isUv
 */
function formatHpbw(m, axis, isUv){
	if (m == null) return '—';
	const clipped = m[`hpbw_${axis}_clipped`];
	const deg = m[`hpbw_${axis}_deg`];
	if (clipped || !Number.isFinite(deg)) return 'clipped';
	if (isUv){
		const native = m[`hpbw_${axis}`];
		const n = Number.isFinite(native) ? native.toFixed(3) : '—';
		return `${deg.toFixed(1)} deg (${n})`;
	}
	return `${deg.toFixed(1)} deg`;
}

/**
 * @param {number} db
 */
function formatSll(db){
	if (!Number.isFinite(db)) return '—';
	return `${db.toFixed(1)} dB`;
}

/**
 * @param {number} theta
 * @param {number} phi
 */
function formatAnglePair(theta, phi){
	if (!Number.isFinite(theta) || !Number.isFinite(phi)) return '—';
	return `(${theta.toFixed(1)}, ${phi.toFixed(1)}) deg`;
}

/**
 * @param {number} pct
 */
function formatHpbwPercent(pct){
	if (pct >= 9.95) return pct.toFixed(0);
	return pct.toFixed(1);
}

/**
 * Squint in degrees, plus each axis offset as a fraction of that axis HPBW.
 * @param {Record<string, number|boolean>|null} m
 * @param {[string, string]} axisNames
 */
function formatSquint(m, axisNames){
	if (m == null || !Number.isFinite(m.squint_deg)) return '—';
	const deg = `${m.squint_deg.toFixed(2)} deg`;
	const parts = [];
	for (const [i, axis] of ['ax1', 'ax2'].entries()){
		const clipped = m[`hpbw_${axis}_clipped`];
		const hpbw = m[`hpbw_${axis}_deg`];
		if (clipped || !Number.isFinite(hpbw) || hpbw <= 0) continue;
		const along = m[`squint_${axis}_deg`];
		if (!Number.isFinite(along)) continue;
		const pct = formatHpbwPercent(100 * along / hpbw);
		parts.push(`${pct}% ${axisNames[i]}`);
	}
	if (parts.length === 0) return deg;
	return `${deg} (${parts.join(', ')})`;
}

/**	 *
 * Create scene for Phased Array simulator.
 *
 * @param {string} prepend - Prepend used on HTML IDs.
 * */
export class PhasedArrayScene extends SceneParent{
	constructor(prepend){
		super(prepend, ['refresh', 'reset'])
		this.create_queue(this.find_element('progress'), this.find_element('status'));
		this.arrayControl = new SceneControlPhasedArray(this);
		this.farfieldControl = new SceneControlFarfieldDomain(this, 'farfield-domain');

		this.plotFF = new ScenePlotFarfield2DPlotly(this, this.find_element('farfield-plot-2d'), 'farfield-2d-colormap');
		this.plot1D = new ScenePlotFarfieldCuts(this, this.find_element('farfield-canvas-1d'), 'farfield-1d-colormap');
		this.plotTaper = new SceneTaperCuts(this, this.find_element('taper-canvas-1d'), 'taper-1d-colormap');
		this.geoPlot1 = new ScenePlot2DGeometryGeneric(this, this.find_element('geo-canvas-1'), "Element", "deg");
		this.geoPlot2 = new ScenePlot2DGeometryGeneric(this, this.find_element('geo-canvas-2'), "Element", "dB");

		this.geoPlot1.bind_phased_array_scene(this.arrayControl);
		this.geoPlot2.bind_phased_array_scene(this.arrayControl);
		this.plot1D.bind_farfield_scene(this.farfieldControl);
		this.plotFF.bind_farfield_scene(this.farfieldControl);
		this.plotTaper.bind_phased_array_scene(this.arrayControl);
		this.plot1D.install_scale_control('farfield-1d-scale');
		this.plotFF.install_scale_control('farfield-2d-scale');
		const update_power_displays = () => {
			const pa = this.arrayControl.pa;
			const pwr = this.arrayControl.powerControl;
			if (pa == null) return;
			const pAnt = pwr.getWatts();
			const pAvail = pa.availablePowerWatts(pAnt);
			const pStim = pa.totalPowerWatts(pAnt);
			const matched = pa.coupling === 'matched';
			const gamma = matched ? Number(pa.reflectionGamma) || 0 : 0;
			const pRef = matched ? pStim * gamma : 0;
			const pAcc = matched ? pStim * (1 - gamma) : pStim;
			this.find_element('available-power').innerHTML = pwr.formatEirp(pAvail, 1);
			this.find_element('stimulated-power').innerHTML = pwr.formatEirp(pStim, 1);
			this.find_element('reflected-power').innerHTML = matched ? pwr.formatEirp(pRef, 1) : '—';
			this.find_element('accepted-power').innerHTML = matched ? pwr.formatEirp(pAcc, 1) : '—';
			const util = pAvail > 0 ? pStim / pAvail : 0;
			const utilDb = util > 0 ? (10 * Math.log10(util)).toFixed(1) : '-Inf';
			this.find_element('power-utilization').innerHTML = `${(util * 100).toFixed(1)} % (${utilDb} dB)`;
			const ff = this.farfieldControl.ff;
			if (ff == null || ff.dirMax == null) return;
			const eirp = ff.dirMax * pAcc;
			this.find_element('peak-eirp').innerHTML = pwr.formatEirp(eirp, 1);
			this.plotFF.set_title_metrics(
				`Directivity: ${(10*Math.log10(ff.dirMax)).toFixed(1)} dB, EIRP: ${pwr.formatEirp(eirp, 1)}`
			);
		};
		this.farfieldControl.add_max_monitor('directivity', (v) => {
			let idir = this.farfieldControl.ff.idealDirectivity;
			this.find_element('calc-directivity').innerHTML = `${(10*Math.log10(v)).toFixed(1)} dB`
			this.find_element('aperture-efficiency').innerHTML = `${((v / idir) * 100).toFixed(1)} % (${(10 * Math.log10(v / idir)).toFixed(1)} dB)`
			update_power_displays();
		});
		this.farfieldControl.add_max_monitor('ideal-directivity', (v) => {
			this.find_element('ideal-directivity').innerHTML = `${(10 * Math.log10(v)).toFixed(1)} dB`
		});
		this.farfieldControl.add_max_monitor('pattern-metrics', (m) => {
			this.update_pattern_metrics_display(m);
		});
		this.arrayControl.powerControl.addEventListener('power-changed', () => {
			update_power_displays();
		});

		this.arrayControl.addEventListener("phased-array-calculation-changed", () => {
			this.find_element('calc-area').innerHTML = `${this.arrayControl.pa.geometry.area.toFixed(1)} λ<sub>0</sub><sup>2</sup>`
			this.find_element('element-count').innerHTML = `${this.arrayControl.pa.geometry.length}`
			update_power_displays();
		});
		this.find_element('refresh').addEventListener('click', () => {
			this.update_url_parameters();
			this.build_queue();
		});
		this.find_element('reset').addEventListener('click', () => {
			this.reset_url_parameters();
			this.build_queue();
		});
		this.create_popup_overlay();
		this.bind_url_elements();
		this.trigger_event('scene-loaded');
	}
	update_pattern_metrics_display(m){
		const ff = this.farfieldControl.ff;
		const labels = {
			spherical: ['HPBW Theta', 'HPBW Phi'],
			uv: ['HPBW U', 'HPBW V'],
			ludwig3: ['HPBW Az', 'HPBW El'],
		};
		const axisNames = {
			spherical: ['Theta', 'Phi'],
			uv: ['U', 'V'],
			ludwig3: ['Az', 'El'],
		};
		const pair = (ff && labels[ff.domain]) ? labels[ff.domain] : ['HPBW Axis 1', 'HPBW Axis 2'];
		const names = (ff && axisNames[ff.domain]) ? axisNames[ff.domain] : ['Axis 1', 'Axis 2'];
		this.find_element('hpbw-ax1-label').textContent = pair[0];
		this.find_element('hpbw-ax2-label').textContent = pair[1];
		const isUv = ff != null && ff.domain === 'uv';
		this.find_element('hpbw-ax1').textContent = formatHpbw(m, 'ax1', isUv);
		this.find_element('hpbw-ax2').textContent = formatHpbw(m, 'ax2', isUv);
		this.find_element('hpbw-large').textContent = formatHpbw(m, 'large', isUv);
		this.find_element('hpbw-small').textContent = formatHpbw(m, 'small', isUv);
		this.find_element('nearest-sll').textContent = formatSll(m && m.nearest_sll_db);
		this.find_element('largest-sll').textContent = formatSll(m && m.largest_sll_db);
		this.find_element('requested-angle').textContent = formatAnglePair(
			m && m.requested_theta_deg, m && m.requested_phi_deg
		);
		this.find_element('achieved-angle').textContent = formatAnglePair(
			m && m.peak_theta_deg, m && m.peak_phi_deg
		);
		this.find_element('squint-angle').textContent = formatSquint(m, names);
	}
	build_queue(){
		this.queue.request(() => {
			this.arrayControl.add_to_queue(this.queue);
			this.farfieldControl.add_to_queue(this.queue);
			this.queue.start();
		});
	}
}
