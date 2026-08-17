import {adjust_theta_phi, linspace, rad2deg} from "../../util.js";
import {SceneObjectABC} from "../scene-abc.js"
/** @import { FarfieldHint } from "../../phasedarray/farfield.js" */
/** @import { SceneControlFarfieldDomain } from "../../index-scenes.js" */
/** @import { SceneParent } from "../scene-abc.js" */

const PLOTLY_CMAPS = ['Viridis', 'Inferno', 'Plasma', 'Hot', 'Jet', 'Rainbow', 'Turbo'];
const GRID_LINE = 'rgba(255, 255, 255, 0.2)';
const BORDER_LINE = 'rgba(255, 255, 255, 0.4)';
const PLOT_CONFIG = {
	responsive: true,
	displayModeBar: true,
	displaylogo: false,
	scrollZoom: true,
};

/**
 * 2-D farfield plot using Plotly heatmap. Canvas engines in
 * scene-plot-2d-farfield.js are left in place and unused.
 */
export class ScenePlotFarfield2DPlotly extends SceneObjectABC{
	/**
	 * @param {SceneParent} parent
	 * @param {HTMLElement} plotDiv
	 * @param {string} cmapKey
	 */
	constructor(parent, plotDiv, cmapKey){
		super(parent.prepend, []);
		parent.addEventListener("scene-loaded", () => {this.trigger_event("scene-loaded")});
		parent.add_child(this);
		this.parent = parent;
		this.el = plotDiv;
		this.add_event_types('data-min-changed');
		this.min = -40;
		this.ff = null;
		this._ready = false;
		this._title = '';
		this._meshKey = '';
		this._hoverBound = false;
		this._sphN = 0;
		this._sphZ = null;
		this._sphX = null;
		this._sphY = null;

		this._init_colormap(cmapKey);
		this.addEventListener('data-min-changed', () => {this._restyle_scale();});
		this._cmapSelect.addEventListener('change', () => {this._restyle_colorscale();});
		if (typeof window.installThemeChanged === 'function'){
			window.installThemeChanged(() => {this._relayout_theme();});
		}

		this._hover = document.createElement('div');
		this._hover.className = 'ff-plotly-hover';
		this.el.parentElement.appendChild(this._hover);

		this._resizeObs = new ResizeObserver(() => {
			if (this._ready && window.Plotly) window.Plotly.Plots.resize(this.el);
		});
		this._resizeObs.observe(this.el);
	}
	install_scale_control(key){
		const ele = this.find_element(key);
		const _val = () => {
			const v = -Math.max(5, Math.abs(ele.value));
			ele.value = Math.abs(v);
			return v;
		}
		ele.addEventListener('change', () => {
			this.min = _val();
			this.trigger_event('data-min-changed', this.min);
		})
		this.min = _val();
	}
	/**
	 * @param {SceneControlFarfieldDomain} scene
	 */
	bind_farfield_scene(scene){
		scene.addEventListener('farfield-calculation-complete', (ff) => {
			this.load_farfield(ff);
			this.update();
		});
	}
	/**
	 * @param {FarfieldHint} ff
	 */
	load_farfield(ff){
		this.ff = ff;
	}
	get isValid(){ return this.ff != null; }
	/**
	 * Update Plotly title with Directivity / EIRP (no z transfer).
	 * @param {string} text
	 */
	set_title_metrics(text){
		this._title = text;
		if (this._ready) this._apply_title();
	}
	/**
	 * Draw or restyle the heatmap from the current farfield.
	 */
	update(){
		if (!this.isValid || !window.Plotly) return;
		const ff = this.ff;
		const key = this._mesh_key(ff);
		if (!this._ready || key !== this._meshKey){
			this._meshKey = key;
			this._react(ff);
			return;
		}
		const z = this._z_data(ff);
		window.Plotly.restyle(this.el, {z: [z]}, [0]);
	}
	_init_colormap(cmapKey){
		const sel = this.find_element(cmapKey);
		this._cmapSelect = sel;
		sel.replaceChildren();
		const add = (name) => {
			const opt = document.createElement('option');
			opt.value = name;
			opt.textContent = name;
			sel.appendChild(opt);
		};
		PLOTLY_CMAPS.forEach((name) => {
			add(name);
			add(name + '_r');
		});
		sel.value = 'Viridis';
		sel.setAttribute('data-default-value', 'Viridis');
	}
	_cmap_style(){
		const v = this._cmapSelect.value || 'Viridis';
		const rev = v.endsWith('_r');
		return {
			colorscale: rev ? v.slice(0, -2) : v,
			reversescale: rev,
		};
	}
	_mesh_key(ff){
		let extra = '';
		if (ff.domain === 'uv') extra = `:${ff.u[0]}:${ff.u[ff.u.length - 1]}`;
		return `${ff.domain}:${ff.meshPoints[0]}:${ff.meshPoints[1]}${extra}`;
	}
	_theme_colors(){
		const s = getComputedStyle(document.documentElement);
		return {
			bg: s.getPropertyValue('--cell-bg').trim() || '#16161f',
			fg: s.getPropertyValue('--text-color').trim() || '#AAA',
		};
	}
	_axis_base(){
		return {
			showticklabels: false,
			ticks: '',
			zeroline: false,
			automargin: false,
		};
	}
	_layout(ff){
		const {bg, fg} = this._theme_colors();
		const layout = {
			title: {
				text: this._title,
				font: {size: 13, color: fg},
				x: 0.5,
				xanchor: 'center',
				y: 0.98,
				yanchor: 'top',
			},
			margin: {t: 48, l: 8, r: 8, b: 8, pad: 0},
			paper_bgcolor: bg,
			plot_bgcolor: bg,
			font: {color: fg},
			dragmode: 'zoom',
			hovermode: 'closest',
			showlegend: false,
			autosize: true,
			uirevision: 'farfield-2d',
		};
		if (ff.domain === 'spherical') this._layout_spherical(layout);
		else if (ff.domain === 'uv') this._layout_uv(layout, ff);
		else this._layout_ludwig3(layout, ff);
		return layout;
	}
	_layout_spherical(layout){
		const ax = {
			...this._axis_base(),
			range: [-1, 1],
			showgrid: false,
			showline: false,
			constrain: 'domain',
		};
		layout.xaxis = {...ax};
		layout.yaxis = {...ax, scaleanchor: 'x', scaleratio: 1};
		layout.shapes = this._spherical_shapes();
	}
	_spherical_shapes(){
		const phiSteps = 13;
		const thetaSteps = 7;
		const shapes = [];
		const start = 1 / (thetaSteps - 1);
		const c = 2 * Math.PI / (phiSteps - 1);
		for (let i = 0; i < phiSteps - 1; i++){
			const ph = i * c;
			const cs = Math.cos(ph);
			const sn = Math.sin(ph);
			shapes.push({
				type: 'line',
				x0: cs * start, y0: sn * start,
				x1: cs, y1: sn,
				line: {color: GRID_LINE, width: 1, dash: 'dash'},
				layer: 'above',
			});
		}
		const rc = 1 / (thetaSteps - 1);
		for (let i = 1; i < thetaSteps - 1; i++){
			const r = i * rc;
			shapes.push({
				type: 'circle',
				xref: 'x', yref: 'y',
				x0: -r, y0: -r, x1: r, y1: r,
				line: {color: GRID_LINE, width: 1, dash: 'dash'},
				fillcolor: 'rgba(0,0,0,0)',
				layer: 'above',
			});
		}
		return shapes;
	}
	_rect_grid_shapes(x0, x1, y0, y1, xSteps, ySteps, circleR){
		const shapes = [];
		const gridLine = {color: GRID_LINE, width: 1, dash: 'dash'};
		for (let i = 1; i < xSteps - 1; i++){
			const x = x0 + (x1 - x0) * i / (xSteps - 1);
			shapes.push({
				type: 'line',
				xref: 'x', yref: 'y',
				x0: x, x1: x, y0, y1,
				line: gridLine,
				layer: 'above',
			});
		}
		for (let i = 1; i < ySteps - 1; i++){
			const y = y0 + (y1 - y0) * i / (ySteps - 1);
			shapes.push({
				type: 'line',
				xref: 'x', yref: 'y',
				x0, x1, y0: y, y1: y,
				line: gridLine,
				layer: 'above',
			});
		}
		shapes.push({
			type: 'rect',
			xref: 'x', yref: 'y',
			x0, x1, y0, y1,
			line: {color: BORDER_LINE, width: 1, dash: 'solid'},
			fillcolor: 'rgba(0,0,0,0)',
			layer: 'above',
		});
		if (circleR != null){
			shapes.push({
				type: 'circle',
				xref: 'x', yref: 'y',
				x0: -circleR, y0: -circleR, x1: circleR, y1: circleR,
				line: {color: BORDER_LINE, width: 1, dash: 'dash'},
				fillcolor: 'rgba(0,0,0,0)',
				layer: 'above',
			});
		}
		return shapes;
	}
	_layout_uv(layout, ff){
		const u0 = ff.u[0];
		const u1 = ff.u[ff.u.length - 1];
		const v0 = ff.v[0];
		const v1 = ff.v[ff.v.length - 1];
		const ax = {
			...this._axis_base(),
			showgrid: false,
			showline: false,
			zeroline: false,
			constrain: 'domain',
		};
		layout.xaxis = {...ax, range: [u0, u1]};
		layout.yaxis = {
			...ax,
			range: [v0, v1],
			scaleanchor: 'x',
			scaleratio: 1,
		};
		layout.shapes = this._rect_grid_shapes(u0, u1, v0, v1, 11, 11, 1);
	}
	_layout_ludwig3(layout, ff){
		const az0 = ff.az[0] * 180 / Math.PI;
		const az1 = ff.az[ff.az.length - 1] * 180 / Math.PI;
		const el0 = ff.el[0] * 180 / Math.PI;
		const el1 = ff.el[ff.el.length - 1] * 180 / Math.PI;
		const ax = {
			...this._axis_base(),
			showgrid: false,
			showline: false,
			zeroline: false,
			constrain: 'domain',
		};
		layout.xaxis = {...ax, range: [az0, az1]};
		layout.yaxis = {
			...ax,
			range: [el0, el1],
			scaleanchor: 'x',
			scaleratio: 1,
		};
		layout.shapes = this._rect_grid_shapes(az0, az1, el0, el1, 13, 13);
	}
	_trace(ff){
		const cmap = this._cmap_style();
		const xy = this._xy(ff);
		return {
			type: 'heatmap',
			x: xy.x,
			y: xy.y,
			z: this._z_data(ff),
			zmin: this.min,
			zmax: 0,
			zauto: false,
			zsmooth: false,
			showscale: false,
			hoverinfo: 'none',
			hoverongaps: false,
			colorscale: cmap.colorscale,
			reversescale: cmap.reversescale,
		};
	}
	_xy(ff){
		if (ff.domain === 'spherical'){
			this._ensure_spherical_grid(ff);
			return {x: this._sphX, y: this._sphY};
		}
		if (ff.domain === 'uv') return {x: ff.u, y: ff.v};
		return {x: rad2deg(ff.az), y: rad2deg(ff.el)};
	}
	_z_data(ff){
		if (ff.domain === 'spherical') return this._project_spherical(ff);
		return ff.farfield_log;
	}
	_ensure_spherical_grid(ff){
		const n = Math.max(ff.thetaPoints, ff.phiPoints);
		if (this._sphN === n && this._sphZ != null) return;
		this._sphN = n;
		this._sphX = linspace(-1, 1, n);
		this._sphY = linspace(-1, 1, n);
		this._sphZ = new Array(n);
		for (let i = 0; i < n; i++) this._sphZ[i] = new Float32Array(n);
	}
	_project_spherical(ff){
		this._ensure_spherical_grid(ff);
		const n = this._sphN;
		const x = this._sphX;
		const y = this._sphY;
		const z = this._sphZ;
		const log = ff.farfield_log;
		const nTheta = ff.thetaPoints;
		const nPhi = ff.phiPoints;
		const thetaStep = Math.PI / (nTheta - 1);
		const phiStep = Math.PI / (nPhi - 1);
		const itMax = nTheta - 1;
		const ipMax = nPhi - 1;
		for (let iy = 0; iy < n; iy++){
			const v = y[iy];
			const row = z[iy];
			for (let ix = 0; ix < n; ix++){
				const u = x[ix];
				const r = Math.hypot(u, v);
				if (r > 1){
					row[ix] = NaN;
					continue;
				}
				const [th, ph] = adjust_theta_phi(r * Math.PI / 2, Math.atan2(v, u), false);
				let itf = (Math.PI / 2 + th) / thetaStep;
				let ipf = (Math.PI / 2 + ph) / phiStep;
				if (itf < 0) itf = 0;
				else if (itf > itMax) itf = itMax;
				if (ipf < 0) ipf = 0;
				else if (ipf > ipMax) ipf = ipMax;
				const it0 = Math.floor(itf);
				const ip0 = Math.floor(ipf);
				const it1 = it0 < itMax ? it0 + 1 : it0;
				const ip1 = ip0 < ipMax ? ip0 + 1 : ip0;
				const ft = itf - it0;
				const fp = ipf - ip0;
				const v00 = log[ip0][it0];
				const v10 = log[ip0][it1];
				const v01 = log[ip1][it0];
				const v11 = log[ip1][it1];
				row[ix] = (1 - ft) * ((1 - fp) * v00 + fp * v01) + ft * ((1 - fp) * v10 + fp * v11);
			}
		}
		return z;
	}
	_react(ff){
		const Plotly = window.Plotly;
		const data = [this._trace(ff)];
		const layout = this._layout(ff);
		const done = () => {
			this._ready = true;
			this._bind_hover();
			this._apply_title();
		};
		if (!this._ready){
			Plotly.newPlot(this.el, data, layout, PLOT_CONFIG).then(done);
		}
		else{
			Plotly.react(this.el, data, layout, PLOT_CONFIG).then(done);
		}
	}
	_restyle_scale(){
		if (!this._ready || !window.Plotly) return;
		window.Plotly.restyle(this.el, {zmin: this.min, zmax: 0}, [0]);
	}
	_restyle_colorscale(){
		if (!this._ready || !window.Plotly) return;
		const cmap = this._cmap_style();
		window.Plotly.restyle(this.el, {
			colorscale: cmap.colorscale,
			reversescale: cmap.reversescale,
		}, [0]);
	}
	_apply_title(){
		if (!this._ready || !window.Plotly) return;
		window.Plotly.relayout(this.el, {'title.text': this._title});
	}
	_relayout_theme(){
		if (!this._ready || !window.Plotly) return;
		const {bg, fg} = this._theme_colors();
		window.Plotly.relayout(this.el, {
			paper_bgcolor: bg,
			plot_bgcolor: bg,
			'font.color': fg,
			'title.font.color': fg,
		});
	}
	_bind_hover(){
		if (this._hoverBound) return;
		this._hoverBound = true;
		this.el.on('plotly_hover', (ev) => {
			const pt = ev.points && ev.points[0];
			const text = this._hover_text(pt);
			if (!text){
				this._hover.style.display = 'none';
				return;
			}
			this._hover.innerHTML = text;
			this._hover.style.display = 'block';
			const wrap = this.el.parentElement;
			const rect = wrap.getBoundingClientRect();
			const mx = ev.event ? ev.event.clientX : rect.left;
			const my = ev.event ? ev.event.clientY : rect.top;
			const pad = 12;
			let left = mx - rect.left + pad;
			let top = my - rect.top + pad;
			const hw = this._hover.offsetWidth;
			const hh = this._hover.offsetHeight;
			if (left + hw + 4 > rect.width) left = mx - rect.left - hw - pad;
			if (top + hh + 4 > rect.height) top = my - rect.top - hh - pad;
			if (left < 4) left = 4;
			if (top < 4) top = 4;
			this._hover.style.left = `${left}px`;
			this._hover.style.top = `${top}px`;
		});
		this.el.on('plotly_unhover', () => {
			this._hover.style.display = 'none';
		});
	}
	_nearest_index(axis, value){
		if (axis == null || axis.length === 0) return 0;
		if (axis.length === 1) return 0;
		const step = axis[1] - axis[0];
		let i = (step === 0) ? 0 : Math.round((value - axis[0]) / step);
		if (i < 0) i = 0;
		if (i >= axis.length) i = axis.length - 1;
		return i;
	}
	_hover_text(pt){
		const ff = this.ff;
		if (pt == null || ff == null || ff.dirMax == null) return '';
		if (pt.x == null || pt.y == null) return '';
		let it;
		let ip;
		let x;
		let y;
		let xLabel;
		let yLabel;
		let unit = '';
		if (ff.domain === 'spherical'){
			const u = Number(pt.x);
			const v = Number(pt.y);
			const r = Math.hypot(u, v);
			if (r > 1) return '';
			const thetaStep = Math.PI / (ff.thetaPoints - 1);
			const phiStep = Math.PI / (ff.phiPoints - 1);
			const [th, ph] = adjust_theta_phi(r * Math.PI / 2, Math.atan2(v, u), false);
			it = Math.round((Math.PI / 2 + th) / thetaStep);
			ip = Math.round((Math.PI / 2 + ph) / phiStep);
			if (it >= ff.thetaPoints) it = ff.thetaPoints - 1;
			if (ip >= ff.phiPoints) ip = ff.phiPoints - 1;
			if (it < 0) it = 0;
			if (ip < 0) ip = 0;
			x = ff.theta[it] * 180 / Math.PI;
			y = ff.phi[ip] * 180 / Math.PI;
			xLabel = 'θ';
			yLabel = 'φ';
			unit = '°';
		}
		else if (ff.domain === 'uv'){
			it = this._nearest_index(ff.u, Number(pt.x));
			ip = this._nearest_index(ff.v, Number(pt.y));
			x = ff.u[it];
			y = ff.v[ip];
			xLabel = 'u';
			yLabel = 'v';
		}
		else{
			it = this._nearest_index(ff.az, Number(pt.x) * Math.PI / 180);
			ip = this._nearest_index(ff.el, Number(pt.y) * Math.PI / 180);
			x = ff.az[it] * 180 / Math.PI;
			y = ff.el[ip] * 180 / Math.PI;
			xLabel = 'Az';
			yLabel = 'El';
			unit = '°';
		}
		const row = ff.farfield_total[ip];
		if (row == null) return '';
		const total = row[it];
		if (total == null || ff.maxValue <= 0) return '';
		const ff1 = 10 * Math.log10(total / ff.maxValue);
		const ff2 = ff1 + 10 * Math.log10(ff.dirMax);
		const unitSuffix = unit ? ` ${unit}` : '';
		let text = `${xLabel} = ${x.toFixed(2)}${unitSuffix}, ${yLabel} = ${y.toFixed(2)}${unitSuffix}`;
		text += `<br>Directivity: ${ff2.toFixed(2)} dBi (${ff1.toFixed(2)} dB)`;
		const arrayControl = this.parent.arrayControl;
		const pwr = (arrayControl === undefined) ? undefined : arrayControl.powerControl;
		const pa = (arrayControl === undefined) ? null : arrayControl.pa;
		if (pwr !== undefined && pa != null){
			const eirp = (total / ff.maxValue) * ff.dirMax * pa.totalPowerWatts(pwr.getWatts());
			text += `<br>EIRP: ${pwr.formatEirp(eirp)}`;
		}
		return text;
	}
}
