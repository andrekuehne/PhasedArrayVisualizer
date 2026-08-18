import {ScenePlotABC} from "../scene-plot-abc.js"
import {GeometryViews, GeometryViewConversions} from "../../phasedarray/geometry-views.js"
import {SceneControlWithSelector} from "../scene-abc.js"
/** @import { SceneControlPhasedArray } from "../../index-scenes.js"*/
/** @import { PhasedArray } from "../../phasedarray/phasedarray.js"*/
/** @import { ElementPhase } from "../../phasedarray/geometry-views.js"*/

export class ScenePlot2DGeometryABC extends ScenePlotABC{
	constructor(parent, canvas, cmapKey, defaultCMAP, strokeColor){
		if (defaultCMAP === undefined) defaultCMAP = 'viridis';
		if (strokeColor === undefined) strokeColor = 'black'
		let cmap = parent.create_mesh_colormap_selector(cmapKey, defaultCMAP);
		super(parent, canvas, cmap);
		this.strokeColor = strokeColor;
		this.cmap.addEventListener('change', () => {this.build_queue();})
		this.create_hover_items();
		this.create_progress_bar();
		this.create_queue();
		this._popup_args = null;
		this._popup_installed = false;
	}
	event_to_id(e){
		let i = null;
		if (e.isTrusted && this.pa !== undefined){
			const f = this.canvas.index_from_event;
			if (f !== undefined) i = f(e);
		}
		return i;
	}
	install_hover_item(callback){
		this.canvas.addEventListener('mousemove', (e) => {
			if (this.queue.running) return;
			let i = this.event_to_id(e);
			let text = "&nbsp;";
			if (i !== null){
				const geo = this.pa.geometry;
				const t = callback(i);
				text = `Element[${i}] (${geo.x[i].toFixed(2)}, ${geo.y[i].toFixed(2)}): ${t}`
			}
			this.canvas.hover_container.innerHTML = text;
		});
	}
	get isValid(){ return (this.pa !== undefined && this.pa !== null); }
	/**
	* Load Phased array object.
	*
	* @param {PhasedArray} pa
	*
	* @return {null}
	* */
	load_phased_array(pa){ this.pa = pa; }
	/**
	* Bind a Phased Array Scene.
	*
	* @param {SceneControlPhasedArray} scene
	*
	* @return {null}
	* */
	bind_phased_array_scene(scene){
		this.arrayControl = scene;
		scene.addEventListener('phased-array-changed', (pa) => this.load_phased_array(pa));
	}
	draw(data){
		if (this.pa === undefined) return;
		const canvas = this.canvas;
		const colormap = this.cmap.cmap();
		const scale = 600;
		const geo = this.pa.geometry;
		canvas.width = scale;
		canvas.height = scale;
		const ctx = canvas.getContext('2d');
		this.cmap.changed = false;
		ctx.reset();

		const maxX = Math.max(...geo.x) + geo.dx/2;
		const minX = Math.min(...geo.x) - geo.dx/2;
		const maxY = Math.max(...geo.y) + geo.dy/2;
		const minY = Math.min(...geo.y) - geo.dy/2;

		const wx = (maxX - minX);
		const wy = (maxY - minY);
		const sc = Math.min(canvas.width/wx, canvas.height/wy);
		const ox = (canvas.width - wx*sc)/2 - minX*sc;
		const oy = (canvas.height - wy*sc)/2 - minY*sc;
		const dx = geo.dx*0.98*sc/2;
		const dy = geo.dy*0.98*sc/2;

		const _xy_to_wh = (x, y) => [x*sc+ox, scale-(y*sc+oy)];
		const _wh_to_xy = (w, h) => [(w-ox)/sc, (scale-h-oy)/sc];

		canvas.transform_to_xy = _wh_to_xy;
		canvas.transform_to_wh = _xy_to_wh;
		canvas.index_from_event = (e) => {
			const rect = canvas.getBoundingClientRect();
			let ex, ey;
			if (e.type == 'touchstart'){
				ex = e.touches[0].clientX;
				ey = e.touches[0].clientY;
			}
			else{
				ex = e.clientX;
				ey = e.clientY;
			}
			const wx = (ex - rect.left)/rect.width*canvas.width;
			const wy = (ey - rect.top)/rect.height*canvas.height;
			const [x, y] = _wh_to_xy(wx, wy);
			const dx = geo.dx/2;
			const dy = geo.dy/2;
			let eleI = null;
			for (let i = 0; i < geo.length; i++){
				if (((x - geo.x[i])/dx)**2 + ((y - geo.y[i])/dy)**2 <= 1) {
					eleI = i;
					break;
				}
			}
			return eleI;
		};
		for (let i = 0; i < geo.length; i++){
			let [x, y] = _xy_to_wh(geo.x[i], geo.y[i]);
			ctx.beginPath();
			ctx.ellipse(x, y, dx, dy, 0.0, 0.0, 2*Math.PI);
			ctx.closePath();
			ctx.fillStyle = colormap(data[i]);
			ctx.fill();
			ctx.strokeStyle = this.strokeColor;
			ctx.stroke();
		}
	}
	build_queue(){ throw Error("Don't call generic build_queue."); }
	uninstall_popup(){
		this._popup_args = null;
	}
	install_popup(dtype, controls, changedCallback, updaterCallback, clearAllCallback){
		this._popup_args = [dtype, controls, changedCallback, updaterCallback, clearAllCallback];
		if (this._popup_installed) return;
		this._popup_installed = true;

		const _show_popup = (e) => {
			if (this._popup_args === null) return;
			if (e.pointerType == 'touch' && e.type == 'click') return;
			if (this.queue.running) return;
			if (!this.isValid) return;
			let i = this.event_to_id(e);
			if (i === null) return;
			const pa = this.pa;

			const dcontrols = [{
				'label': "Element ID",
				'type': 'number',
				'min': 0,
				'max': pa.geometry.length - 1,
				'id': 'index',
				'value': i,
			},{
				'label': "Location:",
				'type': 'span',
				'id': 'loc',
			},{
				'label': `Current ${this._popup_args[0]}:`,
				'type': 'span',
				'id': 'current-value',
			},{
				'label': "Enable Override",
				'type': 'checkbox',
				'id': 'override',
				'value': true,
			}];
			this._popup_args[1].forEach((e) => dcontrols.push(e));
			const popup = this.create_popup("Manually Change " + this._popup_args[0], dcontrols, (config) => {
				if (config === null) return;
				this._popup_args[2](_i(), config);
				this.parent.build_queue();
			});
			const lbl = popup.element('loc');
			const _i = () => Math.max(0, Math.min(pa.geometry.length - 1, popup.element('index').value));
			const _update = () => {
				const i = _i();
				const res = this._popup_args[3](i);
				lbl.innerHTML = `(${pa.geometry.x[i].toFixed(2)}, ${pa.geometry.y[i].toFixed(2)})`;
				for (const [key, value] of Object.entries(res)) popup.set_element_value(key, value);
			}
			popup.element('index').addEventListener('change', _update);
			popup.add_action(`Clear All ${this._popup_args[0]} Overrides`).addEventListener('click', () => {
				this._popup_args[4]();
				this.parent.build_queue();
			});
			popup.add_note(
				'To override phase/attenuation, select "Enable Override" '
				+ 'and enter the desired value.', 'popup-note');
			_update();
			popup.show_from_event(e);
		}
		const onlongpress = (ele, cb) => {
			let tid;
			ele.addEventListener('touchstart', (e) => {
				tid = setTimeout(() => {
					tid = null;
					e.stopPropagation();
					cb(e);
				}, 500);
			});
			ele.addEventListener('contextmenu', (e) => { e.preventDefault(); });
			ele.addEventListener('touchend', (e) => { if (tid) clearTimeout(tid); else e.preventDefault();});
			ele.addEventListener('touchmove', () => { if (tid) clearTimeout(tid);});
		}
		this.canvas.addEventListener('click', _show_popup);
		onlongpress(this.canvas, _show_popup);
		this._popup_enabled = true;
	}
}

export class ScenePlot2DGeometryGeneric extends ScenePlot2DGeometryABC{
	constructor(parent, div, defaultView, defaultUnit){
		if (defaultView === undefined) defaultView = "Element";
		if (defaultUnit === undefined) defaultUnit = "deg";
		let cid = div.id.substring(parent.prepend.length + 1);
		let div_title = document.createElement("div");
		let div_canvas = document.createElement("div");
		let div_footer = document.createElement("div");
		let div_group = document.createElement("div");
		let div_sel = document.createElement("div");
		let div_desc = document.createElement("div");
		const header_title = document.createElement("h2");
		let canvas = document.createElement("canvas");
		let sel_cmap = document.createElement("select");
		let sel_sel = document.createElement("select");
		let sel_unit = document.createElement("select");
		let inp_scale = document.createElement("input");

		let title_span = document.createElement("span");
		title_span.innerHTML = "&nbsp;"

		sel_cmap.id = parent.prepend + "-" + cid + "-cmap";
		sel_sel.id = parent.prepend + "-" + cid + "-view";
		sel_unit.id = parent.prepend + "-" + cid + "-unit";
		canvas.className = "canvas-grid";
		div_canvas.className = "canvas-wrapper";
		div_title.className = "canvas-header";
		div_footer.className = "canvas-footer";
		div_group.className = "footer-group";
		div_sel.className = "footer-group";
		div_desc.className = "footer-group";

		div_title.appendChild(header_title)
		div_title.appendChild(title_span)
		header_title.innerHTML = "TEST";
		div_canvas.appendChild(canvas);

		div.appendChild(div_title);
		div.appendChild(div_canvas);
		div.appendChild(div_footer);

		const cdiv = (ele, txt) => {
			let div = document.createElement("div");
			let lbl = document.createElement("label");
			lbl.innerHTML = txt;
			lbl.htmlFor = ele.id;

			div.appendChild(lbl);
			div.appendChild(ele);
			return div;
		}
		div_footer.appendChild(div_sel)
		div_footer.appendChild(div_desc)
		div_footer.appendChild(div_group)
		div_sel.appendChild(cdiv(sel_sel, "View"));
		div_sel.appendChild(cdiv(sel_unit, ""));
		div_group.appendChild(cdiv(sel_cmap, "Colormap"));

		const div_scale = cdiv(inp_scale, "Scale");
		const find_cmap = () => {
			if (this.active_unit === null) return "hsv";
			return this.active_unit.default_cmap();
		}
		div_group.appendChild(div_scale);

		super(parent, canvas, cid + "-cmap", find_cmap);
		this.active_view = null;
		this.active_unit = null;
		this.unit_selector = this.find_element(cid + "-unit");
		this.add_event_types('data-min-changed');
		this.addEventListener('data-min-changed', () => {this.build_queue();})
		this.view_selector = new SceneControlGeometryViews(this, cid + "-view", defaultView);
		this.unit_selector = new SceneControlGeometryViewConversions(this, cid + "-unit", defaultUnit)
		this._needsNewData = false;
		this._plot_data = null;
		inp_scale.type = "number";
		inp_scale.max = "200";
		inp_scale.min = "5";
		inp_scale.value = "40";
		inp_scale.id = div.id + "-scale";

		const update_title = () => {
			let vkls = this.view_selector.selected_class();
			let ukls = this.unit_selector.selected_class();
			this.build_queue();
			header_title.innerHTML = `${vkls.title} ${ukls.user_title}`;
			if (ukls.show_scale) div_scale.style.display = "";
			else div_scale.style.display = "none";
			div_desc.innerHTML = `<i style='font-size: 0.8em;'>${vkls.desc}</i>`;
		}
		this.view_selector.addEventListener('active-class-changed', () => {
			update_title()
		});
		this.unit_selector.addEventListener('active-class-changed', () => {
			update_title()
		});
		this.install_scale_control(inp_scale);
		this.install_hover_item((i) => {
			if (this.active_view === null || this.active_unit == null || this.pa == null) return "";
			let t = this.active_unit.string_from_view(this.active_view, this.pa, i);
			const pwr = (this.arrayControl === undefined) ? undefined : this.arrayControl.powerControl;
			if (pwr !== undefined){
				t += `, ${pwr.format(this.pa.elementPowerWatts(pwr.getWatts(), i))}`;
			}
			return t;
		});
	}
	/**
	* Bind a Phased Array Scene.
	*
	* @param {SceneControlPhasedArray} scene
	*
	* @return {null}
	* */
	bind_phased_array_scene(scene){
		super.bind_phased_array_scene(scene);
		scene.addEventListener('phased-array-calculation-changed', () => {
			this.build_queue();
		});
		this.build_queue();
	}
	build_queue(){
		this.queue.request(() => {
			if (this.view_selector === undefined) return;
			if (this.pa === undefined) return;
			const title = this.view_selector.selected_class().title;
			this.queue.add(`Creating ${title}...`, () => {
				this.active_view = this.view_selector.build_active_object();
			});
			this.queue.add(`Creating ${title} Unit...`, () => {
				this.active_unit = this.unit_selector.build_active_object();
				if (!this.active_view.constructor.allow_manual || this.active_unit.constructor.popup_type == null) this.uninstall_popup();
				else if (this.active_unit.constructor.popup_type == "phase") this.install_phase_popup();
				else if (this.active_unit.constructor.popup_type == "atten") this.install_atten_popup();
			});
			this.queue.add(`Scaling ${title}...`, () => {
				this._plot_data = this.active_unit.convert_from_view(this.active_view, this.pa, this.min);
			});
			this.queue.add(`Drawing ${title}...`, () => {
				this.draw();
			});
			this.queue.start("&nbsp;");
		});
	}
	draw(){ return super.draw(this._plot_data); }
	install_phase_popup(){
		this.install_popup('Phase', [{
			'label': `Manual Phase (deg)`,
			'type': 'number',
			'min': 0,
			'max': 360,
			'id': 'value',
			'value': 0,
			'focus': true,
		}], (i, config) => {
			this.pa.set_manual_phase(i, config['override'], config['value'] * Math.PI/180);
		}, (i) => {
			let ov;
			if (this.pa.vectorPhaseIsManual[i]) ov = this.pa.vectorPhaseManual[i] / (2 * Math.PI);
			else ov = this.pa.vQuantizePhaseFactor[i] % 1.0;
			const nv = (ov * 360).toFixed(2);
			return {
				'value': nv,
				'override': this.pa.vectorPhaseIsManual[i],
				'current-value': `${nv} deg`
			}
		}, () => { this.pa.clear_all_manual_phase();})
	}
	install_atten_popup(){
		this.install_popup('Attenuation', [{
			'label': `Manual Attenuation (dB)`,
			'type': 'number',
			'min': -100,
			'max': 100,
			'id': 'value',
			'step': 'none',
			'value': 0,
			'focus': true,
		},{
			'label': `Disable Element`,
			'type': 'checkbox',
			'id': 'disabled',
		}], (i, config) => {
			this.pa.set_manual_magnitude(i, config['override'], 10**(-Math.abs(config['value'])/20), config['disabled']);
		}, (i) => {
			let ov;
			if (this.pa.vectorMagIsManual[i]) ov = this.pa.vectorMagManual[i];
			else ov = this.pa.vQuantizeMag[i];
			const nv = (20*Math.log10(Math.abs(ov))).toFixed(2);
			return {
				'value': nv,
				'override': this.pa.vectorMagIsManual[i],
				'current-value': `${nv} dB`,
				'disabled': this.pa.elementDisabled[i],
			}
		}, () => { this.pa.clear_all_manual_magnitude();})
	}
}
export class SceneControlGeometryViews extends SceneControlWithSelector{
	static autoUpdateURL = false;
	constructor(parent, key, defaultValue){
		super(parent, key, GeometryViews, undefined, true, defaultValue);
	}
	/**
	* Add callable objects to queue.
	*
	* @param {SceneQueue} queue
	*
	* @return {null}
	* */
	add_to_queue(queue){
		const arrayControl = this.parent.arrayControl;
		let needsRecalc = arrayControl.farfieldNeedsCalculation;

		this.needsRedraw = needsRecalc;
	}
}

export class SceneControlGeometryViewConversions extends SceneControlWithSelector{
	static autoUpdateURL = false;
	constructor(parent, key, defaultValue){
		super(parent, key, GeometryViewConversions, undefined, true, defaultValue);
	}
	/**
	* Add callable objects to queue.
	*
	* @param {SceneQueue} queue
	*
	* @return {null}
	* */
	add_to_queue(queue){
		const arrayControl = this.parent.arrayControl;
		let needsRecalc = arrayControl.farfieldNeedsCalculation;

		this.needsRedraw = needsRecalc;
	}
}
