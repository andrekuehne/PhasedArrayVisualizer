import {SceneControl, SceneControlWithSelector, SceneControlWithSelectorAutoBuild, SceneParent} from "./scene/scene-abc.js";
import {FindSceneURL} from "./scene/scene-util.js";
import {ScenePlot1D} from "./scene/plot-1d/scene-plot-1d.js";
import {Geometries} from "./phasedarray/geometry.js";
import {PhasedArray} from "./phasedarray/phasedarray.js";
import {FarfieldDomains} from "./phasedarray/farfield.js"
import {GeometryViews} from "./phasedarray/geometry-views.js"
import {SteeringDomains} from "./phasedarray/steering.js"
import {Illuminations} from "./phasedarray/illumination.js"
import {ElementCosN, ElementGreenPec, ElementTypes, exponentFromPeakDbi, MIN_ELEMENT_GAIN_DBI, PATTERN_GREEN_PEC} from "./phasedarray/element.js"
import {Tapers} from "./phasedarray/tapers.js"
import {defaultWatts, formatPower, formatPowerValue, isPowerScope, isPowerUnit, wattsFrom} from "./phasedarray/power.js"
import {GREEN_PEC_AUTO_ISOLATE_N, MATCHED_AUTO_ISOLATE_N, nMuFromGeometry, Z_REF} from "./phasedarray/matched.js"
import {getRadiatedPowerKernel, zSelfPecDipole} from "./wasm/init.js"
import {linspace} from "./util.js";
/** @import { SceneQueue } from "./scene/scene-queue.js" */

export class SceneControlGeometry extends SceneControlWithSelectorAutoBuild{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, 'geometry', Geometries, parent.find_element('geometry-controls'));
		this.activeGeometry = null;
	}
	control_changed(key){
		super.control_changed(key);
		this.activeGeometry = null;
	}
	get calculationWaiting(){ return this.activeGeometry === null; }
	add_to_queue(queue){
		if (this.calculationWaiting){
			queue.add('Building geometry...', () => {
					this.activeGeometry = this.build_active_object();
					this.activeGeometry.build();
				}
			)
		}
	}
}

const CHANGE_ILLUM	= 1 << 0;
const CHANGE_PHASE 	= 1 << 1;
const CHANGE_ATTEN 	= 1 << 2;
const CHANGE_PHASEQ = 1 << 3;
const CHANGE_ATTENQ = 1 << 4;
const CHANGE_PA 	= 1 << 5;

export class SceneControlPhasedArray extends SceneControl{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, ['phase-bits', 'atten-lsb', 'atten-bits', 'atten-manual', 'phase-manual', 'phase-dither', 'coupling', 'coupling-z0-re', 'coupling-z0-im', 'coupling-xnn', 'coupling-alpha', 'coupling-beta', 'coupling-aniso', 'coupling-eps-x', 'coupling-eps-y', 'coupling-att', 'steer-law']);
		this.pa = null;
		this._matchedFreq = NaN;
		this._matchedKind = null;
		this._matchedN = null;
		this._matchedModel = null;
		this._matchedXnn = NaN;
		this._matchedAlpha = NaN;
		this._matchedBeta = NaN;
		this._matchedAniso = NaN;
		this._matchedAtt = NaN;
		this._matchedEpsX = NaN;
		this._matchedEpsY = NaN;
		this._matchedZ0Re = NaN;
		this._matchedZ0Im = NaN;
		this._matchedH = NaN;
		this._matchedEll = NaN;
		this._matchedA = NaN;
		this._matchedGreenZ0 = NaN;
		this._steerFreq = NaN;
		this.geometryControl = new SceneControlGeometry(this);
		this.taperControl = new SceneControlAllTapers(this);
		this.steerControl = new SceneControlSteeringDomain(this);
		this.illumControl = new SceneControlIllumination(this);
		this.elementControl = new SceneControlElement(this);
		this.powerControl = new SceneControlPower(this);
		this._steerLawRadios = {
			geometric: document.getElementById(this.prepend + '-steer-law-geometric'),
			conjugate: document.getElementById(this.prepend + '-steer-law-conjugate'),
		};
		this._bind_mode_radios('steer-law', this._steerLawRadios);
		this.sync_mode_radios();
		this.addEventListener('scene-loaded', () => {
			this.migrate_coupling_url();
			this.sync_mode_radios();
		});
		this.addEventListener('reset', () => {
			this.find_element('coupling').value = 'isolated';
			this.find_element('steer-law').value = 'geometric';
			this.find_element('coupling-xnn').value = '0';
			this.find_element('coupling-alpha').value = '2';
			this.find_element('coupling-beta').value = '0';
			this.find_element('coupling-aniso').value = '0';
			this.find_element('coupling-eps-x').value = '1';
			this.find_element('coupling-eps-y').value = '1';
			this.find_element('coupling-att').value = '0';
			this.find_element('coupling-z0-re').value = '45';
			this.find_element('coupling-z0-im').value = '5';
			this.sync_mode_radios();
		});
		this.add_event_types(
			'phased-array-changed',
			'phased-array-phase-changed',
			'phased-array-attenuation-changed',
			'phased-array-calculation-changed',
		);
	}
	_bind_mode_radios(key, radios){
		for (const radio of Object.values(radios)){
			if (radio === null) continue;
			radio.addEventListener('change', () => {
				if (!radio.checked) return;
				this.find_element(key).value = radio.value;
				this.find_element(key).dispatchEvent(new Event('change'));
			});
		}
	}
	sync_mode_radios(){
		const law = this.steerLaw();
		for (const [value, radio] of Object.entries(this._steerLawRadios || {})){
			if (radio) radio.checked = value === law;
		}
		const couplingDiv = document.getElementById(this.prepend + '-coupling-div');
		const z0ReDiv = document.getElementById(this.prepend + '-coupling-z0-re-div');
		const z0ImDiv = document.getElementById(this.prepend + '-coupling-z0-im-div');
		const xnnDiv = document.getElementById(this.prepend + '-coupling-xnn-div');
		const alphaDiv = document.getElementById(this.prepend + '-coupling-alpha-div');
		const betaDiv = document.getElementById(this.prepend + '-coupling-beta-div');
		const anisoDiv = document.getElementById(this.prepend + '-coupling-aniso-div');
		const epsXDiv = document.getElementById(this.prepend + '-coupling-eps-x-div');
		const epsYDiv = document.getElementById(this.prepend + '-coupling-eps-y-div');
		const attDiv = document.getElementById(this.prepend + '-coupling-att-div');
		const hideMatch = () => {
			if (couplingDiv) couplingDiv.style.display = 'none';
			if (z0ReDiv) z0ReDiv.style.display = 'none';
			if (z0ImDiv) z0ImDiv.style.display = 'none';
			if (xnnDiv) xnnDiv.style.display = 'none';
			if (alphaDiv) alphaDiv.style.display = 'none';
			if (betaDiv) betaDiv.style.display = 'none';
			if (anisoDiv) anisoDiv.style.display = 'none';
			if (epsXDiv) epsXDiv.style.display = 'none';
			if (epsYDiv) epsYDiv.style.display = 'none';
			if (attDiv) attDiv.style.display = 'none';
		};
		if (this.isGreenPec()){
			hideMatch();
			return;
		}
		const model = this.matchModel();
		const matched = model !== 'isolated';
		const powerLaw = model === 'per-port' || model === 'common';
		const prop = model === 'propagation';
		const showZ0 = model === 'common' || prop;
		const showX = matched;
		if (couplingDiv) couplingDiv.style.display = 'flex';
		if (z0ReDiv) z0ReDiv.style.display = showZ0 ? 'flex' : 'none';
		if (z0ImDiv) z0ImDiv.style.display = showX ? 'flex' : 'none';
		if (xnnDiv) xnnDiv.style.display = showX ? 'flex' : 'none';
		if (alphaDiv) alphaDiv.style.display = powerLaw ? 'flex' : 'none';
		if (betaDiv) betaDiv.style.display = powerLaw ? 'flex' : 'none';
		if (anisoDiv) anisoDiv.style.display = powerLaw ? 'flex' : 'none';
		if (epsXDiv) epsXDiv.style.display = prop ? 'flex' : 'none';
		if (epsYDiv) epsYDiv.style.display = prop ? 'flex' : 'none';
		if (attDiv) attDiv.style.display = prop ? 'flex' : 'none';
	}
	/** Isolated | per-port | common | propagation (legacy matched+match-style migrated on load). */
	matchModel(){
		const v = this.find_element('coupling').value;
		if (v === 'per-port' || v === 'common' || v === 'propagation') return v;
		return 'isolated';
	}
	isGreenPec(){
		return this.elementControl != null && this.elementControl.selected_class() === ElementGreenPec;
	}
	/**
	 * Map legacy URL coupling=matched (+ match-style) onto the Matching select.
	 */
	migrate_coupling_url(){
		const ele = this.find_element('coupling');
		const url = FindSceneURL();
		const raw = url.get_param('coupling');
		let v = raw != null ? raw : ele.value;
		if (v === 'matched'){
			const style = url.get_param('match-style');
			v = style === 'per-port' ? 'per-port' : 'common';
		}
		else if (v !== 'per-port' && v !== 'common' && v !== 'propagation'){
			v = 'isolated';
		}
		if (ele.value !== v) ele.value = v;
		if (url.get_param('match-style') != null) url.delete('match-style');
		if (typeof this.parent.update_url_parameters === 'function') this.parent.update_url_parameters();
	}
	couplingMode(){
		if (this.isGreenPec()){
			const n = this.pa ? this.pa.size : 0;
			return n > GREEN_PEC_AUTO_ISOLATE_N ? 'isolated' : 'matched';
		}
		return this.matchModel() === 'isolated' ? 'isolated' : 'matched';
	}
	_matrix_domain_selected(){
		const ff = this.parent.farfieldControl;
		return ff != null && typeof ff.isMatrixDomain === 'function' && ff.isMatrixDomain();
	}
	/**
	 * Force Isolated when a rebuilt array is too large for matched S.
	 * A later manual Coupling click is not a rebuild, so it stays matched.
	 * @param {number} n
	 */
	applyCouplingSizeSafeguard(n){
		if (this.isGreenPec()) return;
		if (!(n > MATCHED_AUTO_ISOLATE_N) || this.couplingMode() !== 'matched') return;
		this.find_element('coupling').value = 'isolated';
		this.sync_mode_radios();
		if (typeof this.parent.update_url_parameters === 'function') this.parent.update_url_parameters();
	}
	steerLaw(){
		return this.find_element('steer-law').value === 'conjugate' ? 'conjugate' : 'geometric';
	}
	couplingXnn(){
		const v = Number(this.find_element('coupling-xnn').value);
		return Number.isFinite(v) ? v : 0;
	}
	couplingAlpha(){
		const v = Number(this.find_element('coupling-alpha').value);
		if (!Number.isFinite(v) || v < 0) return 2;
		return v;
	}
	couplingBeta(){
		const v = Number(this.find_element('coupling-beta').value);
		return Number.isFinite(v) ? v : 0;
	}
	couplingAniso(){
		const v = Number(this.find_element('coupling-aniso').value);
		return Number.isFinite(v) ? v : 0;
	}
	couplingEpsX(){
		const v = Number(this.find_element('coupling-eps-x').value);
		return Number.isFinite(v) && v >= 0 ? v : 1;
	}
	couplingEpsY(){
		const v = Number(this.find_element('coupling-eps-y').value);
		return Number.isFinite(v) && v >= 0 ? v : 1;
	}
	couplingAtt(){
		const v = Number(this.find_element('coupling-att').value);
		if (!Number.isFinite(v) || v < 0) return 0;
		return v;
	}
	couplingZ0Re(){
		const v = Number(this.find_element('coupling-z0-re').value);
		return Number.isFinite(v) && v > 0 ? v : Z_REF;
	}
	couplingXSelf(){
		const v = Number(this.find_element('coupling-z0-im').value);
		return Number.isFinite(v) ? v : 0;
	}
	/** Common-Z0 / Propagation kernel args: [z_c, x_self]. Per-port passes z_c = 0 so the solver runs. */
	commonZ0Args(){
		const xSelf = this.couplingXSelf();
		if (this.matchModel() === 'per-port') return [0, xSelf];
		return [this.couplingZ0Re(), xSelf];
	}
	frequencyScale(){
		const ele = this.parent.find_element('farfield-frequency', false);
		const v = ele ? Number(ele.value) : 1;
		return Number.isFinite(v) && v > 0 ? v : 1;
	}
	control_changed(key){
		super.control_changed(key);
		if (key === 'coupling' || key === 'steer-law'
			|| key === 'coupling-xnn' || key === 'coupling-alpha'
			|| key === 'coupling-beta' || key === 'coupling-aniso'
			|| key === 'coupling-eps-x' || key === 'coupling-eps-y' || key === 'coupling-att'
			|| key === 'coupling-z0-re' || key === 'coupling-z0-im'){
			this.sync_mode_radios();
			if (typeof this.parent.update_url_parameters === 'function') this.parent.update_url_parameters();
			this.request_recompute();
		}
	}
	compute_matched_basis(freq){
		const pa = this.pa;
		const ep = pa.elementPattern;
		const kind = ep ? ep.kind : 0;
		const elemN = ep ? ep.n : 0;
		const n = pa.size;
		if (kind === PATTERN_GREEN_PEC){
			if (n > GREEN_PEC_AUTO_ISOLATE_N) return;
			const h = ep.h;
			const ell = ep.ell;
			const a = ep.a;
			const zc = ep.z0;
			const xSelf = ep.xself;
			if (
				pa.tRe && pa.tRe.length === n * n
				&& pa.zRe && pa.zRe.length === n * n
				&& this._matchedFreq === freq
				&& this._matchedKind === kind
				&& this._matchedH === h
				&& this._matchedEll === ell
				&& this._matchedA === a
				&& this._matchedGreenZ0 === zc
				&& this._matchedXSelf === xSelf
			) return;
			const kernel = getRadiatedPowerKernel();
			kernel.form_green_pec_dipole(pa.geometry.x, pa.geometry.y, freq, h, ell, a, Z_REF, zc, xSelf);
			pa.set_matched_basis(
				kernel.take_z0(),
				kernel.take_s_re(),
				kernel.take_s_im(),
				kernel.take_t_re(),
				kernel.take_t_im(),
				kernel.take_z0_im(),
				kernel.take_z_re(),
				kernel.take_z_im()
			);
			this._matchedFreq = freq;
			this._matchedKind = kind;
			this._matchedH = h;
			this._matchedEll = ell;
			this._matchedA = a;
			this._matchedGreenZ0 = zc;
			this._matchedXSelf = xSelf;
			this._matchedModel = 'green-pec';
			return;
		}
		const model = this.matchModel();
		const xnn = this.couplingXnn();
		const alpha = this.couplingAlpha();
		const beta = this.couplingBeta();
		const aniso = this.couplingAniso();
		const att = this.couplingAtt();
		const epsX = this.couplingEpsX();
		const epsY = this.couplingEpsY();
		const [zcRe, xSelf] = this.commonZ0Args();
		if (
			pa.tRe && pa.tRe.length === n * n
			&& pa.zRe && pa.zRe.length === n * n
			&& this._matchedFreq === freq
			&& this._matchedKind === kind
			&& this._matchedN === elemN
			&& this._matchedModel === model
			&& this._matchedXnn === xnn
			&& this._matchedAlpha === alpha
			&& this._matchedBeta === beta
			&& this._matchedAniso === aniso
			&& this._matchedAtt === att
			&& this._matchedEpsX === epsX
			&& this._matchedEpsY === epsY
			&& this._matchedZ0Re === zcRe
			&& this._matchedZ0Im === xSelf
		) return;
		const kernel = getRadiatedPowerKernel();
		const nMu = nMuFromGeometry(pa.geometry, freq);
		kernel.set_quadrature(nMu, 2);
		kernel.compute_j0(pa.geometry.x, pa.geometry.y, freq, kind, elemN);
		if (model === 'propagation'){
			kernel.form_matched_s_propagation(
				Z_REF, pa.geometry.x, pa.geometry.y, xnn, att, epsX, epsY, freq, zcRe, xSelf
			);
		}
		else{
			kernel.form_matched_s(Z_REF, pa.geometry.x, pa.geometry.y, xnn, alpha, beta * freq, aniso, zcRe, xSelf);
		}
		pa.set_matched_basis(
			kernel.take_z0(),
			kernel.take_s_re(),
			kernel.take_s_im(),
			kernel.take_t_re(),
			kernel.take_t_im(),
			kernel.take_z0_im(),
			kernel.take_z_re(),
			kernel.take_z_im()
		);
		this._matchedFreq = freq;
		this._matchedKind = kind;
		this._matchedN = elemN;
		this._matchedModel = model;
		this._matchedXnn = xnn;
		this._matchedAlpha = alpha;
		this._matchedBeta = beta;
		this._matchedAniso = aniso;
		this._matchedAtt = att;
		this._matchedEpsX = epsX;
		this._matchedEpsY = epsY;
		this._matchedZ0Re = zcRe;
		this._matchedZ0Im = xSelf;
	}
	/**
	* Add callable objects to queue.
	*
	* @param {SceneQueue} queue
	*
	* @return {null}
	* */
	add_to_queue(queue){
		let changeFlag = 0;
		this.farfieldNeedsCalculation = false
		this.geometryControl.add_to_queue(queue);
		this.taperControl.add_to_queue(queue);
		this.illumControl.add_to_queue(queue);
		this.elementControl.add_to_queue(queue);

		if (this.steerControl.calculationWaiting) changeFlag |= CHANGE_PHASE;
		if (this.taperControl.calculationWaiting) changeFlag |= CHANGE_ATTEN;
		if (this.changed['phase-bits'] || this.changed['phase-dither']) changeFlag |= CHANGE_PHASEQ;
		if (this.changed['atten-bits'] || this.changed['atten-lsb']) changeFlag |= CHANGE_ATTENQ;
		if (this.changed['steer-law']) changeFlag |= CHANGE_PHASE;
		if (this.changed['coupling']
			|| this.changed['coupling-xnn'] || this.changed['coupling-alpha']
			|| this.changed['coupling-beta'] || this.changed['coupling-aniso']
			|| this.changed['coupling-eps-x'] || this.changed['coupling-eps-y']
			|| this.changed['coupling-att']
			|| this.changed['coupling-z0-re'] || this.changed['coupling-z0-im']){
			changeFlag |= CHANGE_PHASE;
			this.farfieldNeedsCalculation = true;
		}
		const freq = this.frequencyScale();
		if (this.steerLaw() === 'conjugate' && this._steerFreq !== freq) changeFlag |= CHANGE_PHASE;

		if (this.elementControl.calculationWaiting) this.farfieldNeedsCalculation = true;

		if (this.geometryControl.calculationWaiting || this.pa === null){
			queue.add('Updating array...', () => {
					let first = this.pa === null;
					this.pa = new PhasedArray(this.geometryControl.activeGeometry);
					this.applyCouplingSizeSafeguard(this.pa.size);
					this.pa.coupling = this.couplingMode();
					this.trigger_event('phased-array-changed', this.pa);
					if (first) this.load_hidden_controls();
				}
			)
			changeFlag |= CHANGE_ILLUM | CHANGE_ATTEN | CHANGE_PHASE;
			this.farfieldNeedsCalculation = true;
		}
		if (this.elementControl.calculationWaiting || this.geometryControl.calculationWaiting || this.pa === null){
			queue.add('Setting element pattern...', () => {
				this.pa.elementPattern = this.elementControl.activeElement;
			});
		}
		if (this.pa !== null && this.pa.requestUpdate) changeFlag |= CHANGE_PA;
		if (this.illumControl.calculationWaiting || (changeFlag & CHANGE_ILLUM)){
			changeFlag |= CHANGE_ATTEN | CHANGE_PHASE;
			queue.add('Computing Illumination...', () => {
				this.pa.set_illumination_type(this.illumControl.activeIllumination);
				this.pa.compute_illumination();
			});
		}
		const greenOversize = this.isGreenPec()
			&& this.pa != null
			&& this.pa.size > GREEN_PEC_AUTO_ISOLATE_N
			&& !this.geometryControl.calculationWaiting;
		const coupling = this.couplingMode();
		const wantsMatch = !greenOversize && (coupling === 'matched' || this._matrix_domain_selected());
		if (wantsMatch){
			const n = this.pa ? this.pa.size : 0;
			const hasT = this.pa && this.pa.tRe && this.pa.tRe.length === n * n;
			const hasZ = this.pa && this.pa.zRe && this.pa.zRe.length === n * n;
			const [zcRe, xSelf] = this.commonZ0Args();
			const basisDirty = this.pa === null
				|| !hasT
				|| !hasZ
				|| this.geometryControl.calculationWaiting
				|| this.elementControl.calculationWaiting
				|| this._matchedFreq !== freq
				|| this._matchedModel !== (this.isGreenPec() ? 'green-pec' : this.matchModel())
				|| (!this.isGreenPec() && (
					this._matchedXnn !== this.couplingXnn()
					|| this._matchedAlpha !== this.couplingAlpha()
					|| this._matchedBeta !== this.couplingBeta()
					|| this._matchedAniso !== this.couplingAniso()
					|| this._matchedAtt !== this.couplingAtt()
					|| this._matchedEpsX !== this.couplingEpsX()
					|| this._matchedEpsY !== this.couplingEpsY()
					|| this._matchedZ0Re !== zcRe
					|| this._matchedZ0Im !== xSelf
				));
			if (basisDirty){
				queue.add('Computing matched S...', () => {
					if (this.isGreenPec()){
						if (this.pa.size > GREEN_PEC_AUTO_ISOLATE_N) return;
						this.compute_matched_basis(freq);
						return;
					}
					if (this.couplingMode() !== 'matched' && !this._matrix_domain_selected()) return;
					this.compute_matched_basis(freq);
				});
				if (coupling === 'matched'){
					changeFlag |= CHANGE_PHASE;
					this.farfieldNeedsCalculation = true;
				}
			}
		}
		if (changeFlag & CHANGE_PHASE){
			queue.add('Calculating phase...', () => {
				const [theta, phi] = this.steerControl.get_theta_phi();
				this.pa.set_theta_phi(theta, phi);
				this.pa.coupling = this.couplingMode();
				if (this.steerLaw() === 'conjugate') this.pa.compute_conjugate_phase(freq);
				else this.pa.compute_phase();
				this._steerFreq = freq;
				this.clear_changed('steer-law', 'coupling', 'coupling-xnn', 'coupling-alpha', 'coupling-beta', 'coupling-aniso', 'coupling-eps-x', 'coupling-eps-y', 'coupling-att', 'coupling-z0-re', 'coupling-z0-im');
			});
		}
		if (changeFlag & CHANGE_ATTEN){
			this.taperControl.add_calculator_queue(queue, this);
		}
		if (changeFlag){
			queue.add('Calculating vector...', () => {
				this.pa.calculate_requested_vector();
				this.update_hidden_controls();
			});
			changeFlag |= CHANGE_PHASEQ | CHANGE_ATTENQ;
		}
		if (changeFlag & CHANGE_PHASEQ){
			queue.add('Quantizing phase...', () => {
				const bits = Math.max(0, Math.min(10, this.find_element('phase-bits').value));
				const dither = this.find_element('phase-dither').checked;
				this.pa.quantize_phase(bits, dither);
				this.trigger_event('phased-array-phase-changed', this.pa);
				this.clear_changed('phase-bits', 'phase-dither');
			});
			this.farfieldNeedsCalculation = true;
		}
		if (changeFlag & CHANGE_ATTENQ){
			queue.add('Quantizing attenuation...', () => {
				const bits = Math.max(0, Math.min(10, this.find_element('atten-bits').value));
				const lsb = Math.max(0, Math.min(5, this.find_element('atten-lsb').value));
				this.pa.quantize_attenuation(bits, lsb);
				this.trigger_event('phased-array-attenuation-changed', this.pa);
				this.clear_changed('atten-bits', 'atten-lsb');
			});
			this.farfieldNeedsCalculation = true;
		}
		if (changeFlag){
			queue.add('Calculating farfield vector change...', () => {
				this.trigger_event('phased-array-calculation-changed', this.pa);
			});
		}
	}
	update_hidden_controls(){
		const mconfig = {};
		const pconfig = {};
		const mele = this.find_element('atten-manual');
		const pele = this.find_element('phase-manual');
		const pa = this.pa;
		for (let i = 0; i < pa.size; i++){
			if (!pa.vectorMagIsManual[i]) continue;
			mconfig[i] = [pa.vectorMagManual[i], pa.elementDisabled[i]];
		}
		if (Object.keys(mconfig).length === 0) mele.value = "";
		else mele.value = JSON.stringify(mconfig);
		mele.dispatchEvent(new Event('change'))
		for (let i = 0; i < pa.size; i++){
			if (!pa.vectorPhaseIsManual[i]) continue;
			pconfig[i] = pa.vectorPhaseManual[i];
		}
		if (Object.keys(pconfig).length === 0) pele.value = "";
		else pele.value = JSON.stringify(pconfig);
		pele.dispatchEvent(new Event('change'))
		const url = FindSceneURL();
		url.check_element('atten-manual', mele);
		url.check_element('phase-manual', pele);
	}
	load_hidden_controls(){
		const pa = this.pa;
		const mele = this.find_element('atten-manual');
		const pele = this.find_element('phase-manual');
		try{
			if (mele.value != ""){
				const mconfig = JSON.parse(mele.value);
				for (let i = 0; i < pa.size; i++){
					if (mconfig[i] === undefined) continue;
					const [v, d] = mconfig[i];
					pa.set_manual_magnitude(i, true, v, d);
				}
			}
		}
		catch(error){ console.log(error); }
		try{
			if (pele.value != ""){
				const pconfig = JSON.parse(pele.value);
				for (let i = 0; i < pa.size; i++){
					if (pconfig[i] === undefined) continue;
					pa.set_manual_phase(i, true, pconfig[i]);
				}
			}
		}
		catch(error){ console.log(error); }
	}
}

export class SceneControlTaper extends SceneControlWithSelectorAutoBuild{
	static autoUpdateURL = false;
	constructor(parent, key, htmlElement){
		super(parent, 'taper', Tapers, htmlElement, key);
		const sel = htmlElement.querySelector('select');
		if (sel !== null && htmlElement.querySelector(`label[for='${sel.id}']`) === null){
			const wrap = document.createElement('div');
			wrap.className = 'form-group';
			wrap.id = sel.id + '-div';
			const lbl = document.createElement('label');
			lbl.htmlFor = sel.id;
			lbl.innerHTML = key.toUpperCase() + '-Taper';
			htmlElement.insertBefore(wrap, sel);
			wrap.appendChild(lbl);
			wrap.appendChild(sel);
			sel.style.width = '';
		}
		this._activeTaper = null;
	}
	control_changed(key){
		super.control_changed(key);
		this._activeTaper = null;
		if (key.endsWith('taper')) this.request_recompute();
	}
	get calculationWaiting(){
		return this._activeTaper === null;
	}
	get activeTaper(){
		if (this._activeTaper === null) this._activeTaper = this.build_active_object();
		return this._activeTaper;
	}
	/**
	* Add callable objects to queue.
	*
	* @param {SceneQueue} queue
	*
	* @return {null}
	* */
	add_to_queue(queue){
		if (this.calculationWaiting){
			queue.add('Building Taper...', () => {
					this._activeTaper = this.build_active_object();
				}
			)
		}
	}
	/**
	* Build a taper control object.
	*
	* @param {SceneParent} parent
	* @param {String} key "x" or "y"
	*
	* @return {SceneControlTaper}
	* */
	static build(parent, key){
		const element = parent.find_element('taper-' + key + '-group')
		const k = parent.prepend + "-" + key + "-taper";
		const _create_group = (p) => {
			let kk = k;
			if (p !== undefined) kk += "-" + p;
			kk += "-div";
			var div = document.createElement('div');
			div.className = 'form-group';
			div.id = kk;
			element.appendChild(div);
			return div;
		}
		const _create_lbl = (div, p) => {
			let kk = k;
			if (p !== undefined) kk += "-" + p;
			const lbl = document.createElement("label");
			lbl.setAttribute("for", kk);
			div.appendChild(lbl);
			return lbl;
		}
		const _create_input = (div, p) => {
			let kk = k;
			if (p !== undefined) kk += "-" + p;
			const inp = document.createElement("input");
			inp.id = kk;
			inp.setAttribute('type', 'Number');
			inp.setAttribute('min', "0");
			inp.setAttribute('max', "100");
			inp.setAttribute('name', kk);
			inp.setAttribute('value', "0");
			div.appendChild(inp);
			return inp;
		}

		const div0 = _create_group();
		const div1 = _create_group('par-1');
		const div2 = _create_group('par-2');

		const lbl0 = _create_lbl(div0);
		lbl0.innerHTML = key.toUpperCase() + "-Taper";

		const sel0 = document.createElement("select");
		sel0.id = k;
		div0.appendChild(sel0);

		_create_lbl(div1, 'par-1');
		_create_input(div1, 'par-1');
		_create_lbl(div2, 'par-2');
		_create_input(div2, 'par-2');
		return new SceneControlTaper(parent, key);
	}
}

export class SceneControlAllTapers extends SceneControl{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, ['taper-sampling']);
		this.xControl = new SceneControlTaper(parent, 'x', parent.find_element('taper-x-group'));
		this.yControl = new SceneControlTaper(parent, 'y', parent.find_element('taper-y-group'));
		this.add_event_types('taper-changed');
	}
	get calculationWaiting(){
		return (
			this.xControl.calculationWaiting
			|| this.yControl.calculationWaiting
			|| this.changed['taper-sampling']
		);
	}
	control_changed(key){
		super.control_changed(key);
		const eleX = this.parent.find_element('taper-x-group');
		const eleY = this.parent.find_element('taper-y-group');
		if (this.find_element('taper-sampling')[1].selected){
			eleY.style.display = 'none';
			eleX.querySelector("label").innerHTML = "R-Taper";
		}
		else{
			eleY.style.display = 'block';
			eleX.querySelector("label").innerHTML = "X-Taper";
		}
		if (key === 'taper-sampling') this.request_recompute();
	}
	/**
	* Add callable objects to queue.
	*
	* @param {SceneQueue} queue
	*
	* @return {null}
	* */
	add_to_queue(queue){
		this.xControl.add_to_queue(queue);
		this.yControl.add_to_queue(queue);
	}
	/**
	* Add callable objects to queue AFTER phased array
	* is created.
	*
	* @param {SceneQueue} queue
	* @param {SceneControlPhasedArray} src
	*
	* @return {null}
	* */
	add_calculator_queue(queue, src){
		if (this.find_element('taper-sampling')[0].selected){
			let taperX, taperY;
			// we're doing x/y sampling.
			queue.add("Calculating X taper...", () => {
				this.clear_changed('taper-sampling');
				const t = this.xControl.activeTaper;
				const geo = src.pa.geometry;
				taperX = t.calculate_from_geometry(geo.x, geo.dx);
			});
			queue.add("Calculating Y taper...", () => {
				const t = this.yControl.activeTaper;
				const geo = src.pa.geometry;
				taperY = t.calculate_from_geometry(geo.y, geo.dy);
			});
			queue.add("Multiplying tapers...", () => {
				this.trigger_event('taper-changed');
				src.pa.set_magnitude_weight(Float32Array.from(taperX, (x, i) => x * taperY[i]));
			});
		}
		else{
			// we're doing r sampling.
			queue.add("Calculating taper...", () => {
				const t = this.xControl.activeTaper;
				const geo = src.pa.geometry;
				src.pa.set_magnitude_weight(t.calculate_from_radial_geometry(geo));
				this.trigger_event('taper-changed');
			});
		}
	}
	create_samples(points, axis){
		const x = linspace(-1, 1, points);
		const dx = x[1] - x[0];
		let y;
		if (this.find_element('taper-sampling')[0].selected){
			if (axis == 'x') y = this.xControl.activeTaper.calculate_from_geometry(x, dx);
			else y = this.yControl.activeTaper.calculate_from_geometry(x, dx);
		}
		else y = this.xControl.activeTaper.calculate_from_geometry(x, dx);
		return [x, y];
	}
}

export class SceneControlFarfieldDomain extends SceneControlWithSelector{
	static autoUpdateURL = false;
	constructor(parent, key){
		super(parent, key, FarfieldDomains);
		this.ff = null;
		this._ffDirty = false;
		this._wasMatrix = false;
		this.validMaxMonitors = new Set(['directivity', 'pattern-metrics']);
		this.maxMonitors = {};
		this.add_event_types('farfield-changed', 'farfield-calculation-complete', 'matrix-plot-ready');
	}
	isMatrixDomain(){
		const kls = this.selected_class();
		return kls != null && (kls.domain === 'z' || kls.domain === 's');
	}
	matrixKind(){
		if (!this.isMatrixDomain()) return null;
		return this.selected_class().domain;
	}
	control_changed(key){
		super.control_changed(key);
		if (key === 'farfield-domain'){
			this.request_recompute();
			return;
		}
		if (key !== 'farfield-frequency' || (this.ff === null && !this.isMatrixDomain())) return;
		this.parent.update_url_parameters();
	}
	/**
	* Add callable functions to monitor values.
	*
	* @param {string} key Examples: directivity
	* @param {function(Number):null} callback
	*
	* @return {null}
	* */
	add_max_monitor(key, callback){
		if (!(this.validMaxMonitors.has(key))){
			throw Error(`Invalid monitor ${key}. Expected: ${Array.from(this.validMaxMonitors).join(', ')}`)
		}
		if (!(key in this.maxMonitors)) this.maxMonitors[key] = [];
		this.maxMonitors[key].push(callback);
	}
	_notify_farfield_complete(){
		this.trigger_event('farfield-calculation-complete', this.ff);
		for (const [key, value] of Object.entries(this.maxMonitors)){
			let val;
			if (key == 'directivity') val = this.ff.dirMax;
			else if (key == 'pattern-metrics') val = this.ff.patternMetrics;
			else throw Error(`Unknown max key ${key}.`)
			value.forEach((e) => e(val));
		}
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
		if (this.isMatrixDomain()){
			if (this.changed['farfield-frequency'] || arrayControl.farfieldNeedsCalculation){
				this._ffDirty = true;
			}
			this._wasMatrix = true;
			this.clear_changed('farfield-domain', 'farfield-points', 'farfield-uv-bound', 'farfield-frequency');
			queue.add('Notifying matrix plot...', () => {
				this.trigger_event('matrix-plot-ready', arrayControl.pa);
			});
			this.needsRedraw = false;
			return;
		}

		let needsRecalc = arrayControl.farfieldNeedsCalculation;
		const domainClassChanged = this.ff != null && this.ff.domain !== this.selected_class().domain;
		if (this.changed['farfield-points'] || this.changed['farfield-frequency'] || this.changed['farfield-uv-bound'] || this.ff === null || this._ffDirty || domainClassChanged){
			queue.add('Creating farfield mesh...', () => {
				this.ff = this.build_active_object();
				this.trigger_event('farfield-changed', this.ff);
			});
			needsRecalc = true;
		}
		if (needsRecalc){
			queue.add_iterator('Calculating farfield...', () => {
				return this.ff.calculator_loop(arrayControl.pa)
			});
			queue.add("Notifying farfield change...", () => {
				this._ffDirty = false;
				this._wasMatrix = false;
				this._notify_farfield_complete();
			})
		}
		else if (this._wasMatrix && this.ff != null){
			this._wasMatrix = false;
			queue.add("Notifying farfield change...", () => {
				this._notify_farfield_complete();
			});
		}
		this.needsRedraw = needsRecalc;
	}
}

export class SceneControlSteeringDomain extends SceneControlWithSelector{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, 'steering-domain', SteeringDomains);
		this._last = this.selected_class();
		document.getElementById(this.prepend + '-steer-broadside').addEventListener('click', () => {
			this.reset_to_broadside();
		});
	}
	get calculationWaiting(){return this.changed['theta'] || this.changed['phi'] || this.changed['steering-domain']};
	reset_to_broadside(){
		const c1 = this.find_object_map('theta');
		const c2 = this.find_object_map('phi');
		c1.set_value(0);
		c2.set_value(0);
		c1.ele.dispatchEvent(new Event('change'));
		c2.ele.dispatchEvent(new Event('change'));
		let scene = this.parent;
		while (scene != null && typeof scene.update_url_parameters !== 'function'){
			scene = scene.parent;
		}
		if (scene != null) scene.update_url_parameters();
		this.request_recompute();
	}
	control_changed(key){
		if (key == this.primaryKey){
			if (this._last === undefined) return;
			const c1 = this.find_object_map('theta');
			const c2 = this.find_object_map('phi');
			const p1 = Number(c1.ele.value);
			const p2 = Number(c2.ele.value);
			const obj = this.build_active_object();
			let [n1, n2] = obj.from(this._last.title, p1, p2);
			if (isNaN(n1) || isNaN(n2)){
				n1 = 0.0;
				n2 = 0.0;
			}
			c1.set_value(n1);
			c2.set_value(n2);
			this._last = this.selected_class();
		}
		super.control_changed(key);
	}
	get_theta_phi(){
		const obj = this.build_active_object();
		this.clear_changed('theta', 'phi', 'steering-domain');
		return [obj.theta_deg, obj.phi_deg];
	}
}

export class SceneControlPower extends SceneControl{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, ['power', 'power-unit', 'power-scope', 'eirp-unit']);
		this._updating = false;
		this.enteredWatts = defaultWatts();
		this.unit = 'dBm';
		this.scope = 'element';
		this.eirpUnit = 'dBW';
		this._scopeRadios = {
			element: document.getElementById(this.prepend + '-power-scope-element'),
			array: document.getElementById(this.prepend + '-power-scope-array'),
		};
		this.add_event_types('power-changed');
		this.addEventListener('scene-loaded', () => { this.load_from_dom(); });
		this.addEventListener('reset', () => { this.apply_defaults(); });
		for (const radio of Object.values(this._scopeRadios)){
			if (radio === null) continue;
			radio.addEventListener('change', () => {
				if (this._updating || !radio.checked) return;
				this.find_element('power-scope').value = radio.value;
				this.find_element('power-scope').dispatchEvent(new Event('change'));
			});
		}
	}
	control_changed(key){
		if (this._updating) return;
		super.control_changed(key);
		if (key === 'power'){
			const raw = String(this.find_element('power').value).trim();
			if (raw === '' || raw === '-' || raw === '.' || raw === '-.') return;
			const watts = wattsFrom(raw, this.unit);
			if (!Number.isFinite(watts) || watts < 0) return;
			this.enteredWatts = watts;
			this.trigger_event('power-changed');
			return;
		}
		if (key === 'power-unit'){
			const u = this.find_element('power-unit').value;
			if (!isPowerUnit(u)) return;
			this.unit = u;
			this.write_power_field();
			this.trigger_event('power-changed');
			return;
		}
		if (key === 'power-scope'){
			const next = this.find_element('power-scope').value;
			if (!isPowerScope(next) || next === this.scope){
				this.sync_scope_radios();
				return;
			}
			const wsum = this.powerWeightSum();
			if (this.scope === 'element' && next === 'array') this.enteredWatts *= wsum;
			else if (this.scope === 'array' && next === 'element' && wsum > 0) this.enteredWatts /= wsum;
			this.scope = next;
			this.sync_scope_radios();
			this.write_power_field();
			this.trigger_event('power-changed');
			return;
		}
		if (key === 'eirp-unit'){
			const u = this.find_element('eirp-unit').value;
			if (!isPowerUnit(u)) return;
			this.eirpUnit = u;
			this.trigger_event('power-changed');
		}
	}
	powerWeightSum(){
		const pa = this.parent.pa;
		if (pa === null || pa === undefined) return 1;
		const s = pa.powerWeightSum;
		if (!Number.isFinite(s) || s <= 0) return 1;
		return s;
	}
	write_power_field(){
		this._updating = true;
		this.find_element('power').value = formatPowerValue(this.enteredWatts, this.unit);
		this._updating = false;
	}
	write_all_fields(){
		this._updating = true;
		this.find_element('power').value = formatPowerValue(this.enteredWatts, this.unit);
		this.find_element('power-unit').value = this.unit;
		this.find_element('power-scope').value = this.scope;
		this.find_element('eirp-unit').value = this.eirpUnit;
		this.sync_scope_radios();
		this._updating = false;
	}
	sync_scope_radios(){
		const radio = this._scopeRadios[this.scope];
		if (radio !== undefined && radio !== null) radio.checked = true;
	}
	load_from_dom(){
		const unit = this.find_element('power-unit').value;
		this.unit = isPowerUnit(unit) ? unit : 'dBm';
		const eirpUnit = this.find_element('eirp-unit').value;
		this.eirpUnit = isPowerUnit(eirpUnit) ? eirpUnit : 'dBW';
		const scope = this.find_element('power-scope').value;
		this.scope = isPowerScope(scope) ? scope : 'element';
		let watts = wattsFrom(this.find_element('power').value, this.unit);
		if (!Number.isFinite(watts) || watts < 0) watts = defaultWatts();
		this.enteredWatts = watts;
		this.write_all_fields();
	}
	apply_defaults(){
		this.enteredWatts = defaultWatts();
		this.unit = 'dBm';
		this.scope = 'element';
		this.eirpUnit = 'dBW';
		this.write_all_fields();
	}
	/**
	 * Full-scale per-element power in watts.
	 * If the user entered total array power, this is P_tot / Σ|a_i|².
	 */
	getWatts(){
		if (this.scope !== 'array') return this.enteredWatts;
		const wsum = this.powerWeightSum();
		if (wsum <= 0) return 0;
		return this.enteredWatts / wsum;
	}
	getArrayWatts(){
		if (this.scope === 'array') return this.enteredWatts;
		return this.enteredWatts * this.powerWeightSum();
	}
	getUnit(){ return this.unit; }
	getEirpUnit(){ return this.eirpUnit; }
	format(watts, decimals){ return formatPower(watts, this.unit, decimals); }
	formatEirp(watts, decimals){ return formatPower(watts, this.eirpUnit, decimals); }
}

export class SceneControlIllumination extends SceneControlWithSelectorAutoBuild{
	static autoUpdateURL = false;
	constructor(parent){
		super(parent, 'illumination-type', Illuminations, parent.find_element('illumination-controls'));
		this.activeIllumination = null;
	}
	control_changed(key){
		super.control_changed(key);
		this.activeIllumination = null;
	}
	get calculationWaiting(){ return this.activeIllumination === null; }
	add_to_queue(queue){
		if (this.calculationWaiting){
			queue.add('Building illumination...', () => {
					this.activeIllumination = this.build_active_object();
				}
			)
		}
	}
}

const MATCH_HELP = "Set Z0 to Re(Z11) and Self X to −Im(Z11) of one isolated element at the current frequency. Not a scan-impedance match; array Sii still has mutual leftover.";

export class SceneControlElement extends SceneControlWithSelectorAutoBuild{
	static autoUpdateURL = false;
	constructor(parent){
		const host = parent.find_element('element-controls');
		super(parent, 'element-type', ElementTypes, host);
		this.activeElement = null;

		const div = document.createElement('div');
		div.classList = "form-group";
		div.id = parent.prepend + "-element-n-div";
		const lbl = document.createElement('label');
		lbl.setAttribute('for', parent.prepend + "-element-n");
		lbl.innerHTML = "n";
		const ele = document.createElement('input');
		ele.type = 'number';
		ele.id = parent.prepend + "-element-n";
		ele.name = ele.id;
		ele.readOnly = true;
		ele.tabIndex = -1;
		div.appendChild(lbl);
		div.appendChild(ele);
		host.appendChild(div);
		this.nDiv = div;
		this.nInput = ele;

		const note = document.createElement('div');
		note.className = 'form-note';
		note.id = parent.prepend + '-element-green-note';
		note.textContent = 'S,T from the PEC Green function. First build can take seconds near 32×32.';
		note.title = 'Unique-lag PEC kernel fills Z; Kurokawa S,T from that Z at the common real Z0. The 32×32 LU can take seconds.';
		note.style.display = 'none';
		host.appendChild(note);
		this.greenNote = note;

		const matchDiv = document.createElement('div');
		matchDiv.classList = "form-group";
		matchDiv.id = parent.prepend + "-element-match-div";
		matchDiv.style.display = 'none';
		const matchLbl = document.createElement('label');
		matchLbl.setAttribute('for', parent.prepend + "-element-match");
		matchLbl.textContent = "Match";
		matchLbl.title = MATCH_HELP;
		const matchBtn = document.createElement('button');
		matchBtn.type = 'button';
		matchBtn.id = parent.prepend + "-element-match";
		matchBtn.textContent = "Match";
		matchBtn.title = MATCH_HELP;
		matchDiv.title = MATCH_HELP;
		matchDiv.appendChild(matchLbl);
		matchDiv.appendChild(matchBtn);
		host.appendChild(matchDiv);
		this.matchDiv = matchDiv;
		this.matchBtn = matchBtn;
		matchBtn.addEventListener('click', () => this.apply_green_match());

		this.addEventListener('active-class-changed', () => {
			this.release_gain_html_min();
			this.update_n_display();
			this.update_green_note();
			if (typeof this.parent.sync_mode_radios === 'function') this.parent.sync_mode_radios();
		});
		this.install_gain_editing();
		this.update_n_display();
		this.update_green_note();
	}
	install_gain_editing(){
		const gainEle = this.find_element('element-gain');
		const maxG = ElementCosN.controls['element-gain'].max;
		const typing = (ev) => {
			const t = ev.inputType || '';
			return t.startsWith('insertText') || t.startsWith('insertFromPaste') || t.startsWith('delete');
		};
		gainEle.addEventListener('focus', () => {
			gainEle.min = 0;
		});
		gainEle.addEventListener('input', (ev) => {
			if (typing(ev)){
				this.update_n_display();
				return;
			}
			this.clamp_gain_field();
			this.update_n_display();
		});
		gainEle.addEventListener('blur', () => {
			this.clamp_gain_field(true);
			gainEle.min = 0;
			gainEle.max = maxG;
			this.update_n_display();
		});
		this.release_gain_html_min();
		gainEle.max = maxG;
	}
	release_gain_html_min(){
		const gainEle = this.find_element('element-gain', false);
		if (gainEle == null) return;
		gainEle.min = 0;
	}
	clamp_gain_field(fillEmpty){
		const gainEle = this.find_element('element-gain');
		const maxG = ElementCosN.controls['element-gain'].max;
		const raw = String(gainEle.value).trim();
		if (raw === '' || raw === '-' || raw === '.' || raw === '-.'){
			if (!fillEmpty) return;
			gainEle.value = MIN_ELEMENT_GAIN_DBI.toFixed(2);
			return;
		}
		let v = Number(raw);
		if (!Number.isFinite(v)){
			if (!fillEmpty) return;
			gainEle.value = MIN_ELEMENT_GAIN_DBI.toFixed(2);
			return;
		}
		if (v < MIN_ELEMENT_GAIN_DBI) v = MIN_ELEMENT_GAIN_DBI;
		if (v > maxG) v = maxG;
		const shown = v <= MIN_ELEMENT_GAIN_DBI + 1e-9 ? MIN_ELEMENT_GAIN_DBI.toFixed(2) : String(v);
		if (gainEle.value !== shown) gainEle.value = shown;
	}
	control_changed(key){
		super.control_changed(key);
		this.activeElement = null;
		this.update_n_display();
		if (key === this.primaryKey) this.request_recompute();
	}
	update_n_display(){
		if (!this.nDiv) return;
		const isCos = this.selected_class() === ElementCosN;
		this.nDiv.style.display = isCos ? "flex" : "none";
		if (!isCos) return;
		const gainEle = this.find_element('element-gain');
		const raw = gainEle ? String(gainEle.value).trim() : '';
		const gain = Number(raw);
		const g = Number.isFinite(gain) ? Math.max(gain, MIN_ELEMENT_GAIN_DBI) : ElementCosN.controls['element-gain'].default;
		this.nInput.value = exponentFromPeakDbi(g).toFixed(3);
	}
	update_green_note(){
		const isPec = this.selected_class() === ElementGreenPec;
		if (this.greenNote) this.greenNote.style.display = isPec ? 'block' : 'none';
		if (this.matchDiv) this.matchDiv.style.display = isPec ? 'flex' : 'none';
		const typeEle = this.find_element('element-type', false);
		if (typeEle){
			if (isPec && ElementGreenPec.help) typeEle.title = ElementGreenPec.help;
			else typeEle.removeAttribute('title');
		}
	}
	apply_green_match(){
		if (this.selected_class() !== ElementGreenPec) return;
		const ep = this.build_active_object();
		const freq = typeof this.parent.frequencyScale === 'function' ? this.parent.frequencyScale() : 1;
		let z;
		try {
			z = zSelfPecDipole(ep.h, ep.ell, ep.a, freq);
		}
		catch {
			return;
		}
		const re = z[0];
		const im = z[1];
		if (!Number.isFinite(re) || !Number.isFinite(im) || !(re > 0)) return;
		const z0Ctrl = ElementGreenPec.controls['element-z0'];
		let zc = re;
		if (z0Ctrl.min != null && zc < z0Ctrl.min) zc = z0Ctrl.min;
		if (z0Ctrl.max != null && zc > z0Ctrl.max) zc = z0Ctrl.max;
		const z0Ele = this.find_element('element-z0');
		const xEle = this.find_element('element-xself');
		z0Ele.value = String(zc);
		xEle.value = String(-im);
		z0Ele.dispatchEvent(new Event('input', {bubbles: true}));
		xEle.dispatchEvent(new Event('input', {bubbles: true}));
		z0Ele.dispatchEvent(new Event('change', {bubbles: true}));
		xEle.dispatchEvent(new Event('change', {bubbles: true}));
		let scene = this.parent;
		while (scene != null && typeof scene.update_url_parameters !== 'function'){
			scene = scene.parent;
		}
		if (scene != null) scene.update_url_parameters();
	}
	get calculationWaiting(){ return this.activeElement === null; }
	add_to_queue(queue){
		if (this.calculationWaiting){
			queue.add('Building element pattern...', () => {
					this.activeElement = this.build_active_object();
					this.update_n_display();
					this.update_green_note();
				}
			)
		}
	}
}

export class SceneTaperCuts extends ScenePlot1D{
	draw(){
		this.reset();
		this.set_xlabel('Window');
		this.set_ylabel('Magnitude');
		this.set_xgrid(-0.5, 0.5, 11);
		this.set_xgrid_points(1);

		const pa = this.arrayScene;
		if (pa === undefined || pa == null) return;
		const taper = pa.taperControl;
		if (taper === undefined || taper === null) return;

		let belowZero = false;
		this.legend_items().forEach((e) => {
			const v = e.getAttribute('data-axis');
			if (v !== null){
				const [x, y] = taper.create_samples(101, v);
				const maxV = Math.max(...Float32Array.from(y, (i) => Math.abs(i)));
				const minV = Math.min(...y);
				if (minV < 0) belowZero = true;
				if (x !== null) this.add_data(x, Float32Array.from(y, (i) => i/maxV), e);
			}
		});
		if (belowZero) this.set_ygrid(-1, 1, 11);
		else this.set_ygrid(0, 1, 11);
		super.draw();
	}
	/**
	* Bind a Phased Array Scene.
	*
	* @param {SceneControlPhasedArray} scene
	*
	* @return {null}
	* */
	bind_phased_array_scene(scene){
		this.arrayScene = scene;
		scene.taperControl.addEventListener('taper-changed', () => {
			this.draw();
		});
	}
}
