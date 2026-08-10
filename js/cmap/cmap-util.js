export class ColormapControl{
	constructor(selector, defaultSelection){
		this.useCaller = (typeof defaultSelection === 'function');
		this.caller = defaultSelection;
		this.changed = true;
		this.selector = selector;
		if (defaultSelection === undefined) defaultSelection = 'viridis';
		if (this.useCaller) defaultSelection = "default";
		this.defaultSelection = defaultSelection;

		const _add_ele = cm => {
			const ele = document.createElement('option');
			ele.value = cm;
			ele.innerHTML = cm;
			selector.appendChild(ele);
			if (defaultSelection == cm) ele.selected = true;
		}
		if (this.useCaller) _add_ele(defaultSelection);
		this.constructor.Colormaps.forEach(_add_ele);
		selector.addEventListener('change', () => {
			this.changed = true;
		});
		window.installThemeChanged(() => {
			this.changed = true;
		});
	}
	addEventListener(e, callback){ this.selector.addEventListener(e, callback); }
	cmap(){
		const cms = this.constructor.Colormaps;
		const find_colormap = this.constructor.find_colormap;
		let offset = 0;

		if (this.useCaller){
			let cs = this.caller();
			let ss = this.selector[0].selected;
			for (let i = 0; i < cms.length; i++){
				if (ss && cms[i] == cs) this.selector[i + 1].style.fontStyle = 'italic';
				else this.selector[i + 1].style.fontStyle = '';
			}
			if (ss) return find_colormap(cs);
			offset += 1;
		}
		for (let i = 0; i < cms.length; i++)
			if (this.selector[i + offset].selected)
				return find_colormap(cms[i]);
		if (this.useCaller) return find_colormap(this.caller());
		return find_colormap(this.defaultSelection);
	}
}

/**
 * Convert HSV to RGB.
 * @param {Number} h - Hue (0-1).
 * @param {Number} s - Saturation (0-1).
 * @param {Number} v - Value (0-1).
 *
 * @returns {String} - rbg color.
 */
export function hsv2rgb(h, s, v){
	let r, g, b;

	const i = Math.floor(h * 6);
	const f = h * 6 - i;
	const p = v * (1 - s);
	const q = v * (1 - f * s);
	const t = v * (1 - (1 - f) * s);

	switch (i % 6){
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		case 5:
			r = v;
			g = p;
			b = q;
			break;
	}

	return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
}
