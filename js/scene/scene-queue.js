/**
 * Queue iterator object.
 *
 * @typedef {Object} QueueIteratorResult
 * @property {String} text - Text to display on progress bar.
 * @property {number} progress - Progress value.
 * @property {number} max - Maximum value expected on `progress`.
 */

export class SceneQueue{
	constructor(progressElement, statusElement){
		this.progress = progressElement;
		this.status = statusElement;
		this.generation = 0;
		this.reset();
		this.channel = new MessageChannel();
		this.channel.port1.onmessage = (ev) => {
			if (ev.data !== this.generation) return;
			this.process_queue();
		};
		this.running = false;
		this._pendingBuild = null;
	}
	/**
	 * Run `buildFn` now, or keep only the latest callback if a run is already
	 * in progress. `buildFn` should fill the queue and call `start()`; `reset()`
	 * is applied here so an in-flight run is never cancelled.
	 *
	 * @param {function():void} buildFn
	 */
	request(buildFn){
		this._pendingBuild = buildFn;
		if (!this.running) this._drainPending();
	}
	_drainPending(){
		if (this.running) return;
		const pending = this._pendingBuild;
		if (pending == null) return;
		this._pendingBuild = null;
		this.reset();
		pending();
	}
	/**
	* Add callable object to queue.
	*
	* @param {String} text String to display on status bar
	* @param {function():null} func Callback function.
	*
	* @return {null}
	* */
	add(text, func){
		this.queue.push({
			text: text,
			func: func,
			type: 'function',
		});
	}
	/**
	* Add iterator to queue. Progress bar will be updated to match
	* progress of iterator. Sync and async generators are both accepted.
	*
	* @param {String} text String to display on status bar
	* @param {function():Iterator<QueueIteratorResult>|AsyncIterator<QueueIteratorResult>} func Callback iterator that yields information.
	*
	* @return {null}
	* */
	add_iterator(text, func){
		this.queue.push({
			text: text,
			func: func,
			type: 'iterator',
		});
	}
	dump(){
		this.queue.forEach((entry) => {
			console.log("---- " + entry['text']);
			console.log(entry['func']);
		});
	}
	next(){
		if (this.queue.length == 0) return null;
		return this.queue.shift();
	}
	reset(){
		this.generation++;
		this.queue = [];
		this.startingLength = 0;
		this._current = null;
		this.running = false;
	}
	start(finalText){
		if (finalText === undefined) finalText = "Complete";
		this.finalText = finalText;
		const prog = this.progress;
		this.startingLength = this.queue.length;
		this._current = null;
		prog.value = 0;
		prog.max = this.startingLength;
		if (!this.running) this.process_queue();
	}
	get length(){ return this.queue.length; }

	log(string, toConsole){
		if (string !== undefined){
			if (toConsole !== false) console.log(string)
			this.status.innerHTML = string;
		}
	}
	_continue(gen){
		if (gen !== this.generation) return;
		this.channel.port2.postMessage(gen);
	}
	process_queue(){
		const gen = this.generation;
		const prog = this.progress;
		const cont = () => { this._continue(gen); };
		this.running = true;
		if (this._current === null){
			this._current = this.next();
			prog.value = prog.max - this.length;
			if (this._current === null){
				this.log(this.finalText, false);
				this.running = false;
				this._drainPending();
				return;
			}
			this.log(this._current['text']);
			cont();
			return;
		}
		const c = this._current;
		if (c['type'] == 'next'){
			Promise.resolve(c['func'].next()).then((v) => {
				if (gen !== this.generation) return;
				if (v.done) {
					this._current = null;
					prog.max = this.startingLength;
					prog.value = prog.max - this.length;
				}
				else{
					this.log(v.value['text']);
					prog.max = v.value['max'];
					prog.value = v.value['progress'];
				}
				cont();
			}, (err) => {
				console.error(err);
				if (gen !== this.generation) return;
				this._current = null;
				this.log(String(err), true);
				cont();
			});
			return;
		}
		if (c['type'] == 'iterator'){
			this._current['type'] = 'next';
			this._current['func'] = c['func']();
			cont();
			return;
		}
		try {
			c['func']();
			this._current = null;
			cont();
		} catch (err) {
			console.error(err);
			if (gen !== this.generation) return;
			this._current = null;
			this.log(String(err), true);
			cont();
		}
	}
}
