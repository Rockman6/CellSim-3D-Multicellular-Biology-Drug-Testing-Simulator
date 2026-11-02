// The CellSim engine, loaded only when a culture dish needs it, so the lab
// stays light (about 150 MB less memory) until then.
export class Engine {
  constructor(onStatus) {
    this.onStatus = onStatus;
    this.worker = null;
    this.ready = null;
    this.version = "";
    this.pending = new Map();
    this.seq = 0;
  }
  start() {
    if (this.ready) return this.ready;
    this.onStatus("busy", "Starting the CellSim engine…");
    this.worker = new Worker(new URL("./engine-worker.js", import.meta.url));
    this.ready = new Promise((resolve, reject) => {
      this.worker.onmessage = (ev) => {
        const m = ev.data;
        if (m.type === "status") this.onStatus("busy", m.text);
        else if (m.type === "ready") { this.version = m.version; this.onStatus("ready", `CellSim ${m.version} ready`); resolve(); }
        else if (m.type === "fatal") { this.onStatus("error", `Engine failed: ${m.message}`); reject(new Error(m.message)); }
        else if (m.type === "result" || m.type === "error") {
          const p = this.pending.get(m.id);
          if (!p) return;
          this.pending.delete(m.id);
          if (m.type === "result") p.resolve(JSON.parse(m.data)); else p.reject(new Error(m.message));
        }
      };
      this.worker.onerror = (e) => { this.onStatus("error", "The engine could not start. It needs internet the first time."); reject(e); };
    });
    return this.ready;
  }
  async dish(cfg) {
    await this.start();
    const id = ++this.seq;
    this.onStatus("busy", `Simulating ${cfg.label || "the dish"} to ${cfg.t} h…`);
    const out = await new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ id, cmd: "dish", cfg });
    });
    this.onStatus("ready", `CellSim ${this.version} ready`);
    return out;
  }
}
