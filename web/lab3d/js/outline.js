// Outlines around whole objects: white for what you aim at, blue for what you
// searched for. The objects are drawn flat into a small mask with no depth
// test, so a searched-for object shows through walls; one full-screen pass
// then draws the mask's edges over the finished frame. Works for any shape,
// and costs nothing while nothing is outlined.
import * as THREE from "three";

const HOVER = 2, TRACK = 3;   // render layers used only by the mask pass

const EDGE = /* glsl */ `
  uniform sampler2D mask;
  uniform vec2 texel;
  uniform vec3 hoverColor;
  uniform vec3 trackColor;
  uniform float fill;
  varying vec2 vUv;
  void main() {
    vec2 c = texture2D(mask, vUv).rg;
    vec2 m = vec2(0.0);
    for (int i = 0; i < 8; i++) {
      float a = float(i) * 0.785398;
      vec2 d = vec2(cos(a), sin(a)) * texel;
      m = max(m, texture2D(mask, vUv + d * 1.5).rg);
      m = max(m, texture2D(mask, vUv + d * 3.0).rg * 0.7);
    }
    vec2 edge = clamp(m - c, 0.0, 1.0);
    vec4 col = vec4(trackColor, 1.0) * c.g * fill;
    col = mix(col, vec4(trackColor, 1.0), edge.g);
    col = mix(col, vec4(hoverColor, 0.85), edge.r);
    gl_FragColor = col;
  }`;

export class Outliner {
  constructor(renderer, { scale = 0.5 } = {}) {
    this.r = renderer;
    this.scale = scale;
    this.rt = new THREE.WebGLRenderTarget(2, 2, { depthBuffer: false });
    this.mats = {
      [HOVER]: new THREE.MeshBasicMaterial({ color: 0xff0000, depthTest: false, depthWrite: false }),
      [TRACK]: new THREE.MeshBasicMaterial({ color: 0x00ff00, depthTest: false, depthWrite: false, transparent: true, blending: THREE.AdditiveBlending }),
    };
    this.quadMat = new THREE.ShaderMaterial({
      uniforms: {
        mask: { value: this.rt.texture }, texel: { value: new THREE.Vector2(0.5, 0.5) },
        hoverColor: { value: new THREE.Color(0xffffff) }, trackColor: { value: new THREE.Color(0x5b8cff) }, fill: { value: 0.22 },
      },
      vertexShader: "varying vec2 vUv; void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }",
      fragmentShader: EDGE,
      transparent: true, depthTest: false, depthWrite: false,
    });
    this.quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), this.quadMat);
    this.quad.frustumCulled = false;
    this.ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.obj = { [HOVER]: null, [TRACK]: null };
    this.cc = new THREE.Color();
  }
  get active() { return !!(this.obj[HOVER] || this.obj[TRACK]); }
  // Outline `obj` (or nothing) as the aimed-at or the searched-for object.
  set(kind, obj) {
    const L = kind === "hover" ? HOVER : TRACK;
    if (this.obj[L] === obj) return false;
    if (this.obj[L]) this.obj[L].traverse((m) => m.layers.disable(L));
    this.obj[L] = obj || null;
    if (obj) obj.traverse((m) => { if (m.isMesh && m.name !== "hitbox" && m.name !== "shadow" && m.material && m.material.visible !== false) m.layers.enable(L); });
    return true;
  }
  resize(w, h, pixelRatio) {
    const s = this.scale * pixelRatio;
    this.rt.setSize(Math.max(2, Math.round(w * s)), Math.max(2, Math.round(h * s)));
    this.quadMat.uniforms.texel.value.set(1 / this.rt.width, 1 / this.rt.height);
  }
  // Call right after the frame is drawn.
  draw(scene, camera) {
    if (!this.active) return;
    const r = this.r;
    const bg = scene.background, layers = camera.layers.mask, autoClear = r.autoClear, autoShadow = r.shadowMap.autoUpdate;
    r.getClearColor(this.cc);
    const clearAlpha = r.getClearAlpha();
    scene.background = null;
    r.shadowMap.autoUpdate = false;
    r.autoClear = false;
    r.setRenderTarget(this.rt);
    r.setClearColor(0x000000, 0);
    r.clear(true, false, false);
    for (const L of [HOVER, TRACK]) {
      if (!this.obj[L]) continue;
      camera.layers.set(L);
      scene.overrideMaterial = this.mats[L];
      r.render(scene, camera);
    }
    scene.overrideMaterial = null;
    camera.layers.mask = layers;
    scene.background = bg;
    r.setRenderTarget(null);
    r.setClearColor(this.cc, clearAlpha);
    r.render(this.quad, this.ortho);
    r.autoClear = autoClear;
    r.shadowMap.autoUpdate = autoShadow;
  }
}
