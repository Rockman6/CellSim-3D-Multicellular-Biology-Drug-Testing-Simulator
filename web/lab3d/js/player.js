// First-person movement: WASD, mouse look, jump, sneak and sprint. You speed
// up and slow down over a fraction of a second instead of starting and
// stopping dead, slide along benches instead of sticking to them, and the
// head bob follows your stride, with a footstep at each step.
import * as THREE from "three";
import { settings } from "./settings.js";

const EYE = 1.62, SNEAK_EYE = 1.36, SIT_EYE = 1.22, RADIUS = 0.26;
const SPEED = { sneak: 1.0, walk: 2.1, sprint: 3.9 };   // m/s: a brisk indoor walk
const ACCEL = 11, BRAKE = 14, STRIDE = 0.72;            // m/s², m/s², metres per step

export class Player {
  constructor(camera) {
    this.camera = camera;
    this.pos = new THREE.Vector3(3, 0, 2.9);
    this.vel = new THREE.Vector2();
    this.yaw = 2.35;
    this.pitch = -0.12;
    this.vy = 0;
    this.keys = new Set();
    this.eye = EYE;
    this.phase = 0;      // stride phase, in steps
    this.sway = 0;
    this.dip = 0;        // the knees bending when you land
    this.sitting = null;
    this.zoom = 1;
    this.moved = false;
    this.onStep = null;  // (loudness 0..1) for footstep sounds
  }
  look(dx, dy) {
    const k = 0.0022 * settings.sensitivity;
    this.yaw -= dx * k;
    this.pitch = THREE.MathUtils.clamp(this.pitch - dy * k * (settings.invertY ? -1 : 1), -1.5, 1.5);
  }
  sit(stool) {
    const p = stool.position;
    this.sitting = stool;
    this.pos.set(p.x, 0, p.z);
    this.vel.set(0, 0);
  }
  update(dt, colliders) {
    const k = this.keys;
    let fx = 0, fz = 0;
    if (k.has("KeyW")) fz += 1;
    if (k.has("KeyS")) fz -= 1;
    if (k.has("KeyA")) fx -= 1;
    if (k.has("KeyD")) fx += 1;
    const input = fx || fz;
    if (input && this.sitting) this.sitting = null;   // walking stands you up
    const sneak = k.has("ShiftLeft") || k.has("ShiftRight");
    const sprint = !sneak && fz > 0 && (k.has("ControlLeft") || k.has("ControlRight"));
    const top = sneak ? SPEED.sneak : sprint ? SPEED.sprint : SPEED.walk;
    // The velocity you want, in world space (forward is -z when yaw = 0).
    const want = new THREE.Vector2();
    if (input && !this.sitting) {
      const len = Math.hypot(fx, fz), sin = Math.sin(this.yaw), cos = Math.cos(this.yaw);
      want.set(((fx * cos - fz * sin) / len) * top, ((-fx * sin - fz * cos) / len) * top);
    }
    const ground = this.pos.y <= 0.0001;
    const rate = (want.lengthSq() > 0 ? ACCEL : BRAKE) * (ground ? 1 : 0.25);
    const dv = want.clone().sub(this.vel), step = rate * dt;
    if (dv.length() <= step) this.vel.copy(want); else this.vel.addScaledVector(dv.normalize(), step);
    this.moved = false;
    const speed = this.vel.length();
    if (speed > 0.001) {
      const before = this.pos.clone();
      this.pos.x += this.vel.x * dt;
      this.pos.z += this.vel.y * dt;
      this.collide(colliders);
      // Lose the speed you pushed into a wall, so you slide along it.
      const moved = this.pos.clone().sub(before);
      if (dt > 0) this.vel.set(moved.x / dt, moved.z / dt).clampLength(0, top + 0.01);
      if (ground) this.stride(Math.hypot(moved.x, moved.z), speed / SPEED.walk);
      this.moved = true;
    } else this.vel.set(0, 0);
    // Gravity and jumping.
    if (k.has("Space") && ground && !this.sitting) this.vy = 3.6;
    if (this.vy !== 0 || this.pos.y > 0) {
      const wasUp = this.pos.y > 0;
      this.vy -= 9.8 * dt;
      this.pos.y = Math.max(0, this.pos.y + this.vy * dt);
      if (this.pos.y === 0) {
        if (wasUp) { this.dip = Math.min(0.09, -this.vy * 0.022); this.onStep?.(0.9); }
        this.vy = 0;
      }
      this.moved = true;
    }
    // Eye height: standing, sneaking or sitting, plus the landing dip.
    const targetEye = this.sitting ? SIT_EYE : sneak ? SNEAK_EYE : EYE;
    if (Math.abs(targetEye - this.eye) > 0.001) this.moved = true;
    this.eye += (targetEye - this.eye) * Math.min(1, dt * 9);
    if (this.dip > 0.0005) { this.dip *= Math.exp(-dt * 9); this.moved = true; } else this.dip = 0;
    // Head bob: down at each footfall, a little side to side, scaled by speed.
    const amp = settings.bob && ground ? Math.min(1.3, speed / SPEED.walk) : 0;
    const s = Math.sin(this.phase * Math.PI);
    const bobY = amp * (0.016 * (s * s) - 0.008);
    this.sway += ((amp ? Math.sin(this.phase * Math.PI * 0.5 + 0.4) * 0.011 * amp : 0) - this.sway) * Math.min(1, dt * 10);
    const cam = this.camera;
    const sin = Math.sin(this.yaw), cos = Math.cos(this.yaw);
    cam.position.set(this.pos.x + cos * this.sway, this.pos.y + this.eye + bobY - this.dip, this.pos.z - sin * this.sway);
    cam.rotation.set(this.pitch, this.yaw, amp ? -this.sway * 0.25 : 0, "YXZ");
    const fov = settings.fov / this.zoom;
    if (Math.abs(cam.fov - fov) > 0.01) { cam.fov += (fov - cam.fov) * Math.min(1, dt * 12); cam.updateProjectionMatrix(); this.moved = true; }
  }
  // Advance the stride; a footstep sounds each time a foot comes down.
  stride(dist, loud) {
    const before = Math.floor(this.phase);
    this.phase += dist / STRIDE;
    if (Math.floor(this.phase) !== before) this.onStep?.(Math.min(1, 0.35 + 0.5 * loud));
  }
  // Push the body (a circle) out of every box it overlaps. Doing it twice
  // settles corners, where two boxes push at once.
  collide(colliders) {
    for (let pass = 0; pass < 2; pass++) {
      for (const c of colliders) {
        if (c.top !== undefined && this.pos.y > c.top - 0.05) continue;   // jumping over low things
        const cx = Math.max(c.x0, Math.min(this.pos.x, c.x1)), cz = Math.max(c.z0, Math.min(this.pos.z, c.z1));
        let dx = this.pos.x - cx, dz = this.pos.z - cz;
        const d2 = dx * dx + dz * dz;
        if (d2 >= RADIUS * RADIUS) continue;
        if (d2 < 1e-10) {
          // The centre is inside the box: leave by the nearest side.
          const opts = [[c.x0 - RADIUS - this.pos.x, 0], [c.x1 + RADIUS - this.pos.x, 0], [0, c.z0 - RADIUS - this.pos.z], [0, c.z1 + RADIUS - this.pos.z]];
          opts.sort((a, b) => Math.abs(a[0] + a[1]) - Math.abs(b[0] + b[1]));
          this.pos.x += opts[0][0]; this.pos.z += opts[0][1];
          continue;
        }
        const d = Math.sqrt(d2), push = RADIUS - d;
        this.pos.x += (dx / d) * push;
        this.pos.z += (dz / d) * push;
      }
    }
  }
}
