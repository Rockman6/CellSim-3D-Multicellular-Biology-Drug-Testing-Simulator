// First-person movement: WASD, mouse look, jump, sneak and sprint, with
// simple circle-against-box collision so you can't walk through benches.
import * as THREE from "three";
import { settings } from "./settings.js";

const EYE = 1.62, SNEAK_EYE = 1.38, SIT_EYE = 1.22, RADIUS = 0.26;

export class Player {
  constructor(camera) {
    this.camera = camera;
    this.pos = new THREE.Vector3(3, 0, 2.9);
    this.yaw = 2.35;
    this.pitch = -0.12;
    this.vy = 0;
    this.keys = new Set();
    this.eye = EYE;
    this.bob = 0;
    this.sitting = null;
    this.zoom = 1;
    this.moved = false;
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
  }
  update(dt, colliders) {
    const k = this.keys;
    let fx = 0, fz = 0;
    if (k.has("KeyW")) fz += 1;
    if (k.has("KeyS")) fz -= 1;
    if (k.has("KeyA")) fx -= 1;
    if (k.has("KeyD")) fx += 1;
    const moving = fx || fz;
    if (moving && this.sitting) this.sitting = null;   // walking stands you up
    const sneak = k.has("ShiftLeft") || k.has("ShiftRight");
    const sprint = !sneak && (k.has("ControlLeft") || k.has("ControlRight"));
    const speed = sneak ? 1.2 : sprint ? 4.4 : 2.6;
    this.moved = false;
    if (moving) {
      const len = Math.hypot(fx, fz);
      const sin = Math.sin(this.yaw), cos = Math.cos(this.yaw);
      // Forward is -z when yaw = 0.
      const vx = ((fx * cos - fz * sin) / len) * speed * dt;
      const vz = ((-fx * sin - fz * cos) / len) * speed * dt;
      this.tryMove(vx, 0, colliders);
      this.tryMove(0, vz, colliders);
      this.bob += dt * (sprint ? 13 : 9);
      this.moved = true;
    }
    // Gravity and jumping.
    if (k.has("Space") && this.pos.y <= 0.0001 && !this.sitting) this.vy = 3.6;
    if (this.vy !== 0 || this.pos.y > 0) {
      this.vy -= 9.8 * dt;
      this.pos.y = Math.max(0, this.pos.y + this.vy * dt);
      if (this.pos.y === 0) this.vy = 0;
      this.moved = true;
    }
    const targetEye = this.sitting ? SIT_EYE : sneak ? SNEAK_EYE : EYE;
    this.eye += (targetEye - this.eye) * Math.min(1, dt * 10);
    const bobY = settings.bob && moving && this.pos.y === 0 ? Math.sin(this.bob) * 0.025 : 0;
    const cam = this.camera;
    cam.position.set(this.pos.x, this.pos.y + this.eye + bobY, this.pos.z);
    cam.rotation.set(this.pitch, this.yaw, 0, "YXZ");
    const fov = settings.fov / this.zoom;
    if (Math.abs(cam.fov - fov) > 0.01) { cam.fov += (fov - cam.fov) * Math.min(1, dt * 12); cam.updateProjectionMatrix(); this.moved = true; }
  }
  tryMove(dx, dz, colliders) {
    const nx = this.pos.x + dx, nz = this.pos.z + dz;
    for (const c of colliders) {
      if (c.top !== undefined && this.pos.y > c.top - 0.05) continue;   // jumping over low things
      const cx = Math.max(c.x0, Math.min(nx, c.x1)), cz = Math.max(c.z0, Math.min(nz, c.z1));
      if ((nx - cx) ** 2 + (nz - cz) ** 2 < RADIUS * RADIUS) return;
    }
    this.pos.x = nx; this.pos.z = nz;
  }
}
