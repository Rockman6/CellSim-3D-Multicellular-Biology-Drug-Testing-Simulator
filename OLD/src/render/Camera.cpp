#include "Camera.h"
#include <cmath>

// ══════════════════════════════════════════════════════════════════════════
//  Camera — Hybrid orbit + FPS free-look
//  Based on LearnOpenGL camera tutorial + arcball best practices
//  Left-drag: orbit around target (azimuth + elevation)
//  Right-drag: pan (shift target in screen space)
//  Scroll: zoom (distance)
//  WASD: FPS-style movement in facing direction
//  Mouse sensitivity scales for comfortable control at all zoom levels
// ══════════════════════════════════════════════════════════════════════════

static simd_float4x4 makePerspective(float fovDeg, float aspect, float nearZ, float farZ) {
    float fovRad = fovDeg * (M_PI / 180.0f);
    float ys = 1.0f / tanf(fovRad * 0.5f);
    float xs = ys / aspect;
    float zs = farZ / (nearZ - farZ);
    return (simd_float4x4){{
        {xs,  0,   0,         0},
        {0,   ys,  0,         0},
        {0,   0,   zs,       -1},
        {0,   0,   zs*nearZ,  0}
    }};
}

static simd_float4x4 makeLookAt(simd_float3 eye, simd_float3 center, simd_float3 up) {
    simd_float3 f = simd_normalize(center - eye);
    simd_float3 s = simd_normalize(simd_cross(f, up));
    simd_float3 u = simd_cross(s, f);
    return (simd_float4x4){{
        { s.x,  u.x, -f.x, 0},
        { s.y,  u.y, -f.y, 0},
        { s.z,  u.z, -f.z, 0},
        {-simd_dot(s, eye), -simd_dot(u, eye), simd_dot(f, eye), 1}
    }};
}

Camera::Camera() { updatePosition(); }

void Camera::updatePosition() {
    float cosElev = cosf(elevation_);
    float sinElev = sinf(elevation_);
    float cosAzim = cosf(azimuth_);
    float sinAzim = sinf(azimuth_);

    // Camera position on sphere around target
    position_ = {
        panX_ + distance_ * cosElev * sinAzim,
        panY_ + distance_ * sinElev,
        panZ_ + distance_ * cosElev * cosAzim
    };
    target_ = {panX_, panY_, panZ_};
}

void Camera::onMouseButton(int button, int action, double x, double y) {
    if (button == 0) leftDrag_  = (action == 1);
    if (button == 1) rightDrag_ = (action == 1);
    lastX_ = x; lastY_ = y;
}

void Camera::onMouseMove(double x, double y) {
    double dx = x - lastX_;
    double dy = y - lastY_;
    lastX_ = x; lastY_ = y;

    // ── Orbit rotation (left drag) ──────────────────────────────────
    // Sensitivity: slower when zoomed in for fine control
    if (leftDrag_) {
        float sens = 0.003f * rotateSensitivity;
        azimuth_   -= (float)dx * sens;
        elevation_  = fmaxf(-1.55f, fminf(1.55f, elevation_ - (float)dy * sens));
        updatePosition();
    }

    if (rightDrag_) {
        float panS = distance_ * 0.0008f * panSensitivity;
        float rightX =  cosf(azimuth_);
        float rightZ = -sinf(azimuth_);
        panX_ -= (float)dx * panS * rightX;
        panZ_ -= (float)dx * panS * rightZ;
        panY_ += (float)dy * panS;
        updatePosition();
    }
}

void Camera::onScroll(double yoffset) {
    // Proportional zoom — feels the same at all distances
    float speed = distance_ * 0.10f * zoomSensitivity;
    distance_ = fmaxf(0.3f, fminf(300.0f, distance_ + (float)yoffset * speed));
    updatePosition();
}

void Camera::onResize(int width, int height) {
    if (height > 0) aspect_ = (float)width / (float)height;
}

// ── FPS-style WASD movement ─────────────────────────────────────────────
// Moves the orbit center (target) in the camera's facing direction.
// W/S = forward/back along camera look direction (horizontal component)
// A/D = strafe left/right
// Space = up, Shift = down
void Camera::updateFPS(float dt, bool w, bool a, bool s, bool d,
                       bool space, bool arrowDown) {
    if (!w && !a && !s && !d && !space && !arrowDown) return;

    float speed = fmaxf(distance_ * 0.35f, 0.15f) * dt * moveSensitivity;

    // Forward: from camera toward target, projected to horizontal XZ
    float fwdX = -sinf(azimuth_);
    float fwdZ = -cosf(azimuth_);
    // Right: perpendicular to forward (90° clockwise in XZ)
    float rgtX = -fwdZ;
    float rgtZ =  fwdX;

    float mx = 0, my = 0, mz = 0;
    if (w) { mx += fwdX * speed; mz += fwdZ * speed; }
    if (s) { mx -= fwdX * speed; mz -= fwdZ * speed; }
    if (a) { mx -= rgtX * speed; mz -= rgtZ * speed; }
    if (d) { mx += rgtX * speed; mz += rgtZ * speed; }
    if (space)     my += speed * 0.6f;
    if (arrowDown) my -= speed * 0.6f;

    panX_ += mx;
    panY_ += my;
    panZ_ += mz;
    updatePosition();
}

simd_float4x4 Camera::getViewProjection() const {
    simd_float4x4 proj = makePerspective(fov_, aspect_, nearZ_, farZ_);
    simd_float3 up = {0, 1, 0};
    simd_float4x4 view = makeLookAt(position_, target_, up);
    return simd_mul(proj, view);
}

simd_float3 Camera::getPosition() const { return position_; }

simd_float3 Camera::getForward() const {
    return simd_normalize(target_ - position_);
}

void Camera::setSingleCellView() {
    azimuth_ = 0.3f; elevation_ = 0.5f; distance_ = 12.0f;
    panX_ = 0.0f; panY_ = -5.0f; panZ_ = 0.0f;
    updatePosition();
}

void Camera::setColonyView() {
    azimuth_ = 0.3f; elevation_ = 0.6f; distance_ = 90.0f;
    panX_ = 0.0f; panY_ = -3.0f; panZ_ = 0.0f;
    updatePosition();
}

void Camera::followCell(simd_float3 cellPos, float cellRadius, float zoomFactor) {
    // Smoothly lerp pan toward cell position
    float lerpRate = 0.05f;
    panX_ += (cellPos.x - panX_) * lerpRate;
    panY_ += (cellPos.y - panY_) * lerpRate;
    panZ_ += (cellPos.z - panZ_) * lerpRate;
    // Zoom to show the cell nicely (5× radius = default framing), scaled by zoomFactor.
    if (zoomFactor < 0.05f) zoomFactor = 0.05f;
    float targetDist = cellRadius * 5.0f * zoomFactor;
    distance_ += (targetDist - distance_) * lerpRate;
    updatePosition();
}
