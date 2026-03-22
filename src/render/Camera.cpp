#include "Camera.h"
#include <cmath>

// ══════════════════════════════════════════════════════════════════════════
//  Orbit Camera implementation
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

Camera::Camera() {
    updatePosition();
}

void Camera::updatePosition() {
    float cosElev = cosf(elevation_);
    float sinElev = sinf(elevation_);
    float cosAzim = cosf(azimuth_);
    float sinAzim = sinf(azimuth_);

    position_ = {
        distance_ * cosElev * sinAzim + panX_,
        distance_ * sinElev + panY_,
        distance_ * cosElev * cosAzim
    };
    target_ = {panX_, panY_ - 2.0f, 0.0f};
}

void Camera::onMouseButton(int button, int action, double x, double y) {
    if (button == 0) leftDrag_  = (action == 1);
    if (button == 1) rightDrag_ = (action == 1);
    lastX_ = x;
    lastY_ = y;
}

void Camera::onMouseMove(double x, double y) {
    double dx = x - lastX_;
    double dy = y - lastY_;
    lastX_ = x;
    lastY_ = y;

    if (leftDrag_) {
        azimuth_   -= (float)dx * 0.005f;
        elevation_  = fmaxf(-0.1f, fminf(1.3f, elevation_ - (float)dy * 0.005f));
        updatePosition();
    }
    if (rightDrag_) {
        panX_ -= (float)dx * 0.02f;
        panY_ += (float)dy * 0.02f;
        updatePosition();
    }
}

void Camera::onScroll(double yoffset) {
    distance_ = fmaxf(8.0f, fminf(300.0f, distance_ + (float)yoffset * 2.0f));
    updatePosition();
}

void Camera::onResize(int width, int height) {
    if (height > 0) aspect_ = (float)width / (float)height;
}

simd_float4x4 Camera::getViewProjection() const {
    simd_float4x4 proj = makePerspective(fov_, aspect_, nearZ_, farZ_);
    simd_float3 up = {0, 1, 0};
    simd_float4x4 view = makeLookAt(position_, target_, up);
    return simd_mul(proj, view);
}

simd_float3 Camera::getPosition() const {
    return position_;
}
