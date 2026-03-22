#pragma once
#include <simd/simd.h>

// ══════════════════════════════════════════════════════════════════════════
//  Orbit Camera
// ══════════════════════════════════════════════════════════════════════════

class Camera {
public:
    Camera();

    // Input handlers
    void onMouseButton(int button, int action, double x, double y);
    void onMouseMove(double x, double y);
    void onScroll(double yoffset);
    void onResize(int width, int height);

    // Getters
    simd_float4x4 getViewProjection() const;
    simd_float3   getPosition() const;
    float         getAspect() const { return aspect_; }

private:
    void updatePosition();

    // Orbit parameters
    float azimuth_   = 0.3f;
    float elevation_  = 0.6f;
    float distance_   = 90.0f;
    float panX_       = 0.0f;
    float panY_       = -3.0f; // FLOOR_Y + 4

    // Camera state
    simd_float3 position_;
    simd_float3 target_;
    float fov_    = 55.0f;  // degrees
    float aspect_ = 16.0f / 9.0f;
    float nearZ_  = 0.1f;
    float farZ_   = 500.0f;

    // Mouse state
    bool  leftDrag_  = false;
    bool  rightDrag_ = false;
    double lastX_    = 0.0;
    double lastY_    = 0.0;
};
