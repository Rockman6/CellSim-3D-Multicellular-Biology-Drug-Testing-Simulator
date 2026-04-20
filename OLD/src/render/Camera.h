#pragma once
#include <simd/simd.h>

// ══════════════════════════════════════════════════════════════════════════
//  Camera — Orbit mode + FPS (WASD) mode
// ══════════════════════════════════════════════════════════════════════════

class Camera {
public:
    Camera();

    // Input handlers
    void onMouseButton(int button, int action, double x, double y);
    void onMouseMove(double x, double y);
    void onScroll(double yoffset);
    void onResize(int width, int height);

    // FPS-style WASD movement (call every frame with dt)
    void updateFPS(float dt, bool w, bool a, bool s, bool d,
                   bool space, bool arrowDown);

    // Presets
    void setSingleCellView();
    void setColonyView();

    // Follow a specific cell — smoothly track it.
    // `zoomFactor` ≥ 0 multiplies the framed distance; 1.0 = default framing,
    // 0.5 = zoomed in 2×, 2.0 = zoomed out 2×.
    void followCell(simd_float3 cellPos, float cellRadius, float zoomFactor = 1.0f);

    // Sensitivity (adjustable from UI)
    float rotateSensitivity = 1.0f;  // Multiplier for orbit rotation
    float panSensitivity    = 1.0f;  // Multiplier for right-drag pan
    float zoomSensitivity   = 1.0f;  // Multiplier for scroll zoom
    float moveSensitivity   = 1.0f;  // Multiplier for WASD movement

    // Getters
    simd_float4x4 getViewProjection() const;
    simd_float3   getPosition() const;
    simd_float3   getForward() const;
    float         getAspect() const { return aspect_; }
    float         getDistance() const { return distance_; }

private:
    void updatePosition();

    // Orbit parameters
    float azimuth_   = 0.3f;
    float elevation_  = 0.6f;
    float distance_   = 90.0f;
    float panX_       = 0.0f;
    float panY_       = -3.0f;
    float panZ_       = 0.0f;

    // Camera state
    simd_float3 position_;
    simd_float3 target_;
    float fov_    = 55.0f;
    float aspect_ = 16.0f / 9.0f;
    float nearZ_  = 0.01f;
    float farZ_   = 500.0f;

    // Mouse state
    bool  leftDrag_  = false;
    bool  rightDrag_ = false;
    double lastX_    = 0.0;
    double lastY_    = 0.0;
};
