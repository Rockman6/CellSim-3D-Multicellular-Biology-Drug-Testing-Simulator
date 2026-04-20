cmake_minimum_required(VERSION 3.24)
project(CellSimOLD VERSION 1.0 LANGUAGES CXX OBJCXX)

# ══════════════════════════════════════════════════════════════════════════
#  CellSim / OLD — frozen Campaign-2 prototype
#
#  This subtree is the pre-restart cell simulator (HeLa agent model +
#  p53/MDM2/apoptosis + cisplatin/Nutlin benchmarks) preserved as a
#  regression snapshot. The active Campaign-1 tree lives in the repo's
#  top-level src/ + CMakeLists.txt.
#
#  Build standalone from this directory:
#      cmake -B build && cmake --build build && ctest --test-dir build
#  Or via the parent project's add_subdirectory(OLD).
# ══════════════════════════════════════════════════════════════════════════

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OBJCXX_STANDARD 20)

# ── vcpkg integration ────────────────────────────────────────────────────
if(DEFINED ENV{VCPKG_ROOT})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
elseif(EXISTS "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
endif()

set(OLD_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(OLD_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")

# Assets (Metal shaders) live at the repo root, one level up from OLD/.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../assets/shaders")
    set(OLD_ASSETS "${CMAKE_CURRENT_SOURCE_DIR}/../assets")
else()
    set(OLD_ASSETS "${CMAKE_CURRENT_SOURCE_DIR}/assets")
endif()

# ── Find packages (GUI target only — headless tests do not need these) ───
find_package(glfw3 CONFIG QUIET)
find_package(glm CONFIG QUIET)
find_package(imgui CONFIG QUIET)
find_package(implot CONFIG QUIET)

find_library(METAL_FRAMEWORK Metal)
find_library(METALKIT_FRAMEWORK MetalKit)
find_library(COCOA_FRAMEWORK Cocoa)
find_library(QUARTZCORE_FRAMEWORK QuartzCore)
find_library(IOKIT_FRAMEWORK IOKit)
find_library(UTI_FRAMEWORK UniformTypeIdentifiers)

find_path(CGLTF_INCLUDE_DIRS "cgltf.h")

# ── GUI target (optional — only if Metal + deps are available) ───────────
if(glfw3_FOUND AND glm_FOUND AND imgui_FOUND AND implot_FOUND AND METAL_FRAMEWORK)
    set(SOURCES
        ${OLD_SRC}/main.mm
        ${OLD_SRC}/gpu/MetalContext.mm
        ${OLD_SRC}/render/Camera.cpp
        ${OLD_SRC}/render/MeshLibrary.mm
        ${OLD_SRC}/render/GLBLoader.mm
        ${OLD_SRC}/simulation/biochem/Bioagent.cpp
        ${OLD_SRC}/simulation/biochem/TargetLibrary.cpp
    )

    add_executable(CellSimOLD ${SOURCES})

    target_include_directories(CellSimOLD PRIVATE
        ${OLD_SRC}
        ${CGLTF_INCLUDE_DIRS}
    )

    target_link_libraries(CellSimOLD PRIVATE
        glfw
        glm::glm
        imgui::imgui
        implot::implot
        ${METAL_FRAMEWORK}
        ${METALKIT_FRAMEWORK}
        ${COCOA_FRAMEWORK}
        ${QUARTZCORE_FRAMEWORK}
        ${IOKIT_FRAMEWORK}
        ${UTI_FRAMEWORK}
    )

    target_compile_options(CellSimOLD PRIVATE
        -Wall -Wextra -Wno-unused-parameter
        $<$<CONFIG:Release>:-O2>
        $<$<CONFIG:Debug>:-O0 -g>
    )

    set_source_files_properties(
        ${OLD_SRC}/main.mm
        ${OLD_SRC}/gpu/MetalContext.mm
        ${OLD_SRC}/render/MeshLibrary.mm
        ${OLD_SRC}/render/GLBLoader.mm
        PROPERTIES COMPILE_FLAGS "-x objective-c++"
    )

    # ── Compile Metal shaders ────────────────────────────────────────────
    set(METAL_SHADER_DIR "${OLD_ASSETS}/shaders")
    set(METAL_OUTPUT_DIR "${CMAKE_BINARY_DIR}")
    execute_process(
        COMMAND xcrun --find metal
        RESULT_VARIABLE METAL_COMPILER_FOUND
        OUTPUT_QUIET ERROR_QUIET
    )
    if(METAL_COMPILER_FOUND EQUAL 0 AND EXISTS ${METAL_SHADER_DIR})
        file(GLOB METAL_SOURCES "${METAL_SHADER_DIR}/*.metal")
        set(METAL_AIR_FILES "")
        foreach(METAL_SRC ${METAL_SOURCES})
            get_filename_component(METAL_NAME ${METAL_SRC} NAME_WE)
            set(AIR_FILE "${METAL_OUTPUT_DIR}/${METAL_NAME}.air")
            add_custom_command(
                OUTPUT ${AIR_FILE}
                COMMAND xcrun -sdk macosx metal -c ${METAL_SRC} -o ${AIR_FILE}
                DEPENDS ${METAL_SRC}
                COMMENT "Compiling Metal shader: ${METAL_NAME}.metal"
            )
            list(APPEND METAL_AIR_FILES ${AIR_FILE})
        endforeach()
        add_custom_command(
            OUTPUT "${METAL_OUTPUT_DIR}/default.metallib"
            COMMAND xcrun -sdk macosx metallib ${METAL_AIR_FILES} -o "${METAL_OUTPUT_DIR}/default.metallib"
            DEPENDS ${METAL_AIR_FILES}
            COMMENT "Linking Metal shader library"
        )
        add_custom_target(MetalShadersOLD ALL
            DEPENDS "${METAL_OUTPUT_DIR}/default.metallib"
        )
        add_dependencies(CellSimOLD MetalShadersOLD)
    endif()
else()
    message(STATUS "OLD: GUI deps missing — skipping CellSimOLD binary; headless tests still build.")
endif()

# ── Headless validators (the 6 regression witnesses) ─────────────────────
function(cellsim_old_test name src_file)
    add_executable(${name}
        ${OLD_TESTS}/${src_file}
        ${OLD_SRC}/simulation/biochem/Bioagent.cpp
        ${OLD_SRC}/simulation/biochem/TargetLibrary.cpp)
    target_include_directories(${name} PRIVATE ${OLD_SRC})
    target_compile_options(${name} PRIVATE
        -Wall -Wextra -Wno-unused-parameter -O2)
endfunction()

cellsim_old_test(cellsim_headless               headless_doubling.cpp)
cellsim_old_test(cellsim_infection              headless_infection.cpp)
cellsim_old_test(cellsim_p53_pulsing            headless_p53_pulsing.cpp)
cellsim_old_test(cellsim_cisplatin_arrest       headless_cisplatin_arrest.cpp)
cellsim_old_test(cellsim_cisplatin_apoptosis    headless_cisplatin_apoptosis.cpp)
cellsim_old_test(cellsim_blind_benchmarks       headless_blind_benchmarks.cpp)

message(STATUS "OLD: 6 headless validators configured (Campaign-2 regression snapshot).")
