# CellSim — Computational Cell Biology Simulator & In Silico Drug Testing Platform

A high-performance, real-time 3D cell biology simulator with integrated pharmacology for in silico drug testing. Built with C++20 and Apple Metal for GPU-accelerated rendering.

![Platform](https://img.shields.io/badge/platform-macOS-blue)
![Language](https://img.shields.io/badge/language-C%2B%2B20-orange)
![GPU](https://img.shields.io/badge/GPU-Apple%20Metal-silver)
![License](https://img.shields.io/badge/license-MIT-green)

## Overview

CellSim models individual cells as autonomous agents with full intracellular biology — CDK/Cyclin cell cycle ODE, dual-pathway metabolism, Warburg effect, nutrient diffusion, contact inhibition, telomere erosion, and apoptosis. Users can apply chemotherapy drugs (Cisplatin, Doxorubicin, Paclitaxel, 5-FU) and observe dose-dependent cell kill, phase arrest, and resistance emergence in real time.

## Features

### Cell Biology Engine
- **CDK/Cyclin ODE** — Novak-Tyson 7-variable model (CycD, Rb, E2F, CycE, CycA, CycB, p21) driving G1→S→G2→M phase transitions
- **Fick Diffusion Field** — 64×64 grid for O₂, glucose, CO₂, pH, and drug concentration with Dirichlet boundary conditions
- **Dual-Pathway Metabolism** — glycolysis + oxidative phosphorylation with Warburg switch under hypoxia
- **Contact Inhibition** — Hippo-YAP pathway: mechanical pressure → p21/p27 induction → reversible G0 arrest at confluence
- **Telomere Erosion** — 20bp loss per division → replicative senescence (Hayflick limit)
- **Multi-Timescale Architecture** — FAST (Ca²⁺), MEDIUM (ATP/ROS/stress), SLOW (cell cycle/genome)

### Drug Testing Platform
- **4 Pre-built Drugs**: Cisplatin, Doxorubicin, Paclitaxel, 5-Fluorouracil
- **PhysiPKPD Pharmacodynamics** — Hill equation dose-response (EC50, Hill coefficient)
- **4 Mechanisms of Action**: anti-proliferative, pro-apoptotic, DNA damage, mitochondrial toxin
- **Drug Diffusion** — Fick diffusion of drug through tissue with cellular uptake/efflux
- **Resistance Mutations** — MDR pump upregulation via stochastic mutation at division
- **Application Modes** — uniform bath or point injection with wash-out

### Visualization & Analysis
- **3D Metal Rendering** — translucent bio-luminescent cells with GLB organelle models (nucleus, ER, Golgi, mitochondria)
- **ImGui Research Interface** — population stats, CDK/Cyclin ODE readout, metabolism panel, drug treatment controls
- **ImPlot Time-Series** — live plots of population dynamics, ATP, stress, phase distribution
- **CSV Data Export** — population time-series + per-cell snapshots with native Save dialog

## Screenshots

*(Run the simulator to see: translucent cells with glowing organelles, ImGui panels with real-time ODE data, drug dose-response curves)*

## Build

### Prerequisites (macOS)
```bash
# Xcode Command Line Tools
xcode-select --install

# Homebrew + build tools
brew install cmake

# vcpkg package manager
git clone https://github.com/microsoft/vcpkg.git ~/vcpkg
~/vcpkg/bootstrap-vcpkg.sh

# Dependencies
~/vcpkg/vcpkg install glfw3 glm "imgui[glfw-binding,metal-binding]" implot nlohmann-json stb cgltf
```

### Build & Run
```bash
cd CellSim
cmake -B build -DCMAKE_TOOLCHAIN_FILE=$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/CellSim
```

### Organelle Models
Place GLB files in `assets/models/`:
- `nucleus.glb`, `smooth ER.glb`, `rough ER.glb`, `golgi apparatus.glb`, `mitochondria.glb`

## Usage

### Controls
| Input | Action |
|-------|--------|
| Drag | Rotate camera |
| Right-drag | Pan |
| Scroll | Zoom |
| ESC | Quit |

### Drug Testing Workflow
1. Let colony grow to confluence (~500 cells) at 20× speed
2. Open **Drug Treatment** panel
3. Select drug (e.g., Cisplatin) and set concentration (e.g., 2 µM = EC50)
4. Click **Apply Uniform**
5. Observe viability drop, phase distribution changes, cell death
6. Click **WASH OUT** — survivors resume growth
7. **Export CSV** for quantitative analysis

## Scientific Basis

### Cell Cycle Model
Novak-Tyson minimal CDK/Cyclin oscillator with checkpoint gates:
- **G1/S checkpoint**: CycE/CDK2 > threshold, Rb phosphorylated, p21 < 0.5
- **G2/M checkpoint**: CycB/CDK1 (MPF) > threshold, DNA damage cleared
- **Contact inhibition**: Hippo-YAP → p27 induction at confluence

### Pharmacodynamics (adapted from PhysiPKPD)
```
dC_internal/dt = uptake × C_external − efflux × C_internal
dDamage/dt = Hill(C_int, EC50, n) × maxEffect − repairRate × Damage
Hill(C, EC50, n) = C^n / (EC50^n + C^n)
```

### Drug Parameters

| Drug | EC50 (µM) | Hill | Mechanism | Reference |
|------|-----------|------|-----------|-----------|
| Cisplatin | 2.0 | 1.5 | DNA damage | PMC2751448 |
| Doxorubicin | 0.5 | 2.0 | DNA damage + mito toxin | PMC1501422 |
| Paclitaxel | 0.01 | 2.5 | M-phase arrest | PMC2751448 |
| 5-Fluorouracil | 5.0 | 1.2 | S-phase block + DNA damage | BioModels |

## References

1. Novak B & Tyson JJ (2008). "Design principles of biochemical oscillators." *Nat Rev Mol Cell Biol* 9:981-991.
2. Ghaffarizadeh A et al. (2018). "PhysiCell: An open source physics-based cell simulator for 3-D multicellular systems." *PLoS Comput Biol* 14:e1005991.
3. Bergman D et al. (2023). "PhysiPKPD: A pharmacokinetics and pharmacodynamics module for PhysiCell." *GigaByte*.
4. Green DR & Kroemer G (2004). "The pathophysiology of mitochondrial cell death." *Science* 305:626-629.
5. Delarue M et al. (2018). "Compressive stress inhibits proliferation in tumor spheroids through a volume limitation." *Dev Cell*.
6. Blackburn EH (2001). "Switching and signaling at the telomere." *Cell* 106:661-673.
7. Casciari JJ et al. (1992). "Variations in tumor cell growth rates and metabolism with oxygen concentration, glucose concentration, and extracellular pH." *Biotechnol Bioeng*.
8. Vander Heiden MG et al. (2009). "Understanding the Warburg effect." *Science* 324:1029-1033.

## License

MIT License. See [LICENSE](LICENSE) for details.

Portions adapted from PhysiPKPD (BSD-3), PhysiCell (BSD-3), and BioModels Database (CC0).
