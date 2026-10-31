# CellSim Lab

A biology lab you walk into in first person. Prepare slides of onion skin, cheek cells or a pondweed leaf, then focus a microscope on them. Or grow living A549 lung cancer cells, treat them with a drug and incubate them. The living cells are simulated by the Cell Sim engine (`cellsim.cell.dish`), running in the browser.

## Run it

Browsers block 3D models and the engine on pages opened straight from a file, so serve the `web` folder:

```bash
cd "Cell Sim/web"
python3 -m http.server 8000
```

Then open <http://localhost:8000/lab3d/>.

The engine uses the wheel in `../lab/wheel/` (the same one the Lab page uses). GitHub Pages builds it on every deploy. To build it locally:

```bash
cd "Cell Sim"
python -m pip wheel --no-deps --no-build-isolation -w web/lab/wheel .
cd web/lab/wheel && ls -t cellsim-*-py3-none-any.whl | head -1 > LATEST
```

The first time you look at a culture dish, the browser downloads Python (Pyodide, about 30 s). After that it is cached.

## Controls

| Key | Does |
|---|---|
| W A S D, mouse | Walk, look |
| Left click | Take the thing you look at into your pocket |
| Right click | Use or open the thing you look at; with an item in hand, use it on what you look at, or put it down |
| Hold right click | Lower a cover slip slowly (a quick tap drops it flat and traps air) |
| 1–9, scroll | Pick a hotbar slot |
| Q | Put down the item in your hand |
| E | Inventory |
| J | Lab journal |
| Space, Shift, Ctrl | Jump, sneak, sprint |
| Esc | Pause |
| F3 | Frame rate, draw calls and memory |

## The experiments

- **Onion skin cells:** goggles, wash your hands, slide on the white tile, cut the onion with the scalpel, peel the skin with the forceps, one drop of iodine, cover slip, stage, then focus at 4×, 10× and 40×.
- **Your cheek cells:** swab your cheek, smear it on the slide, put the swab in the disinfectant, then add methylene blue.
- **Pondweed leaf cells:** a drop of water first, then one leaf from the tank. At 400× the chloroplasts stream around each cell.
- **Living cancer cells:** take Dish A from the incubator and look at it, add drops of paclitaxel, doxorubicin or cisplatin, incubate, and compare with Dish B, the untreated control. The doses are real dilutions: two 50 µL drops of 10 µM paclitaxel in 5.1 mL of medium give 0.196 µM. The engine's cell counts appear in the journal.

## Graphics and battery

Settings → Graphics: **Auto** picks a preset from the graphics chip.

| Preset | Detail |
|---|---|
| Low | Older Macs and laptops on battery: 0.75× resolution, no shadows, no smoothing, one ceiling lamp, 30 fps |
| Medium | Shadows from the windows, two lamps |
| High | Sharper shadows, four lamps |
| Ultra | Full resolution, 4096 px shadows, all six lamps, no frame cap |

The game stops drawing when nothing moves, and the Cell Sim engine (about 150 MB of memory) only loads once a culture dish needs it.

## What's where

| Path | |
|---|---|
| `index.html`, `css/game.css` | Screens and HUD |
| `js/main.js` | Game loop, input, inventory, the microscope stage, incubator and drawers |
| `js/room.js` | The room, furniture, lights and colliders |
| `js/items.js` | Every item and what left- and right-click do with it |
| `js/props.js` | Equipment built in code: slides, bottles, beakers, dishes, the pondweed tank, tools |
| `js/liquids.js` | Volumes, colour mixing and drug dilution |
| `js/scope.js`, `js/specimens.js`, `js/culture.js` | The view through the eyepiece |
| `js/engine.js`, `js/engine-worker.js`, `js/lab_engine.py` | The Cell Sim engine, in a background worker |
| `js/tasks.js`, `js/journal.js` | Experiments and the lab journal |
| `assets/` | Models, textures and lighting (`downloads/` holds the raw Sketchfab files and is not committed) |

## Honest notes

- The living cells are the engine's: positions, cell-cycle phase, p53, caspase-3 and DNA damage come from `run_dish_stream`. Only their outlines are drawn here.
- The onion, cheek and pondweed slides are drawn at measured cell sizes. They are not simulated.

## Credits

- **Cell engine:** Cell Sim (MIT), in the browser via Pyodide.
- **3D engine:** three.js r160 (MIT).
- **Poly Haven (CC0):** yellow onion, metal stool, steel shelves, fluorescent lights, drawer cabinet, clipboard, stationery, wall clock, magnifying glass, binder notebook, metal office desk, floor tiles, painted plaster, and the hospital room HDRI.
- **Sketchfab (CC Attribution 4.0):**
  - [Chemistry Lab Table](https://sketchfab.com/3d-models/chemistry-lab-table-fc5951d8af674da4bf55aef0693ff1d0) by Jawahar Yokesh
  - [lab bench](https://sketchfab.com/3d-models/lab-bench-1ac2a62c52a848bbaf746146dc7253f8) by Konstantin_Keller (the incubator)
  - [Compound Microscope](https://sketchfab.com/3d-models/compound-microscope-8709be06c5b24549b19bba7a1988a56f) by kg007
  - [Door Classroom](https://sketchfab.com/3d-models/door-classroom-9-mb-fcc87f41b47f49a8941d2b6847246909) by Mehdi Shahsavan
- **Built in code:** the room, slides, cover slips, bottles, beakers, culture dishes, the pondweed tank, the lab tools, the safety goggles, the hand sanitizer and the window blinds.
