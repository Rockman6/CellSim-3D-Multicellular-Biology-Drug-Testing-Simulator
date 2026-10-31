# CellSim

A biology lab you walk into in first person. Prepare slides of onion skin, cheek cells, a pondweed leaf, a squashed root tip or red onion skin, then focus a microscope on them. Or grow living A549 lung cancer cells, treat them with a drug and incubate them. The living cells are simulated by the CellSim engine (`cellsim.cell.dish`), running in the browser.

## Run it

Browsers block 3D models and the engine on pages opened straight from a file, so serve the `web` folder:

```bash
cd CellSim/web
python3 -m http.server 8000
```

Then open <http://localhost:8000/lab3d/>.

The engine uses the wheel in `../lab/wheel/` (the same one the Lab page uses). GitHub Pages builds it on every deploy. To build it locally:

```bash
cd CellSim
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
| F | Find anything: type its name, typos and all, and follow the blue marker; the thing is outlined, even through walls. F again stops |
| R | With the micropipette in your hand: set its volume (2–200 µL) |
| Space, Shift, Ctrl | Jump, sneak, sprint |
| Esc | Pause |
| F3 | Frame rate, draw calls and memory |

## The experiments

- **Onion skin cells:** goggles, wash your hands, slide on the white tile, cut the onion with the scalpel, peel the skin with the forceps, one drop of iodine, cover slip, stage, then focus at 4×, 10× and 40×.
- **Your cheek cells:** swab your cheek, smear it on the slide, put the swab in the disinfectant, then add methylene blue.
- **Pondweed leaf cells:** a drop of water first, then one leaf from the tank. At 400× the chloroplasts stream around each cell.
- **Mitosis in a root tip:** cut the last 5 mm of a root from the onion by the windows, soften it in 1 M hydrochloric acid in the 60 °C water bath, stain it with toluidine blue, squash it under the cover slip with a paper towel, then click cells at 400× and name their stage. The journal works out the mitotic index.
- **Osmosis in red onion cells:** mount red onion skin in water, then draw salt solution under the cover slip with a paper towel: the protoplasts shrink from their walls. Draw water under and they swell back.
- **Living cancer cells:** take Dish A from the incubator and look at it, add drops of paclitaxel, doxorubicin or cisplatin, incubate, and compare with Dish B, the untreated control. The doses are real dilutions: two 50 µL drops of 10 µM paclitaxel in 5.1 mL of medium give 0.196 µM. The engine's cell counts appear in the journal.
- **Dose–response curve:** with the micropipette (a fresh tip for every well), put a different volume of one drug into each well of the 6-well plate, keep one well as the control, incubate 72 h and read the plate. The journal plots viability against dose and reads off the IC50. Each well is its own engine run.
- **Counting cells:** lift the cells off a dish with trypsin, mix 10 µL with 10 µL of trypan blue, load the hemocytometer and click every cell in the four corner squares. The count is checked against the engine's dish, scaled from the simulated patch to the whole 60 mm dish.

## Graphics and battery

Settings → Graphics: **Auto** picks a preset from the graphics chip.

| Preset | Detail |
|---|---|
| Low | Older Macs and laptops on battery: 0.75× resolution, no shadows, no smoothing, one ceiling lamp, 30 fps |
| Medium | Shadows from the windows, two lamps |
| High | Sharper shadows, four lamps |
| Ultra | Full resolution, 4096 px shadows, all six lamps, no frame cap |

The game stops drawing when nothing moves, the sun's shadows are redrawn only when something in the lab moves, and the CellSim engine (about 150 MB of memory) only loads once a culture dish needs it. Small parts that share a material are merged into one mesh, so the whole room is about 600 draw calls.

## What's where

| Path | |
|---|---|
| `index.html`, `css/game.css` | Screens and HUD |
| `js/main.js` | Game loop, input, inventory, the microscope stage, incubator, plate reader and drawers |
| `js/player.js`, `js/audio.js` | Walking (speeding up, slowing down, sliding along benches, head bob by stride) and the sound effects, made in code |
| `js/finder.js`, `js/outline.js` | Find anything (typo-tolerant search) and the outlines round what you aim at or look for |
| `js/room.js` | The room, furniture, lights and colliders |
| `js/items.js` | Every item and what left- and right-click do with it |
| `js/props.js` | Equipment built in code: glassware with real wall thickness, bottles, the micropipette, the 6-well plate, the hemocytometer, the water bath, tools |
| `js/liquids.js` | Volumes, colour mixing and drug dilution |
| `js/scope.js`, `js/specimens.js`, `js/culture.js` | The view through the eyepiece |
| `js/engine.js`, `js/engine-worker.js`, `js/lab_engine.py` | The CellSim engine, in a background worker |
| `js/tasks.js`, `js/journal.js` | Experiments and the lab journal |
| `assets/` | Models, textures and lighting (`downloads/` holds the raw Sketchfab files and is not committed) |

## Honest notes

- The living cells are the engine's: positions, cell-cycle phase, p53, caspase-3 and DNA damage come from `run_dish_stream`. Only their outlines are drawn here.
- The onion, cheek, pondweed, root tip and red onion slides are drawn at measured cell sizes and, for the root tip, measured stage proportions (about 12 % of meristem cells dividing). They are not simulated.
- Osmosis follows a simple osmotic model (each cell's sap about 0.55 osmol/L, varying cell to cell; a protoplast shrinks to sap ÷ outside), sped up so it happens in seconds. It is not the CellSim engine.
- Softening the root tip takes 20 s here instead of 5 minutes.

## Credits

- **Cell engine:** CellSim (MIT), in the browser via Pyodide.
- **3D engine:** three.js r160 (MIT).
- **Poly Haven (CC0):** yellow onion, metal stool, steel shelves, fluorescent lights, drawer cabinet, clipboard, stationery, wall clock, magnifying glass, binder notebook, metal office desk, floor tiles, painted plaster, and the hospital room HDRI.
- **Sketchfab (CC Attribution 4.0):**
  - [Chemistry Lab Table](https://sketchfab.com/3d-models/chemistry-lab-table-fc5951d8af674da4bf55aef0693ff1d0) by Jawahar Yokesh
  - [lab bench](https://sketchfab.com/3d-models/lab-bench-1ac2a62c52a848bbaf746146dc7253f8) by Konstantin_Keller (the incubator)
  - [Compound Microscope](https://sketchfab.com/3d-models/compound-microscope-8709be06c5b24549b19bba7a1988a56f) by kg007
  - [Door Classroom](https://sketchfab.com/3d-models/door-classroom-9-mb-fcc87f41b47f49a8941d2b6847246909) by Mehdi Shahsavan
- **Built in code:** the room, slides, cover slips, bottles, beakers, culture dishes, the 6-well plate, the micropipette and tips, microtubes, the hemocytometer, the plate reader, the water bath, the pondweed tank, the lab tools, the safety goggles, the hand sanitizer, the window blinds and the sound effects.
