// The experiments you can choose, as steps that tick off when the game sees
// you do them. Free roam has no steps.
export const EXPERIMENTS = {
  onion: {
    title: "Onion skin cells", tag: "Plant cells", blurb: "Peel one layer of cells, stain it with iodine and find the nuclei.",
    steps: [
      ["goggles", "Put on the safety goggles (on the microscope bench)"],
      ["hands", "Wash your hands at the sink"],
      ["slide-on-tile", "Take a clean slide and put it on the white tile"],
      ["onion-cut", "Hold the scalpel and cut a piece of onion"],
      ["sample:onion", "With the forceps, peel a thin layer of skin and lay it on the slide"],
      ["drop:iodine", "Add one drop of iodine"],
      ["cover:angle", "Lower a cover slip slowly (hold right-click)"],
      ["stage", "Clip the slide onto the microscope stage"],
      ["focus:4", "Look through the eyepiece and focus with the 4× lens"],
      ["focus:10", "Switch to the 10× lens and sharpen the image"],
      ["focus:40", "Switch to the 40× lens. Fine focus only"],
      ["dispose:slide", "Clean up: the slide goes in the glass waste bin"],
    ],
  },
  cheek: {
    title: "Your cheek cells", tag: "Animal cells", blurb: "Swab your own cheek, stain the cells blue and compare them with plant cells.",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands at the sink"],
      ["slide-on-tile", "Take a clean slide and put it on the white tile"],
      ["swabbed", "Take a cotton swab and rub the inside of your cheek"],
      ["sample:cheek", "Smear the swab across the middle of the slide"],
      ["dispose:swab", "Drop the used swab in the disinfectant beaker"],
      ["drop:mb", "Add one drop of methylene blue"],
      ["cover:angle", "Lower a cover slip slowly (hold right-click)"],
      ["stage", "Clip the slide onto the microscope stage"],
      ["focus:4", "Focus with the 4× lens"],
      ["focus:10", "Switch to 10×"],
      ["focus:40", "Switch to 40×. Find a dark-blue nucleus"],
      ["dispose:slide", "Clean up: slide into the glass waste"],
    ],
  },
  elodea: {
    title: "Pondweed leaf cells", tag: "Plant cells, alive", blurb: "Mount one living leaf and watch its chloroplasts stream around each cell.",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands at the sink"],
      ["slide-on-tile", "Take a clean slide and put it on the white tile"],
      ["drop:water", "Put one drop of distilled water on the slide"],
      ["sample:elodea", "With the forceps, take one leaf from the tank and lay it in the drop"],
      ["cover:angle", "Lower a cover slip slowly (hold right-click)"],
      ["stage", "Clip the slide onto the microscope stage"],
      ["focus:4", "Focus with the 4× lens"],
      ["focus:10", "Switch to 10×"],
      ["focus:40", "Switch to 40× and watch the chloroplasts move"],
      ["dispose:slide", "Clean up: slide into the glass waste"],
    ],
  },
  culture: {
    title: "Living cancer cells", tag: "Simulated by CellSim", blurb: "Treat a dish of A549 lung cancer cells with a drug, incubate it, and compare it with an untreated control.",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands"],
      ["dish-taken", "Open the incubator and take out Dish A"],
      ["dish-stage", "Clip Dish A onto the microscope stage"],
      ["observe:culture", "Look at the living cells with the 10× lens"],
      ["drug", "Take the dish off the stage and add 2 drops of one drug from the vials"],
      ["incubate", "Put Dish A back in the incubator and incubate it for 24 h"],
      ["observe:after", "Look at Dish A again. What changed?"],
      ["compare", "Incubate Dish B (no drug) for the same time and look at it too"],
    ],
  },
  dose: {
    title: "Dose–response curve", tag: "Simulated by CellSim", blurb: "Treat the six wells of a plate with rising doses of a drug, incubate for 72 h and find the dose that kills half the cells (the IC50).",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands"],
      ["plate-taken", "Open the incubator and take out Plate 1"],
      ["tip", "Hold the micropipette and fit a tip from the tip box"],
      ["dose:3", "Dose 3 or more wells, each with a different volume of one drug (R sets the volume; a fresh tip for every well)"],
      ["control", "Leave at least one well untreated: it is the control"],
      ["incubate:72", "Put the plate back in the incubator and let 72 h pass"],
      ["read", "Slide the plate into the plate reader on the microscope bench"],
      ["dr-curve", "Open the journal (J), Cell culture tab: find the IC50 on the curve"],
    ],
  },
  count: {
    title: "Counting cells", tag: "Simulated by CellSim", blurb: "Lift living cells off a dish, stain them with trypan blue and count them in a hemocytometer: how many cells, and how many alive?",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands"],
      ["dish-taken", "Open the incubator and take out a dish"],
      ["trypsin", "Hold the trypsin bottle and right-click the dish to lift the cells off"],
      ["tip", "Fit a tip on the micropipette and set it to 10 µL (press R)"],
      ["trypan-mixed", "Draw up 10 µL of cells and add them to the trypan blue tube in the rack"],
      ["hemo-loaded", "With a fresh tip, draw up 10 µL of the mix and load the hemocytometer"],
      ["hemo-stage", "Clip the hemocytometer onto the microscope stage"],
      ["focus:10", "Focus at 10× on the grid"],
      ["count:4", "Click every cell in the four corner squares to count it"],
      ["count:done", "Work out the cells per mL and how many are alive"],
    ],
  },
  mitosis: {
    title: "Mitosis in a root tip", tag: "Plant cells dividing", blurb: "Soften an onion root tip in warm acid, stain and squash it, then find cells in every stage of mitosis.",
    steps: [
      ["goggles", "Put on the safety goggles: you'll use acid"],
      ["hands", "Wash your hands"],
      ["root-cut", "By the windows: cut the last 5 mm off a root with the scalpel, and pick it up with the forceps"],
      ["acid", "Add 3 drops of 1 M hydrochloric acid to the watch glass in the water bath"],
      ["root-acid", "Put the root tip in the warm acid with the forceps"],
      ["root-soft", "When it's soft (5 minutes, 20 s here), take it out with the forceps"],
      ["sample:root", "Lay it on a clean slide"],
      ["drop:toluidine", "Add one drop of toluidine blue"],
      ["cover:angle", "Lower a cover slip slowly (hold right-click)"],
      ["squash", "Hold the paper towels and press down hard on the cover slip (hold right-click)"],
      ["stage", "Clip the slide onto the microscope stage"],
      ["focus:40", "Find the small cells just behind the tip at 400×"],
      ["classify:25", "Click 25 cells and say which stage of mitosis each is in"],
      ["stages:all", "Find one cell in each stage: prophase, metaphase, anaphase and telophase"],
    ],
  },
  osmosis: {
    title: "Osmosis in red onion cells", tag: "Plant cells, water", blurb: "Watch red onion cells shrink away from their walls in salt water, then swell back in pure water.",
    steps: [
      ["goggles", "Put on the safety goggles"],
      ["hands", "Wash your hands"],
      ["slide-on-tile", "Take a clean slide and put it on the white tile"],
      ["redonion-cut", "Cut the red onion with the scalpel"],
      ["sample:redonion", "Peel the thin red skin with the forceps and lay it on the slide"],
      ["drop:water", "Add one drop of distilled water"],
      ["cover:angle", "Lower a cover slip slowly (hold right-click)"],
      ["stage", "Clip it onto the stage and focus at 100×: red sap fills each cell"],
      ["irrigate:salt", "Off the stage, add a drop of salt solution at one edge, then draw it under with a paper towel"],
      ["osmo:plasmolysed", "Look again: the cells shrink away from their walls"],
      ["irrigate:water", "Draw distilled water under the same way"],
      ["osmo:recovered", "Watch the cells swell back against their walls"],
    ],
  },
};

export class Tasks {
  constructor(onChange) {
    this.exp = null;
    this.done = new Set();
    this.counts = {};
    this.onChange = onChange;
  }
  start(id) {
    this.exp = id ? EXPERIMENTS[id] : null;
    this.id = id;
    this.done = new Set();
    this.counts = {};
    this.onChange();
  }
  // Events arrive from the game; drops are counted for "add 2 drops".
  event(name) {
    if (!this.exp) return;
    this.counts[name] = (this.counts[name] || 0) + 1;
    const map = { "focus-any": null };
    let key = map[name] === undefined ? name : map[name];
    if (name === "drug-drop" && this.counts[name] >= 2) key = "drug";
    if (!key || this.done.has(key)) return;
    if (!this.exp.steps.some(([k]) => k === key)) return;
    this.done.add(key);
    this.onChange(key);
  }
  view() {
    if (!this.exp) return null;
    const steps = this.exp.steps.map(([k, text]) => ({ k, text, done: this.done.has(k) }));
    const now = steps.findIndex((s) => !s.done);
    return { title: this.exp.title, steps, now, finished: now === -1 };
  }
}
