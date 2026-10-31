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
    title: "Living cancer cells", tag: "Simulated by Cell Sim", blurb: "Treat a dish of A549 lung cancer cells with a drug, incubate it, and compare it with an untreated control.",
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
