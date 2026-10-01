**Yes — you can build a strong, non-AI cell simulator that meaningfully approaches real cell behavior.**
It will never reach "exact perfect" replacement for novel drugs or unknown viruses (the scale gap and open-ended biology make that physically impossible), but you can build the best physics-first, interpretable, uncertainty-aware hybrid simulator possible. This is exactly what leading non-AI groups (Covert lab, Karr lab, and others) have done successfully for minimal cells and specific pathways.

This is the **best non-AI approach** — a clean, hierarchical, physics-grounded roadmap tailored to a "no black-box" discipline and the existing CellSim foundation. Deterministic / mechanistic methods, rule-based integration, rigorous validation. No neural nets, no end-to-end learned surrogates.

> Amendment: 2026-04-20. This non-AI roadmap supersedes the earlier 5-campaign roadmap which included an ML-potential layer (archived at `docs/goal_v1_2026-04-20.txt`). Layer 1.5 MACE-OFF23 has been removed. Campaign 1 has 7 layers. All UQ is Sobol / Monte-Carlo / parameter sweeps + MAPIE conformal as post-hoc statistical wrapper. Docking uses AutoDock Vina as primary.

### Core Philosophy (non-AI only)
- **Physics-first at every scale**: quantum mechanics → classical MD → coarse-grained models → agent-based rules → PDEs. Never replace a module with a black-box predictor.
- **Modular hybrid integration**: each biological process gets its own mathematical representation (ODEs for metabolism / signaling, stochastic simulation for gene expression, agent-based for cytoskeleton / traffic, PDE for diffusion). Link them through shared state variables and time-step synchronisation (Covert / Karr style).
- **Strict interpretability and mass conservation**: every rate constant must cite literature or physics-derived priors. Enforce mass-balance invariants at every integration step.
- **Uncertainty via ensembles and sensitivity**: parameter sweeps, Monte-Carlo sampling, and Sobol analysis (SALib) instead of neural uncertainty.

### Recommended Non-AI Roadmap (5 Campaigns)

**Campaign 1: Atomic → Molecular Foundation (Years 1–2)**
Accurate small-system physics without AI shortcuts. Seven layers:

- **1.1 Chem foundation**: RDKit for parsing + 3D embedding + OpenFF/AMBER force fields for parametrisation (AM1-BCC charges where possible, MMFF fallback). Cache all parameters.
- **1.2 Classical MD**: OpenMM with ff14SB / ff19SB + TIP3P (or better AMBER/CHARMM variants). Langevin integrator with fixed seeds for reproducibility. Support vacuum, solvated boxes, and short trajectories (10–100 ns on GPU).
- **1.3 Docking + free energy**: AutoDock Vina for poses (empirical, auditable scoring). Alchemical FEP (perses) for ΔG on small sets. Cache poses and energies. GNINA CNN-scored mode may ship only as an explicitly labeled "fast-guess" option alongside Vina, never as the sole value.
- **1.4 Quantum basics**: xTB or PySCF for reactive fragments, transition states, small-molecule electronics (HOMO / LUMO, partial charges). Use for parametrising unusual chemistries.
- **1.5 Coarse-grained**: Martini 3 for membranes, lipids, and large complexes. Run 1–10 µs trajectories.
- **1.6 UQ without AI**: parameter sweeps, Sobol global sensitivity (SALib), Monte-Carlo ensembles over rate constants. Every output gets mean ± SD or confidence intervals from sampling. MAPIE conformal prediction acceptable only as a post-hoc non-parametric wrapper around physics predictions.
- **1.7 Validation harness**: PDBBind pose recovery, CASF benchmarks, ChEMBL IC50 ranking on held-out sets. Mass-balance and energy-conservation gates.

Exit criteria: reproduce known ligand / protein MD equilibration, docking poses within 2 Å RMSD on ≥ 75 % of test set, basic ΔG trends, UQ calibration error ≤ 10 %.

**Campaign 2: Mesoscale Cellular Modules (Years 2–4)**
Connect physics to biology using mechanistic rules (no learned models).

- Central Dogma at scale: 2 000–5 000 real genes (GENCODE + ReMap GRN). Gillespie / tau-leaping stochastic simulation for transcription / translation with pausing, splicing, ribosome competition. Deterministic fallback for reproducibility.
- Metabolism: full Recon3D / 4 core with dynamic flux-balance analysis (COBRApy, no ML). Compartmental redox (cytosol vs mito NADH / GSH), lipid metabolism, one-carbon, nucleotide pools.
- Signalling & cell cycle: hybrid ODE + stochastic modules for MAPK, PI3K / AKT, Wnt, TGF-β, JAK / STAT, Notch, Hedgehog, Rho GTPases. p53-MDM2-PUMA-BAX axis with real delay topology (Geva-Zatorsky).
- Apoptosis & DNA repair: rule-based cascades with mass-action kinetics, stochastic bursts for single-cell variability. Fix caspase-3 execution so real cell death occurs.
- Vesicular traffic, cytoskeleton, membrane transport: agent-based rules with physics-informed forces (Hertz contact, Langevin).
- 5–8 cell types: parameter files per type loaded from `data/`.

Validation: reproduce literature phenotypes (p53 pulsing period 5–6 h, cisplatin G2/M arrest, apoptosis commitment) across multiple cell types with external blind tests.

**Campaign 3: Whole-Cell Integration + 3D Tissue (Years 3–5)**
- Hybrid solver: synchronise modules with shared state variables and adaptive time-stepping (RK45 with error control for stiff ODEs).
- 3D spatial model: lattice or off-lattice agent-based for cell shape, ECM fibres, focal adhesions, mechanotransduction (YAP / TAZ via force input).
- Hypoxic cores, angiogenesis, immune co-culture (macrophages, T-cells).
- Mass conservation and energy balance enforced globally at every step.

Validation: spheroid growth curves, drug response in 3D, immune infiltration benchmarks from partner labs.

**Campaign 4: Adversarial Validation + Active Learning (Years 4–6)**
- Red-team suite: external chemists design compounds to break the model; add failures as new benchmarks.
- Active-learning loop (non-AI): use sensitivity analysis to propose the most informative wet-lab experiment (e.g., "measure this rate under these conditions"). Partner labs run it; result refines parameters.
- Patient variability: sample from 1000 Genomes + GTEx for genetic backgrounds.
- Uncertainty: full ensembles + Sobol indices for every prediction.

**Campaign 5: Production & Regulatory (Years 5–7+)**
- Open Docker / K8s service + YAML extension marketplace for community pathways.
- SBML export for interoperability.
- First MIDD-style package for a narrow use case (e.g., p53-axis ranking).
- Cross-lab reproducibility study.

### Key Non-AI Techniques
- **Stochastic simulation**: Gillespie SSA or tau-leaping for gene expression, replication, apoptosis bursts.
- **Rule-based / hybrid integration**: Covert / Karr style — each process has its own math (ODE, Boolean, agent), linked by shared variables.
- **Flux-balance analysis**: deterministic LP for metabolism (COBRApy).
- **Agent-based modelling**: cytoskeleton, vesicular traffic, cell migration, immune interactions.
- **PDE solvers**: nutrient / oxygen diffusion, pH gradients, morphogen fields.
- **Parameter sweeps & sensitivity**: Sobol, Monte-Carlo for uncertainty instead of neural ensembles.

### Practical Discipline
- Start with Campaign 1 layers (chem + MD + docking). Add rigorous gates: mass conservation, energy conservation, reproduction of known short MD trajectories.
- For every new module, require a literature-cited mechanistic representation and at least one external blind validation point.
- Use the team for curation (literature rates, pathway rules) and wet-lab validation, not for writing AI code.
- Keep everything interpretable: every rate must have a PMID or physics derivation.

This approach yields a simulator that is **transparent, physics-grounded, and useful for triage**. It will not replace the dish for novel entities, but it can dramatically reduce wasted experiments and focus wet-lab effort on the hardest cases.
