# Third-party notices

CellSim itself is MIT licensed (see LICENSE). Portions of this software incorporate or adapt code, models or data from the projects below; their licenses apply to those portions.


- PhysiPKPD (BSD-3-Clause License)
  Bergman et al., "PhysiPKPD: A pharmacokinetics and pharmacodynamics
  module for PhysiCell", GigaByte 2023.
  https://github.com/drbergman/PhysiPKPD

- PhysiCell (BSD-3-Clause License)
  Ghaffarizadeh et al., PLoS Computational Biology 2018.
  https://github.com/MathCancer/PhysiCell

- Drug parameters derived from BioModels Database (CC0 Public Domain)
  https://www.ebi.ac.uk/biomodels/

- `benchmarks/cell/gdsc_reference.csv` is a 30-row extract of fitted IC50
  and AUC values from Genomics of Drug Sensitivity in Cancer (GDSC1 and
  GDSC2, release 8.4), Wellcome Sanger Institute.
  Yang et al., Nucleic Acids Research 41:D955 (2013);
  Iorio et al., Cell 166:740 (2016). https://www.cancerrxgene.org/
  The full tables are not redistributed; `scripts/gdsc_reference.py`
  downloads them from the Sanger FTP site.

- Benchmark structures under `benchmarks/` are from the RCSB Protein
  Data Bank (https://www.rcsb.org/), whose data are in the public domain
  (CC0); FreeSolv hydration free energies are from Mobley & Guthrie,
  J Comput Aided Mol Des 28:711 (2014), CC-BY 4.0.
