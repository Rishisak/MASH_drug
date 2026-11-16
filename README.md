# ProteinPrep MVP

A local Flask website that retrieves a PDB structure from RCSB, extracts one protein chain, repairs missing heavy atoms, adds hydrogens at pH 7.4, runs basic output checks, and attempts Meeko PDBQT generation.

## Important scientific scope

- This MVP is protein-chain focused. It does not yet automatically reconstruct or retain every cofactor, metal, crystallographic ligand, alternate conformation, or covalent modification.
- Missing residue blocks are reported but deliberately not rebuilt automatically. Loop modelling requires additional choices and validation.
- Nonstandard residues are reported, not silently substituted.
- A PDBQT file is only generated if Meeko is available and accepts the repaired protein.
- Passing these checks does not establish docking accuracy or biological activity.

## Run on Windows (inside your existing project environment)

1. Extract this folder into `C:\Users\<you>\OneDrive\Desktop\MASH_Drug_Discovery\ProteinPrep`.
2. Activate the existing environment from the parent project:
   ```powershell
   cd C:\Users\<you>\OneDrive\Desktop\MASH_Drug_Discovery
   .\mash_env\Scripts\Activate.ps1
   ```
3. Install Flask if needed:
   ```powershell
   python -m pip install Flask
   ```
4. Start the app:
   ```powershell
   cd .\ProteinPrep
   python app.py
   ```
5. Open `http://127.0.0.1:5000` in your browser.

The first example is PDB ID `9AV4`, chain `A`.

## Output

Each run creates a uniquely named folder under `outputs/` containing the downloaded mmCIF, extracted chain PDB, repaired protein PDB, optional receptor PDBQT, and JSON validation report.

## Next development steps

1. Add a robust ligand/cofactor inventory using the mmCIF Chemical Component Dictionary.
2. Add explicit cofactor templates and preserve NAD for HSD17B13 9AV4.
3. Add PDBQT atom-type/charge validation and geometric checks.
4. Add file upload for PDB and mmCIF.
5. Add downloadable ZIP bundles and run history.
