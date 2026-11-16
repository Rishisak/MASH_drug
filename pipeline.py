from pathlib import Path
import json, re, shutil, subprocess, urllib.request, uuid
from datetime import datetime, timezone

from Bio.PDB import MMCIFParser, MMCIF2Dict, PDBIO, Select, is_aa
from pdbfixer import PDBFixer
from openmm.app import PDBFile

class ProteinChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    def accept_residue(self, residue):
        return is_aa(residue, standard=False)

def _count_pdb(path):
    atoms, residues, chains = 0, set(), set()
    with open(path, "r", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms += 1
                chains.add(line[21:22].strip())
                residues.add((line[21:22], line[22:27], line[17:20]))
    return {"atom_records": atoms, "residues": len(residues), "chains": sorted(chains)}

def run_preprocessing(pdb_id, chain_id, output_root):
    job = f"{pdb_id}_{chain_id}_{uuid.uuid4().hex[:8]}"
    out = Path(output_root) / job
    out.mkdir(parents=True)
    report = {
        "job_id": job, "pdb_id": pdb_id, "chain_id": chain_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "checks": [], "warnings": []
    }

    cif_path = out / f"{pdb_id}.cif"
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    urllib.request.urlretrieve(cif_url, cif_path)
    report["checks"].append({"name": "PDB/mmCIF download", "status": "PASS", "detail": f"Downloaded {pdb_id}.cif from RCSB PDB."})

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cif_path))
    model = next(structure.get_models(), None)
    if model is None or chain_id not in model:
        raise ValueError(f"Chain {chain_id} was not found in the first model of {pdb_id}.")
    chain_file = out / "extracted_chain_protein.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(chain_file), ProteinChainSelect(chain_id))
    counts_in = _count_pdb(chain_file)
    if counts_in["atom_records"] == 0:
        raise ValueError("No standard amino-acid residues were extracted for this chain.")
    report["checks"].append({"name": "Chain extraction", "status": "PASS", "detail": f"Extracted chain {chain_id}: {counts_in['residues']} amino-acid residues and {counts_in['atom_records']} atom records."})

    # Inspect the original CIF for non-protein components. They are deliberately
    # not merged into the generic protein-only PDBFixer repair step.
    try:
        cif_data = MMCIF2Dict(str(cif_path))
        groups = cif_data.get("_atom_site.group_PDB", [])
        comp_ids = cif_data.get("_atom_site.label_comp_id", [])
        auth_chains = cif_data.get("_atom_site.auth_asym_id", [])
        label_chains = cif_data.get("_atom_site.label_asym_id", [])
        hetero_names = set()
        for i, group in enumerate(groups):
            comp = comp_ids[i] if i < len(comp_ids) else ""
            auth_chain = auth_chains[i] if i < len(auth_chains) else ""
            label_chain = label_chains[i] if i < len(label_chains) else ""
            if group == "HETATM" and comp not in {"HOH", "WAT"} and chain_id in {auth_chain, label_chain}:
                hetero_names.add(comp)
        report["detected_nonprotein_components"] = sorted(hetero_names)
        if hetero_names:
            report["warnings"].append(
                "Non-protein components were detected in the source structure: "
                + ", ".join(sorted(hetero_names))
                + ". They are not automatically inserted into the generic repaired PDB. "
                  "Cofactor/metal handling must be reviewed separately."
            )
    except Exception:
        report["warnings"].append("Could not automatically inventory all non-protein components.")

    repaired = out / "repaired_protein.pdb"
    fixer = PDBFixer(filename=str(chain_file))
    fixer.findMissingResidues()
    missing_residue_blocks = len(fixer.missingResidues)
    # Avoid inventing unresolved loops automatically in this first MVP.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstandard = len(fixer.nonstandardResidues)
    fixer.findMissingAtoms()
    missing_atom_count = sum(len(v) for v in fixer.missingAtoms.values()) + sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(repaired, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    counts_out = _count_pdb(repaired)
    report["repair_summary"] = {
        "missing_residue_blocks_detected_not_built": missing_residue_blocks,
        "nonstandard_residues_detected": nonstandard,
        "missing_heavy_atoms_added": missing_atom_count,
        "hydrogen_pH": 7.4
    }
    report["checks"].append({"name": "Missing heavy-atom repair", "status": "PASS", "detail": f"PDBFixer added {missing_atom_count} missing heavy atoms/terminal atoms where identified."})
    report["checks"].append({"name": "Hydrogen addition", "status": "PASS", "detail": "Hydrogens added using PDBFixer at pH 7.4."})
    status = "WARNING" if missing_residue_blocks or nonstandard else "PASS"
    report["checks"].append({"name": "Missing residues / nonstandard residues", "status": status, "detail": f"{missing_residue_blocks} unresolved residue blocks detected (not rebuilt); {nonstandard} nonstandard residues detected (not automatically substituted)."})

    report["output_structure"] = counts_out
    report["checks"].append({"name": "Output structure integrity", "status": "PASS" if counts_out["atom_records"] > 0 else "FAIL", "detail": f"Output contains {counts_out['residues']} residues and {counts_out['atom_records']} atom records."})

    # Meeko receptor conversion. Use the installed command-line executable if available.
    receptor_base = out / "prepared_receptor"
    pdbqt = out / "prepared_receptor.pdbqt"
    meeko_exe = shutil.which("mk_prepare_receptor")
    if not meeko_exe:
        # Windows venv executable fallback
        candidates = [
            Path.cwd() / "mash_env" / "Scripts" / "mk_prepare_receptor.exe",
            Path(__file__).resolve().parent.parent / "mash_env" / "Scripts" / "mk_prepare_receptor.exe"
        ]
        meeko_exe = str(next((p for p in candidates if p.exists()), "")) or None

    if meeko_exe:
        cmd = [meeko_exe, "--read_pdb", str(repaired), "-o", str(receptor_base), "-p"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0 and pdbqt.exists():
            report["checks"].append({"name": "PDBQT receptor generation", "status": "PASS", "detail": "Meeko generated prepared_receptor.pdbqt."})
        else:
            report["checks"].append({"name": "PDBQT receptor generation", "status": "WARNING", "detail": (proc.stderr or proc.stdout or "Meeko did not create the expected file.")[-1200:]})
    else:
        report["checks"].append({"name": "PDBQT receptor generation", "status": "WARNING", "detail": "Meeko executable not found on PATH. Activate the mash_env environment or add mk_prepare_receptor.exe to PATH."})

    report["files"] = []
    for p in out.iterdir():
        if p.is_file() and p.name != "report.json":
            report["files"].append({"name": p.name, "url": f"/download/{job}/{p.name}", "bytes": p.stat().st_size})
    report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    report["files"].append({"name": "report.json", "url": f"/download/{job}/report.json", "bytes": (out / "report.json").stat().st_size})
    return report
