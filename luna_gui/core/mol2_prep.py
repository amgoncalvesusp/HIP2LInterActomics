"""Split MOL2/PDB complexes (protein + ligand) into separate files.

Consolidated from the project's original molecular-preparation prototype. Given
a folder of
MOL2 files where each file contains a protein (atoms 1..last_pa) followed by a
ligand (atoms > last_pa), or PDB complex files containing protein, ligand and
optional waters, produces:
  <out_dir>/proteinas_pdb/<name>.pdb
  <out_dir>/ligantes_mol2/<name>_ligand.mol2
  <out_dir>/ligantes_sdf/<name>_ligand.sdf
"""
from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# Standard residue names recognised as belonging to the protein/solvent.
# Ligand atoms are anything NOT in this set.
# ---------------------------------------------------------------------------
_PROTEIN_RESIDUES: frozenset[str] = frozenset({
    # Standard amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Histidine / cysteine protonation variants (common in GOLD/CHARMM)
    "ASH", "GLH", "LYN", "HSD", "HSE", "HSP", "HIE", "HID", "HIP", "CYX", "CYM",
    # N/C terminus caps
    "ACE", "NME", "NHE",
    # Water / solvent
    "HOH", "WAT", "WTM", "TIP", "SOL", "T3P",
    # Selenium / non-standard but common
    "MSE", "SEC", "PYL", "SEP", "TPO", "PTR", "CSO", "HYP",
})

_WATER_RESIDUES: frozenset[str] = frozenset({
    "HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD",
})

_VALID_ELEMENTS: frozenset[str] = frozenset({
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi",
})

_COVALENT_RADII: dict[str, float] = {
    "H": 0.31,
    "B": 0.85,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


_EXTERNAL_PDB_TO_SDF_SCRIPT = r"""
import json
import sys


def normalize_sdf_record(sdf_text, molecule_name, properties=None):
    text = sdf_text.replace("\r\n", "\n").replace("\r", "\n").strip()
    lines = text.splitlines()
    if lines:
        lines[0] = molecule_name
    else:
        lines = [molecule_name, "  HIP2LInterActomics", "", "M  END"]
    text = "\n".join(lines).rstrip()
    if text.endswith("$$$$"):
        text = text[:-4].rstrip()
    if properties:
        for key, value in properties.items():
            text += f"\n>  <{key}>\n{value}\n"
    text += "\n$$$$"
    return text + "\n"


def apply_formal_charges(mol, formal_charges):
    for atom, charge in zip(mol.GetAtoms(), formal_charges):
        if charge:
            atom.SetFormalCharge(int(charge))
    mol.UpdatePropertyCache(strict=False)


def gasteiger_charges(mol):
    values = []
    for atom in mol.GetAtoms():
        try:
            charge = atom.GetProp("_GasteigerCharge")
        except Exception:
            charge = ""
        if not charge or charge.lower() in {"nan", "-nan", "inf", "-inf"}:
            charge = "0.000000"
        values.append(f"{atom.GetIdx() + 1}:{float(charge):.6f}")
    return " ".join(values)


def main():
    payload = json.loads(sys.stdin.read() or "{}")
    molecule_name = payload["molecule_name"]
    pdb_block = payload["pdb_block"]
    formal_charges = [int(value or 0) for value in payload.get("formal_charges") or []]
    warnings = []
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdDetermineBonds
    except Exception as exc:
        print(json.dumps({"ok": False, "error": f"RDKit indisponivel no Python externo: {type(exc).__name__}: {exc}"}))
        return

    mol = Chem.MolFromPDBBlock(
        pdb_block,
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    if mol is None or mol.GetNumAtoms() == 0:
        print(json.dumps({"ok": False, "error": "RDKit nao conseguiu ler o bloco PDB do ligante."}))
        return

    mol.SetProp("_Name", molecule_name)
    total_charge = int(sum(formal_charges))
    try:
        rdDetermineBonds.DetermineBonds(
            mol,
            charge=total_charge,
            allowChargedFragments=True,
            embedChiral=True,
        )
    except Exception as exc:
        warnings.append(
            f"{molecule_name}: RDKit DetermineBonds falhou ({type(exc).__name__}); usando conectividade PDB/proximidade."
        )
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is None or mol.GetNumAtoms() == 0:
            print(json.dumps({"ok": False, "error": "RDKit nao conseguiu recuperar conectividade por proximidade."}))
            return
        mol.SetProp("_Name", molecule_name)

    try:
        apply_formal_charges(mol, formal_charges)
        Chem.SanitizeMol(mol)
        Chem.SetAromaticity(mol)
    except Exception as exc:
        print(json.dumps({"ok": False, "error": f"Falha na sanitizacao/aromaticidade RDKit: {type(exc).__name__}: {exc}"}))
        return

    charge_text = ""
    try:
        AllChem.ComputeGasteigerCharges(mol)
        charge_text = gasteiger_charges(mol)
    except Exception as exc:
        warnings.append(
            f"{molecule_name}: RDKit converteu o ligante, mas nao calculou cargas Gasteiger ({type(exc).__name__})."
        )

    sdf_text = Chem.MolToMolBlock(mol, kekulize=False)
    if not sdf_text or "M  END" not in sdf_text:
        print(json.dumps({"ok": False, "error": "RDKit nao gerou um bloco SDF valido."}))
        return

    props = {
        "HIP2L_ChargeModel": "formal PDB charges + RDKit Gasteiger",
        "HIP2L_TotalFormalCharge": str(total_charge),
    }
    if charge_text:
        props["HIP2L_GasteigerCharges"] = charge_text
    print(json.dumps({
        "ok": True,
        "sdf": normalize_sdf_record(sdf_text, molecule_name, props),
        "warnings": warnings,
    }))


main()
"""


@dataclass
class DetectionResult:
    last_pa: int               # detected last protein atom index
    method: str                # "substructure" | "atom_scan" | "fallback"
    n_protein_atoms: int
    n_ligand_atoms: int
    ligand_names: list[str]    # subst_names detected as ligand


def detect_last_protein_atom(mol2_path: str | Path) -> DetectionResult:
    """Heuristically determine the last protein atom index in a combined MOL2.

    Strategy (tried in order):
    1. **SUBSTRUCTURE section** — most reliable: reads the @<TRIPOS>SUBSTRUCTURE
       block, identifies which substructures are ligands by name, and returns
       the atom index just before the first ligand substructure.
    2. **ATOM scan** — fallback: walks ATOM records comparing the first 3 chars
       of each atom's subst_name against the known protein-residue set.
    """
    lines = Path(mol2_path).read_text(errors="replace").splitlines()
    result = _detect_via_substructure(lines)
    if result is not None:
        return result
    return _detect_via_atom_scan(lines)


def _detect_via_substructure(lines: list[str]) -> DetectionResult | None:
    """Parse @<TRIPOS>SUBSTRUCTURE to find the ligand boundary."""
    # Collect substructures: list of (first_atom_id, subst_name)
    in_sub = False
    substrs: list[tuple[int, str]] = []
    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>SUBSTRUCTURE":
            in_sub = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_sub:
            break
        if not in_sub or not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 3:
            continue
        try:
            first_atom = int(parts[2])
            name = parts[1]
            substrs.append((first_atom, name))
        except ValueError:
            continue

    if not substrs:
        return None

    # Sort by first atom index
    substrs.sort(key=lambda t: t[0])

    # Determine which substructures are protein vs ligand
    # A substructure is "protein" if its 3-letter prefix is in _PROTEIN_RESIDUES
    lig_names: list[str] = []
    lig_start: int | None = None
    for first_atom, name in substrs:
        if not _is_protein_or_solvent_substructure(name):
            if lig_start is None:
                lig_start = first_atom
            lig_names.append(name)

    if lig_start is None:
        # No ligand found — treat everything as protein
        last_atom = substrs[-1][0]  # rough; no real ligand boundary
        return DetectionResult(
            last_pa=last_atom, method="substructure",
            n_protein_atoms=last_atom, n_ligand_atoms=0, ligand_names=[],
        )

    last_pa = lig_start - 1
    n_protein, n_ligand = _count_atoms_by_role(lines)
    return DetectionResult(
        last_pa=last_pa, method="substructure",
        n_protein_atoms=n_protein,
        n_ligand_atoms=n_ligand,
        ligand_names=list(dict.fromkeys(lig_names)),  # deduplicated, ordered
    )


def _detect_via_atom_scan(lines: list[str]) -> DetectionResult:
    """Walk ATOM records; last atom whose subst_name prefix is in protein set."""
    in_atom = False
    first_lig = None
    last_prot_before_lig = 0
    total = 0
    n_protein = 0
    n_ligand = 0
    lig_names: list[str] = []
    seen_lig: set[str] = set()

    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atom = True; continue
        if stripped.startswith("@<TRIPOS>") and in_atom:
            break
        if not in_atom or not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        try:
            idx = int(parts[0])
        except ValueError:
            continue
        total = max(total, idx)
        atom_name = parts[1] if len(parts) > 1 else ""
        atom_type = parts[5] if len(parts) > 5 else ""
        subst_name = parts[7] if len(parts) > 7 else ""
        if _is_water_substructure(subst_name) and _is_lone_pair(atom_name, atom_type):
            continue
        if _is_protein_or_solvent_substructure(subst_name):
            n_protein += 1
            if first_lig is None:
                last_prot_before_lig = max(last_prot_before_lig, idx)
        else:
            n_ligand += 1
            if first_lig is None:
                first_lig = idx
            if subst_name and subst_name not in seen_lig:
                seen_lig.add(subst_name)
                lig_names.append(subst_name)

    last_pa = (first_lig - 1) if first_lig is not None else (last_prot_before_lig if last_prot_before_lig > 0 else max(1, total - 1))
    return DetectionResult(
        last_pa=last_pa, method="atom_scan" if last_prot_before_lig > 0 else "fallback",
        n_protein_atoms=n_protein if n_protein else last_pa,
        n_ligand_atoms=n_ligand if n_ligand else max(0, total - last_pa),
        ligand_names=lig_names,
    )


def _count_atoms(lines: list[str]) -> int:
    in_atom = False
    count = 0
    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atom = True; continue
        if stripped.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            if parts and parts[0].isdigit():
                count += 1
    return count


def _count_atoms_by_role(lines: list[str]) -> tuple[int, int]:
    in_atom = False
    n_protein = 0
    n_ligand = 0
    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atom = True; continue
        if stripped.startswith("@<TRIPOS>") and in_atom:
            break
        if not in_atom or not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 6 or not parts[0].isdigit():
            continue
        atom_name = parts[1] if len(parts) > 1 else ""
        atom_type = parts[5] if len(parts) > 5 else ""
        subst_name = parts[7] if len(parts) > 7 else ""
        if _is_water_substructure(subst_name) and _is_lone_pair(atom_name, atom_type):
            continue
        if _is_protein_or_solvent_substructure(subst_name):
            n_protein += 1
        else:
            n_ligand += 1
    return n_protein, n_ligand


@dataclass
class PrepResult:
    files_processed: int
    proteins_written: int
    ligands_written: int
    water_molecules_detected: int
    protein_dir: str
    ligand_dir: str
    errors: list[str]
    log_file: str = ""


def split_docking_folder(
    src_folder: str | Path,
    last_pa: int,
    out_folder: str | Path | None = None,
    protein_subdir: str = "proteinas_pdb",
    ligand_subdir: str = "ligantes_mol2",
    progress_cb=None,
) -> PrepResult:
    """Split every .mol2 in `src_folder` into PDB (protein) + MOL2 (ligand)."""
    src = Path(src_folder)
    out = Path(out_folder) if out_folder else src
    prot_dir = out / protein_subdir
    lig_dir = out / ligand_subdir
    prot_dir.mkdir(parents=True, exist_ok=True)
    lig_dir.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    nfiles = nprot = nlig = nwater = 0
    first_lig_id = last_pa + 1

    mol2_files = sorted(src.glob("*.mol2"))
    total_files = len(mol2_files)

    if progress_cb is not None:
        progress_cb(0, total_files, "", True, "")

    for mol2 in mol2_files:
        nfiles += 1
        ok = True
        error_message = ""
        try:
            with mol2.open("r", errors="replace") as fh:
                lines = fh.readlines()
            prot_atoms, lig_atoms, lig_bonds, water_count = _split_one(lines, last_pa, first_lig_id)
            nwater += water_count
        except Exception as e:
            ok = False
            error_message = str(e)
            errors.append(f"{mol2.name}: {e}")
            if progress_cb is not None:
                progress_cb(nfiles, total_files, mol2.name, ok, error_message)
            continue

        name = mol2.stem
        if prot_atoms:
            _write_pdb(prot_dir / f"{name}.pdb", prot_atoms)
            nprot += 1
        if lig_atoms:
            _write_mol2(lig_dir / f"{name}_ligand.mol2", name, lig_atoms, lig_bonds)
            nlig += 1
        if progress_cb is not None:
            progress_cb(nfiles, total_files, mol2.name, ok, error_message)

    log_file = _write_prep_log(
        out,
        "MOL2",
        nfiles,
        nprot,
        nlig,
        nwater,
        errors,
    )
    return PrepResult(
        files_processed=nfiles, proteins_written=nprot, ligands_written=nlig,
        water_molecules_detected=nwater,
        protein_dir=str(prot_dir), ligand_dir=str(lig_dir), errors=errors,
        log_file=str(log_file),
    )


def split_complex_folder(
    src_folder: str | Path,
    last_pa: int | None = None,
    out_folder: str | Path | None = None,
    progress_cb=None,
    chemistry_python: str | Path | None = None,
) -> PrepResult:
    """Split a homogeneous folder of MOL2 or PDB complexes.

    MOL2 inputs use the historical docking splitter. PDB/ENT inputs are split by
    record/residue role: ATOM/protein HETATM records go to the receptor PDB,
    water records stay with the receptor, and non-protein ATOM/HETATM records
    become one ligand SDF per complex.
    """
    src = Path(src_folder)
    mol2_files = sorted(src.glob("*.mol2"))
    pdb_files = sorted(
        candidate
        for candidate in src.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".pdb", ".ent"}
    )
    if mol2_files and pdb_files:
        raise ValueError(
            "A pasta contem arquivos MOL2 e PDB/ENT. Use uma pasta com apenas "
            "um formato de complexo por preparacao."
        )
    if mol2_files:
        if last_pa is None:
            last_pa = detect_last_protein_atom(mol2_files[0]).last_pa
        return split_docking_folder(src, int(last_pa), out_folder, progress_cb=progress_cb)
    if pdb_files:
        return split_pdb_complex_folder(
            src,
            out_folder,
            progress_cb=progress_cb,
            chemistry_python=chemistry_python,
        )
    raise ValueError("Nenhum arquivo .mol2, .pdb ou .ent encontrado na pasta de origem.")


def split_pdb_complex_folder(
    src_folder: str | Path,
    out_folder: str | Path | None = None,
    protein_subdir: str = "proteinas_pdb",
    ligand_subdir: str = "ligantes_sdf",
    progress_cb=None,
    chemistry_python: str | Path | None = None,
) -> PrepResult:
    """Split every PDB/ENT complex into receptor PDB + ligand SDF.

    PDB is useful for separating records, but it is a poor final ligand format
    for LUNA because bond perception, charge and aromaticity are ambiguous.  The
    separated ligand is therefore written as SDF, which LUNA reads through the
    RDKit backend.  Conversion uses Open Babel or RDKit chemistry perception so
    charge and aromaticity are handled explicitly instead of silently writing a
    geometry-only ligand PDB.
    """
    src = Path(src_folder)
    out = Path(out_folder) if out_folder else src
    prot_dir = out / protein_subdir
    lig_dir = out / ligand_subdir
    prot_dir.mkdir(parents=True, exist_ok=True)
    lig_dir.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    conversion_warnings: list[str] = []
    nfiles = nprot = nlig = nwater = 0
    pdb_files = sorted(
        candidate
        for candidate in src.iterdir()
        if candidate.is_file() and candidate.suffix.lower() in {".pdb", ".ent"}
    )
    total_files = len(pdb_files)

    if progress_cb is not None:
        progress_cb(0, total_files, "", True, "")

    for pdb_file in pdb_files:
        nfiles += 1
        ok = True
        error_message = ""
        try:
            lines = pdb_file.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
            protein_lines, ligand_lines, water_count, split_warnings = _split_pdb_complex(lines)
            nwater += water_count
            conversion_warnings.extend(f"{pdb_file.name}: {warning}" for warning in split_warnings)
        except Exception as exc:
            ok = False
            error_message = str(exc)
            errors.append(f"{pdb_file.name}: {exc}")
            if progress_cb is not None:
                progress_cb(nfiles, total_files, pdb_file.name, ok, error_message)
            continue

        name = pdb_file.stem
        if protein_lines:
            _write_pdb_records(prot_dir / f"{name}.pdb", protein_lines)
            nprot += 1
        if ligand_lines:
            _write_ligand_sdf_from_pdb_records(
                lig_dir / f"{name}_ligand.sdf",
                f"{name}_LIG",
                ligand_lines,
                warnings=conversion_warnings,
                chemistry_python=chemistry_python,
            )
            nlig += 1
        if progress_cb is not None:
            progress_cb(nfiles, total_files, pdb_file.name, ok, error_message)

    all_messages = errors + conversion_warnings
    log_file = _write_prep_log(
        out,
        "PDB/ENT",
        nfiles,
        nprot,
        nlig,
        nwater,
        all_messages,
    )
    return PrepResult(
        files_processed=nfiles,
        proteins_written=nprot,
        ligands_written=nlig,
        water_molecules_detected=nwater,
        protein_dir=str(prot_dir),
        ligand_dir=str(lig_dir),
        errors=all_messages,
        log_file=str(log_file),
    )


def _write_prep_log(
    out_folder: Path,
    input_format: str,
    files_processed: int,
    proteins_written: int,
    ligands_written: int,
    water_molecules_detected: int,
    messages: list[str],
) -> Path:
    path = Path(out_folder) / "preprocess.log"
    lines = [
        "HIP2LInterActomics preprocessing log\n",
        f"Formato: {input_format}\n",
        f"Arquivos lidos: {files_processed}\n",
        f"Proteinas escritas: {proteins_written}\n",
        f"Ligantes escritos: {ligands_written}\n",
        f"Aguas detectadas: {water_molecules_detected}\n",
    ]
    if messages:
        lines.append("\nMensagens:\n")
        lines.extend(f"- {message}\n" for message in messages)
    else:
        lines.append("\nMensagens: nenhuma\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _split_pdb_complex(lines: list[str]) -> tuple[list[str], list[str], int, list[str]]:
    protein_lines: list[str] = []
    ligand_lines: list[str] = []
    waters: set[tuple[str, str, str, str]] = set()
    ligand_serials: set[int] = set()
    warnings: set[str] = set()

    for raw_line in lines:
        if not (raw_line.startswith("ATOM") or raw_line.startswith("HETATM")):
            continue
        line = raw_line if raw_line.endswith("\n") else raw_line + "\n"
        record = line[:6].strip().upper()
        resname = line[17:20].strip() if len(line) >= 20 else ""
        atom_name = line[12:16].strip() if len(line) >= 16 else ""
        element = line[76:78].strip() if len(line) >= 78 else ""
        serial = _pdb_atom_serial(line)
        is_water = _is_water_substructure(resname)
        if is_water and _is_lone_pair(atom_name, element):
            continue
        if is_water:
            waters.add(_pdb_water_key(line))
            if record != "HETATM":
                warnings.add(
                    f"residuo de agua {resname or 'UNK'} encontrado como {record}; mantido junto da proteina."
                )
            protein_lines.append(line)
            continue
        is_protein_residue = _is_protein_residue_name(resname)
        if is_protein_residue:
            if record == "HETATM":
                warnings.add(
                    f"residuo de aminoacido especial {resname or 'UNK'} encontrado como HETATM; mantido na proteina."
                )
            protein_lines.append(line)
            continue
        if _is_lone_pair(atom_name, element):
            continue
        if record == "ATOM":
            warnings.add(
                f"residuo {resname or 'UNK'} encontrado como ATOM, mas nao e aminoacido/agua; tratado como ligante."
            )
        if serial is not None:
            ligand_serials.add(serial)
        ligand_lines.append(line)

    if ligand_serials:
        ligand_lines.extend(_ligand_conect_records(lines, ligand_serials))

    return protein_lines, ligand_lines, len(waters), sorted(warnings)


def _split_one(
    lines: list[str], last_pa: int, first_lig_id: int,
) -> tuple[list[str], list[tuple[int, str]], list[tuple[int, int, int, str]], int]:
    """Split raw MOL2 lines into protein PDB records + ligand atoms/bonds."""
    section: str | None = None
    prot_atoms: list[str] = []
    lig_atoms: list[tuple[int, str]] = []
    lig_bonds: list[tuple[int, int, int, str]] = []
    ligand_atom_map: dict[int, int] = {}
    water_residues: set[str] = set()
    chain_by_subst = _parse_substructure_chains(lines)

    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            section = "ATOM"; continue
        if "@<TRIPOS>BOND" in line:
            section = "BOND"; continue
        if "@<TRIPOS>SUBSTRUCTURE" in line:
            section = "SUB"; continue
        if line.startswith("@<TRIPOS>"):
            section = None; continue

        if section == "ATOM":
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                idx = int(parts[0])
            except ValueError:
                continue

            atom_name = parts[1]
            atom_type = parts[5] if len(parts) > 5 else ""
            subst_id = int(parts[6]) if len(parts) > 6 and parts[6].lstrip("-").isdigit() else 0
            subst_name = parts[7] if len(parts) > 7 else "UNK1"
            is_water = _is_water_substructure(subst_name)
            if is_water:
                water_residues.add(_mol2_water_key(subst_id, subst_name))
            if is_water and _is_lone_pair(atom_name, atom_type):
                continue

            if idx <= last_pa or is_water:
                # Protein/solvent atom -> PDB format.
                try:
                    resname, resseq, icode = _parse_subst_name(subst_name)
                    chain_id = _clean_pdb_chain_id(
                        chain_by_subst.get(subst_id, "A"),
                        fallback="W" if is_water else "A",
                    )
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    element = _infer_element(atom_name, atom_type)
                    if element is None:
                        continue
                    pdb_line = _format_pdb_atom_line(
                        record_name="HETATM" if is_water else "ATOM",
                        idx=idx,
                        atom_name=atom_name,
                        resname=resname,
                        chain_id=chain_id,
                        resseq=resseq,
                        icode=icode,
                        x=x,
                        y=y,
                        z=z,
                        element=element,
                    )
                    prot_atoms.append(pdb_line)
                except (ValueError, IndexError):
                    continue
            else:
                # Ligand atom — keep the line intact after remapping the ID.
                new_idx = len(ligand_atom_map) + 1
                ligand_atom_map[idx] = new_idx
                idx_str = parts[0]
                end_pos = line.find(idx_str) + len(idx_str)
                rest = line[end_pos:]
                lig_atoms.append((new_idx, rest))

        elif section == "BOND":
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                a1, a2 = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if a1 in ligand_atom_map and a2 in ligand_atom_map:
                lig_bonds.append((
                    len(lig_bonds) + 1,
                    ligand_atom_map[a1],
                    ligand_atom_map[a2],
                    parts[3],
                ))

    return prot_atoms, lig_atoms, lig_bonds, len(water_residues)


def count_water_molecules(path: str | Path) -> int:
    """Count water residues/molecules in a PDB/MOL2 file or input folder."""
    p = Path(path)
    if not str(path or "").strip() or not p.exists():
        return 0
    if p.is_dir():
        total = 0
        for candidate in sorted(p.iterdir()):
            if candidate.is_file() and candidate.suffix.lower() in {".pdb", ".ent", ".mol2"}:
                total += count_water_molecules(candidate)
        return total

    suffix = p.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        return _count_water_molecules_in_pdb(p)
    if suffix == ".mol2":
        return _count_water_molecules_in_mol2(p)
    return 0


def count_water_molecules_in_inputs(*paths: str | Path) -> int:
    """Count waters across multiple input paths, ignoring empty/nonexistent paths."""
    return sum(count_water_molecules(path) for path in paths if str(path or "").strip())


def _count_water_molecules_in_pdb(path: Path) -> int:
    waters: set[tuple[str, str, str, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip() if len(line) >= 20 else ""
            if not _is_water_substructure(resname):
                continue
            chain_id = line[21:22].strip() if len(line) >= 22 else ""
            resseq = line[22:26].strip() if len(line) >= 26 else ""
            icode = line[26:27].strip() if len(line) >= 27 else ""
            waters.add((resname.upper(), chain_id, resseq, icode))
    return len(waters)


def _count_water_molecules_in_mol2(path: Path) -> int:
    waters: set[tuple[int, str]] = set()
    molecule_idx = -1
    in_atom = False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped == "@<TRIPOS>MOLECULE":
                molecule_idx += 1
                in_atom = False
                continue
            if stripped == "@<TRIPOS>ATOM":
                in_atom = True
                continue
            if stripped.startswith("@<TRIPOS>") and in_atom:
                in_atom = False
                continue
            if not in_atom or not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 8:
                continue
            subst_id = int(parts[6]) if parts[6].lstrip("-").isdigit() else 0
            subst_name = parts[7]
            if _is_water_substructure(subst_name):
                waters.add((max(0, molecule_idx), _mol2_water_key(subst_id, subst_name)))
    return len(waters)


def _substructure_prefix(subst_name: str) -> str:
    return re.sub(r"[^A-Za-z]", "", str(subst_name or "")).upper()[:3]


def _water_residue_name(subst_name: str) -> str | None:
    raw = str(subst_name or "").strip().upper()
    compact = re.sub(r"[^A-Z0-9]", "", raw)
    for water_name in sorted(_WATER_RESIDUES, key=len, reverse=True):
        if compact.startswith(water_name):
            return water_name
    return None


def _is_water_substructure(subst_name: str) -> bool:
    return _water_residue_name(subst_name) is not None


def _is_protein_residue_name(subst_name: str) -> bool:
    """Return True only for residues that should stay in the receptor PDB.

    Some docking/MD exporters write ligands with ATOM records instead of
    HETATM.  For PDB complexes the record name alone is therefore not enough:
    the residue name is the safer discriminator.  Waters are handled separately
    so they can be counted and stripped of lone-pair pseudo atoms.
    """
    return _substructure_prefix(subst_name) in _PROTEIN_RESIDUES and not _is_water_substructure(subst_name)


def _mol2_water_key(subst_id: int, subst_name: str) -> str:
    name = str(subst_name or "").strip().upper()
    if subst_id:
        return f"id:{subst_id}:{name}"
    return f"name:{name}"


def _is_protein_or_solvent_substructure(subst_name: str) -> bool:
    return _substructure_prefix(subst_name) in _PROTEIN_RESIDUES or _is_water_substructure(subst_name)


def _is_lone_pair(atom_name: str, atom_type: str = "") -> bool:
    name = str(atom_name or "").strip().upper()
    atype = str(atom_type or "").strip().upper()
    return name == "LP" or atype == "LP" or atype.startswith("LP.")


def _parse_substructure_chains(lines: list[str]) -> dict[int, str]:
    """Map MOL2 substructure IDs to PDB chain IDs."""
    in_sub = False
    chain_by_subst: dict[int, str] = {}
    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>SUBSTRUCTURE":
            in_sub = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_sub:
            break
        if not in_sub or not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        try:
            subst_id = int(parts[0])
        except ValueError:
            continue
        chain_id = (parts[5].strip() or "A")[:1]
        chain_by_subst[subst_id] = chain_id
    return chain_by_subst


def _clean_pdb_chain_id(chain_id: str, fallback: str = "A") -> str:
    value = str(chain_id or "").strip()[:1]
    if value and re.match(r"[A-Za-z0-9]", value):
        return value
    return fallback


def _parse_subst_name(subst_name: str) -> tuple[str, int, str]:
    """Parse a MOL2 subst_name into (resname, resseq, icode) for PDB output.

    GOLD and other docking programs concatenate residue name + number into one
    field, e.g. "ALA12", "SER100A", "HOH400".  The PDB format requires them
    in separate fixed columns, and the residue sequence number MUST be an
    integer — otherwise LUNA's PDB parser crashes with 'NoneType has no
    attribute has_id'.

    Examples
    --------
    "ALA12"   → ("ALA", 12, " ")
    "SER100A" → ("SER", 100, "A")
    "HOH400"  → ("HOH", 400, " ")
    "LIG"     → ("LIG",   1, " ")
    "U1"      → ("U  ",   1, " ")
    """
    subst_name = subst_name.strip()
    if not subst_name:
        return ("UNK", 1, " ")

    water_name = _water_residue_name(subst_name)
    if water_name:
        digits = re.sub(r"[^0-9]", "", subst_name) or "1"
        return (water_name.ljust(3)[:3], int(digits[:4]), " ")

    # Pattern: 1-3 letters → resname, then digits → resseq, then optional letter → icode
    m = re.match(r'^([A-Za-z]{1,3})(\d+)([A-Za-z]?)$', subst_name)
    if m:
        resname = m.group(1).upper().ljust(3)[:3]
        resseq  = int(m.group(2)[:4])          # cap at 9999
        icode   = m.group(3) or " "
        return (resname, resseq, icode)

    # Fallback: letters-only (e.g. "LIG", "UNK")
    m2 = re.match(r'^([A-Za-z]{1,3})$', subst_name)
    if m2:
        return (m2.group(1).upper().ljust(3)[:3], 1, " ")

    # Last resort: strip non-alphanumeric, pull any digits found
    letters = re.sub(r'[^A-Za-z]', '', subst_name)[:3].upper().ljust(3)
    digits  = re.sub(r'[^0-9]', '', subst_name) or "1"
    return (letters, int(digits[:4]), " ")


def _write_pdb(path: Path, atoms: list[str]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("REMARK   Separated Protein\n")
        f.writelines(atoms)
        f.write("END\n")


def _pdb_water_key(line: str) -> tuple[str, str, str, str]:
    resname = line[17:20].strip().upper() if len(line) >= 20 else ""
    chain_id = line[21:22].strip() if len(line) >= 22 else ""
    resseq = line[22:26].strip() if len(line) >= 26 else ""
    icode = line[26:27].strip() if len(line) >= 27 else ""
    return (resname, chain_id, resseq, icode)


def _pdb_atom_serial(line: str) -> int | None:
    try:
        return int(line[6:11])
    except ValueError:
        return None


def _parse_conect_numbers(line: str) -> list[int]:
    numbers: list[int] = []
    for start in range(6, len(line), 5):
        token = line[start:start + 5].strip()
        if not token:
            continue
        try:
            numbers.append(int(token))
        except ValueError:
            continue
    return numbers


def _format_conect_line(source: int, targets: list[int]) -> str:
    fields = [source] + targets[:4]
    return "CONECT" + "".join(f"{value:5d}" for value in fields) + "\n"


def _ligand_conect_records(lines: list[str], ligand_serials: set[int]) -> list[str]:
    records: list[str] = []
    seen: set[tuple[int, tuple[int, ...]]] = set()
    for raw_line in lines:
        if not raw_line.startswith("CONECT"):
            continue
        numbers = _parse_conect_numbers(raw_line)
        if len(numbers) < 2 or numbers[0] not in ligand_serials:
            continue
        targets = [number for number in numbers[1:] if number in ligand_serials and number != numbers[0]]
        if not targets:
            continue
        key = (numbers[0], tuple(targets))
        if key in seen:
            continue
        seen.add(key)
        records.append(_format_conect_line(numbers[0], targets))
    return records


def _write_pdb_records(path: Path, atoms: list[str]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("REMARK   Separated Protein\n")
        for line in atoms:
            f.write(line if line.endswith("\n") else line + "\n")
        f.write("END\n")


def _write_ligand_pdb(path: Path, molecule_name: str, atoms: list[str]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"COMPND    {molecule_name}\n")
        f.write("REMARK   Separated Ligand\n")
        for line in atoms:
            f.write(line if line.endswith("\n") else line + "\n")
        f.write("END\n")


def _write_ligand_sdf_from_pdb_records(
    path: Path,
    molecule_name: str,
    records: list[str],
    warnings: list[str] | None = None,
    chemistry_python: str | Path | None = None,
) -> None:
    atom_records = [
        line if line.endswith("\n") else line + "\n"
        for line in records
        if line.startswith(("ATOM", "HETATM"))
        and not _is_lone_pair(
            line[12:16].strip() if len(line) >= 16 else "",
            line[76:78].strip() if len(line) >= 78 else "",
        )
    ]
    if not atom_records:
        raise ValueError(f"Nenhum atomo real encontrado para o ligante {molecule_name}.")

    conect_records = [
        line if line.endswith("\n") else line + "\n"
        for line in records
        if line.startswith("CONECT")
    ]
    pdb_block = _ligand_pdb_block(molecule_name, atom_records, conect_records)
    formal_charges = [_pdb_formal_charge(line) for line in atom_records]
    path.parent.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    for writer in (_write_sdf_with_openbabel, _write_sdf_with_rdkit):
        try:
            if writer(path, molecule_name, pdb_block, formal_charges, warnings):
                return
        except Exception as exc:
            errors.append(f"{writer.__name__}: {type(exc).__name__}: {exc}")

    try:
        if _write_sdf_with_external_openbabel(
            path,
            molecule_name,
            pdb_block,
            formal_charges,
            chemistry_python,
            warnings,
        ):
            return
    except Exception as exc:
        errors.append(f"_write_sdf_with_external_openbabel: {type(exc).__name__}: {exc}")

    try:
        if _write_sdf_with_external_python(
            path,
            molecule_name,
            pdb_block,
            formal_charges,
            chemistry_python,
            warnings,
        ):
            return
    except Exception as exc:
        errors.append(f"_write_sdf_with_external_python: {type(exc).__name__}: {exc}")

    details = "; ".join(errors[:3])
    if details:
        details = f" Tentativas: {details}."
    raise ValueError(
        f"Falha ao converter {molecule_name} para SDF com percepcao de cargas/aromaticidade. "
        "Verifique a aba Inicio > Instalar LUNA para garantir RDKit/Open Babel no luna-env."
        + details
    )


def _ligand_pdb_block(molecule_name: str, atom_records: list[str], conect_records: list[str]) -> str:
    lines = [f"COMPND    {molecule_name}\n", "REMARK   Separated Ligand\n"]
    lines.extend(atom_records)
    lines.extend(conect_records)
    lines.append("END\n")
    return "".join(lines)


def _write_sdf_with_openbabel(
    path: Path,
    molecule_name: str,
    pdb_block: str,
    formal_charges: list[int],
    warnings: list[str] | None = None,
) -> bool:
    try:
        from openbabel import pybel  # type: ignore
        from openbabel import openbabel as ob  # type: ignore
    except Exception:
        return False

    mol = pybel.readstring("pdb", pdb_block)
    if mol is None:
        return False
    mol.title = molecule_name
    for idx, charge in enumerate(formal_charges, start=1):
        if charge == 0:
            continue
        atom = mol.OBMol.GetAtom(idx)
        if atom:
            atom.SetFormalCharge(int(charge))
    try:
        mol.OBMol.PerceiveBondOrders()
    except Exception:
        pass
    try:
        ob.OBAtomTyper().AssignTypes(mol.OBMol)
    except Exception:
        pass
    try:
        ob.OBBondTyper().AssignTypes(mol.OBMol)
    except Exception:
        pass
    try:
        ob.OBAromaticTyper().AssignAromaticFlags(mol.OBMol)
    except Exception:
        try:
            mol.OBMol.SetAromaticPerceived(False)
            mol.OBMol.PerceiveAromaticity()
        except Exception:
            pass
    gasteiger_charges = ""
    try:
        mol.calccharges("gasteiger")
        gasteiger_charges = _openbabel_gasteiger_charges(mol)
    except Exception as exc:
        if warnings is not None:
            warnings.append(
                f"{molecule_name}: Open Babel converteu o ligante, mas nao calculou cargas Gasteiger ({type(exc).__name__})."
            )
    sdf_text = mol.write("sdf")
    if not sdf_text or "M  END" not in sdf_text:
        return False
    props = {
        "HIP2L_ChargeModel": "formal PDB charges + OpenBabel Gasteiger",
        "HIP2L_TotalFormalCharge": str(sum(formal_charges)),
    }
    if gasteiger_charges:
        props["HIP2L_GasteigerCharges"] = gasteiger_charges
    path.write_text(_normalize_sdf_record(sdf_text, molecule_name, props), encoding="utf-8")
    return True


def _write_sdf_with_rdkit(
    path: Path,
    molecule_name: str,
    pdb_block: str,
    formal_charges: list[int],
    warnings: list[str] | None = None,
) -> bool:
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem, rdDetermineBonds  # type: ignore
    except Exception:
        return False

    mol = Chem.MolFromPDBBlock(
        pdb_block,
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    if mol is None or mol.GetNumAtoms() == 0:
        return False
    mol.SetProp("_Name", molecule_name)
    total_charge = int(sum(formal_charges))
    try:
        rdDetermineBonds.DetermineBonds(
            mol,
            charge=total_charge,
            allowChargedFragments=True,
            embedChiral=True,
        )
    except Exception as exc:
        if warnings is not None:
            warnings.append(
                f"{molecule_name}: RDKit DetermineBonds falhou ({type(exc).__name__}); usando conectividade PDB/proximidade."
            )
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is None or mol.GetNumAtoms() == 0:
            return False
        mol.SetProp("_Name", molecule_name)
    _apply_rdkit_formal_charges(mol, formal_charges)
    Chem.SanitizeMol(mol)
    Chem.SetAromaticity(mol)
    gasteiger_charges = ""
    try:
        AllChem.ComputeGasteigerCharges(mol)
        gasteiger_charges = _rdkit_gasteiger_charges(mol)
    except Exception as exc:
        if warnings is not None:
            warnings.append(
                f"{molecule_name}: RDKit converteu o ligante, mas nao calculou cargas Gasteiger ({type(exc).__name__})."
            )
    sdf_text = Chem.MolToMolBlock(mol, kekulize=False)
    if not sdf_text or "M  END" not in sdf_text:
        return False
    props = {
        "HIP2L_ChargeModel": "formal PDB charges + RDKit Gasteiger",
        "HIP2L_TotalFormalCharge": str(total_charge),
    }
    if gasteiger_charges:
        props["HIP2L_GasteigerCharges"] = gasteiger_charges
    path.write_text(_normalize_sdf_record(sdf_text, molecule_name, props), encoding="utf-8")
    return True


def _write_sdf_with_external_openbabel(
    path: Path,
    molecule_name: str,
    pdb_block: str,
    formal_charges: list[int],
    chemistry_python: str | Path | None,
    warnings: list[str] | None = None,
) -> bool:
    if not chemistry_python:
        return False
    py_path = Path(chemistry_python)
    if not py_path.exists():
        return False

    from .env_manager import (
        chemistry_process_env,
        external_program_runtime,
        python_prefix,
    )

    prefix = python_prefix(py_path)
    candidates = (
        prefix / "Library" / "bin" / "obabel.exe",
        prefix / "Scripts" / "obabel.exe",
        prefix / "bin" / "obabel",
    )
    obabel = next((candidate for candidate in candidates if candidate.is_file()), None)
    if obabel is None:
        return False

    env = chemistry_process_env(py_path)
    with external_program_runtime():
        mol2_result = subprocess.run(
            [str(obabel), "-ipdb", "-omol2", "--partialcharge", "gasteiger"],
            input=pdb_block,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            cwd=str(path.parent),
        )
    if mol2_result.returncode != 0:
        raise RuntimeError(_external_failure("Open Babel PDB->MOL2", mol2_result))
    mol2_text = mol2_result.stdout
    if "@<TRIPOS>ATOM" not in mol2_text or "@<TRIPOS>BOND" not in mol2_text:
        raise RuntimeError("Open Babel nao retornou um MOL2 valido")

    with external_program_runtime():
        sdf_result = subprocess.run(
            [str(obabel), "-imol2", "-osdf"],
            input=mol2_text,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            cwd=str(path.parent),
        )
    if sdf_result.returncode != 0:
        raise RuntimeError(_external_failure("Open Babel MOL2->SDF", sdf_result))
    if "M  END" not in sdf_result.stdout:
        raise RuntimeError("Open Babel nao retornou um SDF valido")

    charge_text = _mol2_partial_charges(mol2_text)
    props = {
        "HIP2L_ChargeModel": "formal PDB charges + OpenBabel Gasteiger",
        "HIP2L_TotalFormalCharge": str(sum(formal_charges)),
    }
    if charge_text:
        props["HIP2L_GasteigerCharges"] = charge_text
    elif warnings is not None:
        warnings.append(
            f"{molecule_name}: Open Babel converteu o ligante, mas nao retornou cargas Gasteiger."
        )
    path.write_text(
        _normalize_sdf_record(sdf_result.stdout, molecule_name, props),
        encoding="utf-8",
    )
    return True


def _mol2_partial_charges(mol2_text: str) -> str:
    in_atoms = False
    values: list[str] = []
    for line in mol2_text.splitlines():
        marker = line.strip().upper()
        if marker == "@<TRIPOS>ATOM":
            in_atoms = True
            continue
        if marker.startswith("@<TRIPOS>"):
            if in_atoms:
                break
            continue
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            values.append(f"{int(parts[0])}:{float(parts[-1]):.6f}")
        except ValueError:
            continue
    return " ".join(values)


def _external_failure(label: str, result: subprocess.CompletedProcess[str]) -> str:
    detail = (result.stderr.strip() or result.stdout.strip() or "sem saida")
    if result.returncode == -9:
        detail = "processo encerrado pelo sistema (sinal 9; possivel limite de memoria)"
    return f"{label} falhou ({result.returncode}): {detail[:1200]}"


def _write_sdf_with_external_python(
    path: Path,
    molecule_name: str,
    pdb_block: str,
    formal_charges: list[int],
    chemistry_python: str | Path | None,
    warnings: list[str] | None = None,
) -> bool:
    if not chemistry_python:
        return False
    py_path = Path(chemistry_python)
    if not py_path.exists():
        return False
    from .env_manager import chemistry_process_env, external_program_runtime

    payload = {
        "molecule_name": molecule_name,
        "pdb_block": pdb_block,
        "formal_charges": formal_charges,
    }
    env = chemistry_process_env(py_path)
    with external_program_runtime():
        result = subprocess.run(
            [str(py_path), "-I", "-u", "-c", _EXTERNAL_PDB_TO_SDF_SCRIPT],
            input=json.dumps(payload),
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            cwd=str(path.parent),
        )
    if result.returncode != 0:
        raise RuntimeError(_external_failure("Python externo", result))
    try:
        parsed = json.loads(result.stdout.strip().splitlines()[-1])
    except Exception as exc:
        raise RuntimeError(f"Saida invalida do Python externo: {result.stdout[:1200]}") from exc
    if not parsed.get("ok"):
        raise RuntimeError(str(parsed.get("error") or "conversao externa falhou"))
    sdf_text = str(parsed.get("sdf") or "")
    if "M  END" not in sdf_text:
        raise RuntimeError("conversao externa nao retornou SDF valido")
    path.write_text(sdf_text, encoding="utf-8")
    if warnings is not None:
        for message in parsed.get("warnings") or []:
            warnings.append(str(message))
    return True


def _normalize_sdf_record(
    sdf_text: str,
    molecule_name: str,
    properties: dict[str, str] | None = None,
) -> str:
    text = sdf_text.replace("\r\n", "\n").replace("\r", "\n").strip()
    lines = text.splitlines()
    if lines:
        lines[0] = molecule_name
    else:
        lines = [molecule_name, "  HIP2LInterActomics", "", "M  END"]
    text = "\n".join(lines).rstrip()
    if text.endswith("$$$$"):
        text = text[:-4].rstrip()
    if properties:
        for key, value in properties.items():
            text += f"\n>  <{key}>\n{value}\n"
    text += "\n$$$$"
    return text + "\n"


def _apply_rdkit_formal_charges(mol, formal_charges: list[int]) -> None:
    for atom, charge in zip(mol.GetAtoms(), formal_charges):
        if charge:
            atom.SetFormalCharge(int(charge))
    mol.UpdatePropertyCache(strict=False)


def _rdkit_gasteiger_charges(mol) -> str:
    values: list[str] = []
    for atom in mol.GetAtoms():
        try:
            charge = atom.GetProp("_GasteigerCharge")
        except Exception:
            charge = ""
        if not charge or charge.lower() in {"nan", "-nan", "inf", "-inf"}:
            charge = "0.000000"
        values.append(f"{atom.GetIdx() + 1}:{float(charge):.6f}")
    return " ".join(values)


def _openbabel_gasteiger_charges(mol) -> str:
    values: list[str] = []
    for idx, atom in enumerate(mol.atoms, start=1):
        try:
            charge = float(atom.partialcharge)
        except Exception:
            charge = 0.0
        values.append(f"{idx}:{charge:.6f}")
    return " ".join(values)


def _pdb_formal_charge(line: str) -> int:
    token = line[78:80].strip() if len(line) >= 80 else ""
    if not token:
        return 0
    if re.fullmatch(r"\d[+-]", token):
        sign = 1 if token[1] == "+" else -1
        return sign * int(token[0])
    if re.fullmatch(r"[+-]\d", token):
        sign = 1 if token[0] == "+" else -1
        return sign * int(token[1])
    if token in {"+", "-"}:
        return 1 if token == "+" else -1
    return 0


def _write_sdf_fallback(
    path: Path,
    molecule_name: str,
    atom_records: list[str],
    conect_records: list[str],
) -> None:
    atoms = _parse_pdb_ligand_atoms(atom_records)
    if not atoms:
        raise ValueError("nenhum atomo PDB valido para escrever SDF")
    bonds = _pdb_conect_bonds(conect_records, {atom["serial"]: idx for idx, atom in enumerate(atoms, start=1)})
    if not bonds:
        bonds = _infer_ligand_bonds(atoms)
    if len(atoms) > 999 or len(bonds) > 999:
        raise ValueError("SDF V2000 suporta ate 999 atomos e 999 ligacoes neste conversor")

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{molecule_name}\n")
        f.write("  HIP2LInterActomics\n\n")
        f.write(f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0            999 V2000\n")
        for atom in atoms:
            f.write(
                f"{atom['x']:10.4f}{atom['y']:10.4f}{atom['z']:10.4f} "
                f"{atom['element']:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
            )
        for a1, a2, order in bonds:
            f.write(f"{a1:>3}{a2:>3}{order:>3}  0  0  0  0\n")
        f.write("M  END\n$$$$\n")


def _parse_pdb_ligand_atoms(atom_records: list[str]) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for fallback_serial, line in enumerate(atom_records, start=1):
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atom_name = line[12:16].strip() if len(line) >= 16 else ""
        element = line[76:78].strip() if len(line) >= 78 else ""
        element = _normalize_element(element) or _element_from_atom_name(atom_name)
        if element is None:
            continue
        serial = _pdb_atom_serial(line) or fallback_serial
        atoms.append({"serial": serial, "x": x, "y": y, "z": z, "element": element})
    return atoms


def _pdb_conect_bonds(conect_records: list[str], serial_to_index: dict[int, int]) -> list[tuple[int, int, int]]:
    bonds: set[tuple[int, int]] = set()
    for line in conect_records:
        numbers = _parse_conect_numbers(line)
        if len(numbers) < 2:
            continue
        source = serial_to_index.get(numbers[0])
        if source is None:
            continue
        for target_serial in numbers[1:]:
            target = serial_to_index.get(target_serial)
            if target is None or target == source:
                continue
            bonds.add(tuple(sorted((source, target))))
    return [(a1, a2, 1) for a1, a2 in sorted(bonds)]


def _infer_ligand_bonds(atoms: list[dict[str, object]]) -> list[tuple[int, int, int]]:
    bonds: list[tuple[int, int, int]] = []
    for i, atom_a in enumerate(atoms):
        elem_a = str(atom_a["element"])
        radius_a = _COVALENT_RADII.get(elem_a, 0.77)
        for j in range(i + 1, len(atoms)):
            atom_b = atoms[j]
            elem_b = str(atom_b["element"])
            if elem_a == "H" and elem_b == "H":
                continue
            radius_b = _COVALENT_RADII.get(elem_b, 0.77)
            dx = float(atom_a["x"]) - float(atom_b["x"])
            dy = float(atom_a["y"]) - float(atom_b["y"])
            dz = float(atom_a["z"]) - float(atom_b["z"])
            dist2 = dx * dx + dy * dy + dz * dz
            max_dist = radius_a + radius_b + 0.45
            if 0.16 <= dist2 <= max_dist * max_dist:
                bonds.append((i + 1, j + 1, 1))
    return bonds


def _write_mol2(
    path: Path, name: str,
    atoms: list[tuple[int, str]], bonds: list[tuple[int, int, int, str]],
) -> None:
    na = len(atoms)
    nb = len(bonds)
    eoia = len(str(na))
    eola = len(str(nb)) if nb else 1

    with path.open("w", encoding="utf-8") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{name}_LIG\n")
        f.write(f"{na} {nb} 1 0 0\n")
        f.write("SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for new_idx, rest in atoms:
            f.write(f"  {str(new_idx).rjust(eoia)}{rest}")
        f.write("@<TRIPOS>BOND\n")
        for b_idx, a1, a2, btype in bonds:
            f.write(
                f"  {str(b_idx).rjust(eola)} {str(a1).rjust(eoia)} "
                f"{str(a2).rjust(eoia)} {btype:>4}\n"
            )


def _format_pdb_atom_line(
    record_name: str,
    idx: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    icode: str,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    atom_field = _format_pdb_atom_name(atom_name, element)
    record = str(record_name or "ATOM").strip().upper()[:6] or "ATOM"
    return (
        f"{record:<6}{idx:5d} {atom_field} {resname:>3} {chain_id}{resseq:4d}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}\n"
    )


def _format_pdb_atom_name(atom_name: str, element: str) -> str:
    atom_name = re.sub(r"\s+", "", atom_name)[:4]
    if len(element.strip()) == 1 and len(atom_name) < 4:
        return atom_name.rjust(4)
    return atom_name.ljust(4)


def _infer_element(atom_name: str, atom_type: str) -> str | None:
    for candidate in (_element_from_atom_type(atom_type), _element_from_atom_name(atom_name)):
        if candidate:
            return candidate
    return None


def _element_from_atom_type(atom_type: str) -> str | None:
    token = atom_type.split(".", 1)[0].strip()
    if not token:
        return None
    return _normalize_element(token)


def _element_from_atom_name(atom_name: str) -> str | None:
    cleaned = re.sub(r"[^A-Za-z]", "", atom_name or "")
    if not cleaned:
        return None
    one = _normalize_element(cleaned[:1])
    two = _normalize_element(cleaned[:2]) if len(cleaned) >= 2 else None
    if len(cleaned) <= 2 and two:
        return two
    return one or two


def _normalize_element(token: str) -> str | None:
    letters = re.sub(r"[^A-Za-z]", "", token or "")
    if not letters:
        return None
    if len(letters) >= 2:
        candidate = letters[0].upper() + letters[1].lower()
        if candidate in _VALID_ELEMENTS:
            return candidate
    candidate = letters[0].upper()
    if candidate in _VALID_ELEMENTS:
        return candidate
    return None
