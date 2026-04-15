"""Split docking MOL2 complexes (protein + ligand) into separate files.

Adapted from `Scripts_Daniel/0-reorganizar_moleculas.py`. Given a folder of
MOL2 files where each file contains a protein (atoms 1..last_pa) followed by a
ligand (atoms > last_pa), produces:
  <out_dir>/proteinas_pdb/<name>.pdb
  <out_dir>/ligantes_mol2/<name>_ligand.mol2
"""
from __future__ import annotations

import re
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
    "HSD", "HSE", "HSP", "HIE", "HID", "HIP", "CYX", "CYM",
    # N/C terminus caps
    "ACE", "NME", "NHE",
    # Water / solvent
    "HOH", "WAT", "TIP", "SOL", "T3P",
    # Selenium / non-standard but common
    "MSE", "SEC",
})


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
        prefix = re.sub(r'[^A-Za-z]', '', name).upper()[:3]
        if prefix not in _PROTEIN_RESIDUES:
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
    # Count atoms
    total = _count_atoms(lines)
    return DetectionResult(
        last_pa=last_pa, method="substructure",
        n_protein_atoms=last_pa,
        n_ligand_atoms=total - last_pa,
        ligand_names=list(dict.fromkeys(lig_names)),  # deduplicated, ordered
    )


def _detect_via_atom_scan(lines: list[str]) -> DetectionResult:
    """Walk ATOM records; last atom whose subst_name prefix is in protein set."""
    in_atom = False
    last_prot = 0
    total = 0
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
        subst_name = parts[7] if len(parts) > 7 else ""
        prefix = re.sub(r'[^A-Za-z]', '', subst_name).upper()[:3]
        if prefix in _PROTEIN_RESIDUES:
            last_prot = max(last_prot, idx)
        else:
            if subst_name and subst_name not in seen_lig:
                seen_lig.add(subst_name)
                lig_names.append(subst_name)

    last_pa = last_prot if last_prot > 0 else max(1, total - 1)
    return DetectionResult(
        last_pa=last_pa, method="atom_scan" if last_prot > 0 else "fallback",
        n_protein_atoms=last_pa,
        n_ligand_atoms=total - last_pa,
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


@dataclass
class PrepResult:
    files_processed: int
    proteins_written: int
    ligands_written: int
    protein_dir: str
    ligand_dir: str
    errors: list[str]


def split_docking_folder(
    src_folder: str | Path,
    last_pa: int,
    out_folder: str | Path | None = None,
    protein_subdir: str = "proteinas_pdb",
    ligand_subdir: str = "ligantes_mol2",
) -> PrepResult:
    """Split every .mol2 in `src_folder` into PDB (protein) + MOL2 (ligand)."""
    src = Path(src_folder)
    out = Path(out_folder) if out_folder else src
    prot_dir = out / protein_subdir
    lig_dir = out / ligand_subdir
    prot_dir.mkdir(parents=True, exist_ok=True)
    lig_dir.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    nfiles = nprot = nlig = 0
    first_lig_id = last_pa + 1

    for mol2 in sorted(src.glob("*.mol2")):
        nfiles += 1
        try:
            with mol2.open("r", errors="replace") as fh:
                lines = fh.readlines()
            prot_atoms, lig_atoms, lig_bonds = _split_one(lines, last_pa, first_lig_id)
        except Exception as e:
            errors.append(f"{mol2.name}: {e}")
            continue

        name = mol2.stem
        if prot_atoms:
            _write_pdb(prot_dir / f"{name}.pdb", prot_atoms)
            nprot += 1
        if lig_atoms:
            _write_mol2(lig_dir / f"{name}_ligand.mol2", name, lig_atoms, lig_bonds)
            nlig += 1

    return PrepResult(
        files_processed=nfiles, proteins_written=nprot, ligands_written=nlig,
        protein_dir=str(prot_dir), ligand_dir=str(lig_dir), errors=errors,
    )


def _split_one(
    lines: list[str], last_pa: int, first_lig_id: int,
) -> tuple[list[str], list[tuple[int, str]], list[tuple[int, int, int, str]]]:
    """Split raw MOL2 lines into protein PDB records + ligand atoms/bonds."""
    section: str | None = None
    prot_atoms: list[str] = []
    lig_atoms: list[tuple[int, str]] = []
    lig_bonds: list[tuple[int, int, int, str]] = []

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

            if idx <= last_pa:
                # Protein atom → PDB format
                try:
                    atom_name = parts[1]
                    subst_name = parts[7] if len(parts) > 7 else "UNK1"
                    resname, resseq, icode = _parse_subst_name(subst_name)
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    # Standard PDB fixed-column format (LUNA parser is strict):
                    # cols  1- 6: record type
                    # cols  7-11: atom serial
                    # col  12:    space
                    # cols 13-16: atom name (left-justified)
                    # col  17:    alt loc (space)
                    # cols 18-20: residue name
                    # col  21:    space
                    # col  22:    chain ID
                    # cols 23-26: residue seq number (INTEGER — must be numeric)
                    # col  27:    insertion code
                    # cols 28-30: spaces
                    # cols 31-38: X   cols 39-46: Y   cols 47-54: Z
                    # cols 55-60: occupancy   cols 61-66: B-factor
                    pdb_line = (
                        f"ATOM  {idx:5d} {atom_name:<4} {resname:<3} A{resseq:4d}{icode}   "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
                    )
                    prot_atoms.append(pdb_line)
                except (ValueError, IndexError):
                    continue
            else:
                # Ligand atom — keep the line intact after remapping the ID
                new_idx = idx - last_pa
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
            if a1 >= first_lig_id and a2 >= first_lig_id:
                lig_bonds.append((
                    len(lig_bonds) + 1,
                    a1 - last_pa,
                    a2 - last_pa,
                    parts[3],
                ))

    return prot_atoms, lig_atoms, lig_bonds


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
