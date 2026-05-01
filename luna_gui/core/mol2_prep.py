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

_WATER_RESIDUES: frozenset[str] = frozenset({
    "HOH", "WAT", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD",
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

    return PrepResult(
        files_processed=nfiles, proteins_written=nprot, ligands_written=nlig,
        water_molecules_detected=nwater,
        protein_dir=str(prot_dir), ligand_dir=str(lig_dir), errors=errors,
    )


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
                    chain_id = chain_by_subst.get(subst_id, "A")
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    element = _infer_element(atom_name, atom_type)
                    if element is None:
                        continue
                    pdb_line = _format_pdb_atom_line(
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
    return (
        f"ATOM  {idx:5d} {atom_field} {resname:>3} {chain_id}{resseq:4d}{icode}   "
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
