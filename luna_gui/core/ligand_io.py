"""Lightweight ligand-file parsing — extracts molecule names without RDKit.

Supported formats: MOL2, SDF/MOL. The molecule name is what LUNA's
`MolFileEntry.from_file` expects, one per line, in entries.txt.
"""
from __future__ import annotations

from pathlib import Path


def parse_ligand_file(path: str | Path) -> list[str]:
    """Return the list of molecule names found in `path`.

    Order is preserved. Empty / unnamed molecules are kept as ''.
    """
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".mol2":
        return _parse_mol2(p)
    if ext in (".sdf", ".mol"):
        return _parse_sdf(p)
    raise ValueError(f"Unsupported ligand file format: {ext}")


def _parse_mol2(path: Path) -> list[str]:
    names: list[str] = []
    with path.open("r", errors="replace") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("@<TRIPOS>MOLECULE"):
            # Next non-empty line is the molecule name.
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            names.append(lines[j].strip() if j < len(lines) else "")
            i = j + 1
        else:
            i += 1
    return names


def _parse_sdf(path: Path) -> list[str]:
    names: list[str] = []
    with path.open("r", errors="replace") as f:
        block: list[str] = []
        for line in f:
            if line.startswith("$$$$"):
                if block:
                    names.append(block[0].strip())
                block = []
            else:
                block.append(line)
        if block and any(b.strip() for b in block):
            names.append(block[0].strip())
    return names


def write_entries_file(path: str | Path, names: list[str]) -> None:
    Path(path).write_text("\n".join(names) + "\n", encoding="utf-8")


def consolidate_folder(folder: str | Path, output_mol2: str | Path) -> int:
    """Concatenate all MOL2/SDF files in `folder` into a single MOL2.

    SDF molecules are not converted (would require Open Babel/RDKit) — only
    MOL2 files in the folder are merged. Returns number of molecules written.

    For mixed folders, the user should pre-convert via Open Babel.
    """
    folder = Path(folder)
    out = Path(output_mol2)
    n = 0
    with out.open("w", encoding="utf-8") as fout:
        for f in sorted(folder.iterdir()):
            if f.suffix.lower() != ".mol2":
                continue
            content = f.read_text(errors="replace")
            fout.write(content)
            if not content.endswith("\n"):
                fout.write("\n")
            n += content.count("@<TRIPOS>MOLECULE")
    return n


def consolidate_folder_clean(
    folder: str | Path,
    output_mol2: str | Path,
    drop_lp: bool = True,
) -> tuple[int, list[str]]:
    """Robust consolidation: per-file clean + renumber atoms/bonds from 1.

    Inspired by the Daniel preprocessing notebook. For each .mol2 in `folder`:
      - drops atoms whose name is 'LP' (lone pairs) when drop_lp=True
      - renumbers atoms 1..N (keeps names, coords, types, charges)
      - renumbers bonds 1..M referencing the new atom IDs
      - rewrites the MOLECULE counts line

    Returns (number_of_molecules, list_of_molecule_names).
    """
    folder = Path(folder)
    out = Path(output_mol2)
    total = 0
    names: list[str] = []
    with out.open("w", encoding="utf-8") as fout:
        for f in sorted(folder.iterdir()):
            if f.suffix.lower() != ".mol2":
                continue
            for block, name in _split_mol2_molecules(f):
                cleaned = _clean_mol2_block(block, drop_lp=drop_lp)
                if cleaned is None:
                    continue
                fout.write(cleaned)
                if not cleaned.endswith("\n"):
                    fout.write("\n")
                total += 1
                names.append(name)
    return total, names


def _split_mol2_molecules(path: Path) -> list[tuple[list[str], str]]:
    """Split a multi-mol2 into (lines, name) chunks."""
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    out: list[tuple[list[str], str]] = []
    cur: list[str] = []
    name = ""
    for line in lines:
        if line.startswith("@<TRIPOS>MOLECULE"):
            if cur:
                out.append((cur, name))
            cur = [line]
            name = ""
        else:
            cur.append(line)
            if name == "" and len(cur) == 2:
                # line right after the MOLECULE header is the name
                name = line.strip()
    if cur:
        out.append((cur, name))
    return out


def _clean_mol2_block(lines: list[str], drop_lp: bool) -> str | None:
    """Return a cleaned MOL2 block or None if it has no atoms."""
    # Locate section indices
    sections: dict[str, int] = {}
    for i, l in enumerate(lines):
        if l.startswith("@<TRIPOS>"):
            sections[l.strip()] = i
    if "@<TRIPOS>MOLECULE" not in sections or "@<TRIPOS>ATOM" not in sections:
        return None

    mol_start = sections["@<TRIPOS>MOLECULE"]
    atom_start = sections["@<TRIPOS>ATOM"]
    bond_start = sections.get("@<TRIPOS>BOND")
    # Cap for ATOM/BOND sections = next section after it
    atom_end = bond_start if bond_start is not None else _next_section(lines, atom_start)
    bond_end = _next_section(lines, bond_start) if bond_start is not None else None

    # Parse atoms, remap IDs 1..N, drop LP
    atom_lines = lines[atom_start + 1:atom_end]
    new_atoms: list[str] = []
    id_map: dict[int, int] = {}
    new_idx = 0
    for al in atom_lines:
        parts = al.split()
        if len(parts) < 6:
            continue
        try:
            old_id = int(parts[0])
        except ValueError:
            continue
        atom_name = parts[1]
        if drop_lp and atom_name.upper() == "LP":
            continue
        new_idx += 1
        id_map[old_id] = new_idx
        parts[0] = str(new_idx)
        new_atoms.append(" ".join(parts) + "\n")
    if not new_atoms:
        return None

    # Parse bonds with remapped atom IDs
    new_bonds: list[str] = []
    if bond_start is not None:
        bond_lines = lines[bond_start + 1:bond_end]
        bidx = 0
        for bl in bond_lines:
            parts = bl.split()
            if len(parts) < 4:
                continue
            try:
                a1, a2 = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if a1 not in id_map or a2 not in id_map:
                continue  # skip bonds that reference dropped atoms
            bidx += 1
            new_bonds.append(f"  {bidx:>5} {id_map[a1]:>5} {id_map[a2]:>5} {parts[3]:>4}\n")

    # Rewrite MOLECULE counts line (line mol_start + 2): "<n_atoms> <n_bonds> <n_subst> ..."
    # Header format: [MOLECULE tag, name, counts, mol_type, charge_type, ...]
    header = lines[mol_start:atom_start]
    # Locate the counts line (first line after the name with 3+ integers)
    counts_idx = None
    for i in range(mol_start + 2, atom_start):
        toks = lines[i].split()
        if len(toks) >= 3 and all(t.lstrip("-").isdigit() for t in toks[:3]):
            counts_idx = i - mol_start
            break
    if counts_idx is not None:
        toks = header[counts_idx].split()
        toks[0] = str(len(new_atoms))
        toks[1] = str(len(new_bonds))
        header[counts_idx] = " ".join(toks) + "\n"

    # Reassemble block
    out: list[str] = list(header)
    out.append("@<TRIPOS>ATOM\n")
    out.extend(new_atoms)
    if new_bonds:
        out.append("@<TRIPOS>BOND\n")
        out.extend(new_bonds)
    # Preserve any trailing sections (SUBSTRUCTURE etc.) untouched
    tail_start = None
    for tag, idx in sections.items():
        if tag not in ("@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND"):
            if tail_start is None or idx < tail_start:
                tail_start = idx
    if tail_start is not None:
        out.extend(lines[tail_start:])
    return "".join(out)


def _next_section(lines: list[str], after: int) -> int:
    for j in range(after + 1, len(lines)):
        if lines[j].startswith("@<TRIPOS>"):
            return j
    return len(lines)
