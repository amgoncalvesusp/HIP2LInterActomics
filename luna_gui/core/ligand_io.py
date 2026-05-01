"""Lightweight ligand-file parsing — extracts molecule names without RDKit.

Supported formats: MOL2, SDF/MOL. The molecule name is what LUNA's
`MolFileEntry.from_file` expects, one per line, in entries.txt.
"""
from __future__ import annotations

import os
import tempfile
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


def consolidate_folder(
    folder: str | Path,
    output_mol2: str | Path,
    use_file_stem_as_name: bool = False,
) -> int:
    """Concatenate all MOL2/SDF files in `folder` into a single MOL2.

    SDF molecules are not converted (would require Open Babel/RDKit) — only
    MOL2 files in the folder are merged. Returns number of molecules written.

    For mixed folders, the user should pre-convert via Open Babel.
    """
    folder = Path(folder)
    out = Path(output_mol2)
    inputs = _iter_input_mol2_files(folder, out)
    n = 0
    with _atomic_writer(out) as fout:
        for f in inputs:
            blocks = _split_mol2_molecules(f)
            for idx, (block, _name) in enumerate(blocks, start=1):
                if use_file_stem_as_name:
                    content = _replace_mol2_name(
                        block, _name_from_source_file(f, idx, len(blocks))
                    )
                else:
                    content = "".join(block)
                fout.write(content)
                if not content.endswith("\n"):
                    fout.write("\n")
                n += 1
    return n


def consolidate_sdf_folder(
    folder: str | Path,
    output_sdf: str | Path,
    use_file_stem_as_name: bool = False,
) -> tuple[int, list[str]]:
    """Concatenate all SDF/MOL files in `folder` into a single SDF file."""
    folder = Path(folder)
    out = Path(output_sdf)
    inputs = _iter_input_sdf_files(folder, out)
    total = 0
    names: list[str] = []
    with _atomic_writer(out) as fout:
        for f in inputs:
            blocks = _split_sdf_molecules(f)
            for idx, block in enumerate(blocks, start=1):
                name = _sdf_block_name(block)
                if use_file_stem_as_name:
                    name = _name_from_source_file(f, idx, len(blocks))
                    block = _replace_sdf_name(block, name)
                content = "".join(block).rstrip()
                if not content:
                    continue
                fout.write(content)
                fout.write("\n$$$$\n")
                total += 1
                names.append(name)
    return total, names


def consolidate_ligand_folder(
    folder: str | Path,
    output_dir: str | Path | None = None,
    use_file_stem_as_name: bool = True,
) -> tuple[int, list[str], Path]:
    """Consolidate a folder of MOL2 or SDF/MOL ligands into one ligand file.

    A single output format is produced because LUNA receives one ligand file.
    Mixed MOL2 + SDF folders are rejected so files are not silently ignored.
    """
    folder = Path(folder)
    output_dir = Path(output_dir) if output_dir is not None else folder
    mol2_files = _iter_input_mol2_files(folder, output_dir / "_consolidated_ligands.mol2")
    sdf_files = _iter_input_sdf_files(folder, output_dir / "_consolidated_ligands.sdf")
    if mol2_files and sdf_files:
        raise ValueError(
            "A pasta contém arquivos MOL2 e SDF/MOL. Use uma pasta com apenas um "
            "formato por vez para gerar um arquivo consolidado compatível com o LUNA."
        )
    if mol2_files:
        out = output_dir / "_consolidated_ligands.mol2"
        return (*consolidate_folder_clean(
            folder,
            out,
            drop_lp=True,
            use_file_stem_as_name=use_file_stem_as_name,
        ), out)
    if sdf_files:
        out = output_dir / "_consolidated_ligands.sdf"
        n, names = consolidate_sdf_folder(
            folder,
            out,
            use_file_stem_as_name=use_file_stem_as_name,
        )
        return n, names, out
    return 0, [], output_dir / "_consolidated_ligands.mol2"


def consolidate_folder_clean(
    folder: str | Path,
    output_mol2: str | Path,
    drop_lp: bool = True,
    use_file_stem_as_name: bool = False,
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
    inputs = _iter_input_mol2_files(folder, out)
    total = 0
    names: list[str] = []
    with _atomic_writer(out) as fout:
        for f in inputs:
            blocks = _split_mol2_molecules(f)
            for idx, (block, name) in enumerate(blocks, start=1):
                if use_file_stem_as_name:
                    name = _name_from_source_file(f, idx, len(blocks))
                cleaned = _clean_mol2_block(
                    block,
                    drop_lp=drop_lp,
                    molecule_name=name,
                )
                if cleaned is None:
                    continue
                fout.write(cleaned)
                if not cleaned.endswith("\n"):
                    fout.write("\n")
                total += 1
                names.append(name)
    return total, names


def strip_hydrogens_from_mol2_file(
    input_mol2: str | Path,
    output_mol2: str | Path,
    drop_lp: bool = True,
) -> tuple[int, list[str]]:
    """Write a hydrogen-free MOL2 copy while preserving molecule names/order."""
    src = Path(input_mol2)
    out = Path(output_mol2)
    total = 0
    names: list[str] = []
    with _atomic_writer(out) as fout:
        for block, name in _split_mol2_molecules(src):
            cleaned = _clean_mol2_block(
                block,
                drop_lp=drop_lp,
                drop_hydrogens=True,
                molecule_name=name,
            )
            if cleaned is None:
                continue
            fout.write(cleaned)
            if not cleaned.endswith("\n"):
                fout.write("\n")
            total += 1
            names.append(name)
    return total, names


def _iter_input_mol2_files(folder: Path, output_mol2: Path) -> list[Path]:
    """Return the MOL2 inputs, excluding the output file itself."""
    try:
        out_resolved = output_mol2.resolve(strict=False)
    except Exception:
        out_resolved = output_mol2

    inputs: list[Path] = []
    for candidate in sorted(folder.iterdir()):
        if candidate.suffix.lower() != ".mol2":
            continue
        try:
            same_file = candidate.resolve(strict=False) == out_resolved
        except Exception:
            same_file = candidate == output_mol2
        if same_file:
            continue
        inputs.append(candidate)
    return inputs


def _iter_input_sdf_files(folder: Path, output_sdf: Path) -> list[Path]:
    """Return SDF/MOL inputs, excluding the output file itself."""
    try:
        out_resolved = output_sdf.resolve(strict=False)
    except Exception:
        out_resolved = output_sdf

    inputs: list[Path] = []
    for candidate in sorted(folder.iterdir()):
        if candidate.suffix.lower() not in {".sdf", ".mol"}:
            continue
        try:
            same_file = candidate.resolve(strict=False) == out_resolved
        except Exception:
            same_file = candidate == output_sdf
        if same_file:
            continue
        inputs.append(candidate)
    return inputs


class _AtomicWriter:
    def __init__(self, final_path: Path) -> None:
        self.final_path = final_path
        self.tmp_path: Path | None = None
        self.handle = None

    def __enter__(self):
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(
            prefix=self.final_path.stem + "_",
            suffix=self.final_path.suffix,
            dir=str(self.final_path.parent),
        )
        os.close(fd)
        self.tmp_path = Path(tmp_name)
        self.handle = self.tmp_path.open("w", encoding="utf-8")
        return self.handle

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.handle is not None:
            self.handle.close()
        if self.tmp_path is None:
            return
        if exc_type is None:
            os.replace(self.tmp_path, self.final_path)
        elif self.tmp_path.exists():
            self.tmp_path.unlink()


def _atomic_writer(final_path: Path) -> _AtomicWriter:
    return _AtomicWriter(final_path)


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


def _split_sdf_molecules(path: Path) -> list[list[str]]:
    """Split an SDF/MOL file into molecule blocks without the $$$$ marker."""
    blocks: list[list[str]] = []
    cur: list[str] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                if cur and any(part.strip() for part in cur):
                    blocks.append(cur)
                cur = []
            else:
                cur.append(line)
    if cur and any(part.strip() for part in cur):
        blocks.append(cur)
    return blocks


def _sdf_block_name(block: list[str]) -> str:
    return block[0].strip() if block else ""


def _replace_sdf_name(block: list[str], molecule_name: str) -> list[str]:
    updated = list(block)
    if updated:
        updated[0] = molecule_name + "\n"
    else:
        updated = [molecule_name + "\n"]
    return updated


def _name_from_source_file(path: Path, block_idx: int, total_blocks: int) -> str:
    """Return the filename stem or a stable per-block variant."""
    if total_blocks <= 1:
        return path.stem
    return f"{path.stem}__{block_idx}"


def _replace_mol2_name(lines: list[str], molecule_name: str) -> str:
    """Replace the molecule-name line in a raw MOL2 block."""
    block = list(lines)
    if len(block) >= 2 and block[0].startswith("@<TRIPOS>MOLECULE"):
        block[1] = molecule_name + "\n"
    return "".join(block)


def _clean_mol2_block(
    lines: list[str],
    drop_lp: bool,
    drop_hydrogens: bool = False,
    molecule_name: str | None = None,
) -> str | None:
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
        if drop_hydrogens and _is_hydrogen_atom(parts):
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
    if molecule_name and len(header) >= 2:
        header[1] = molecule_name + "\n"
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


def _is_hydrogen_atom(parts: list[str]) -> bool:
    atom_name = parts[1].strip().upper() if len(parts) > 1 else ""
    atom_type = parts[5].strip().upper() if len(parts) > 5 else ""
    trimmed_name = atom_name.lstrip("0123456789")
    return atom_type.startswith("H") or trimmed_name.startswith("H")


def _next_section(lines: list[str], after: int) -> int:
    for j in range(after + 1, len(lines)):
        if lines[j].startswith("@<TRIPOS>"):
            return j
    return len(lines)
