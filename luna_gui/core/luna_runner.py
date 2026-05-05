"""Build the LUNA `run.py` command line from a ProjectConfig."""
from __future__ import annotations

from pathlib import Path

from .project import ProjectConfig


def build_command(
    py_exe: str,
    run_py: str,
    cfg: ProjectConfig,
    entries_file: str,
) -> list[str]:
    """Assemble argv for `python run.py ...` based on the project config.

    Targets the modern LUNA CLI (luna>=0.13):
      -p PDB_ID, --pdbdir DIR, -l MOL_FILE, -e ENTRIES, -w WORKDIR
      --ifp, --ifp-out, --ifp-matrix, --pse, --psedir, --bind, --fork
    """
    # LUNA now expects a PDB id + a directory instead of a full file path.
    protein_path = Path(cfg.protein_file)
    pdb_id = protein_path.stem
    pdb_dir = str(protein_path.parent)

    cmd: list[str] = [
        py_exe, run_py,
        "-p", pdb_id,
        "--pdbdir", pdb_dir,
        "-l", cfg.ligand_file,
        "-e", entries_file,
        "-w", cfg.workdir,
    ]

    if cfg.out_ifp:
        cmd += [
            "--ifp",
            "-T", cfg.ifp_type,
            "-L", str(cfg.ifp_levels),
            "-R", str(cfg.ifp_radius),
            "-S", str(cfg.ifp_length),
        ]
        if cfg.ifp_bit:
            cmd.append("-B")
        if cfg.ifp_output:
            cmd += ["--ifp-out", cfg.ifp_output]

    if cfg.sim_matrix:
        if not cfg.out_ifp:
            cmd.append("--ifp")
        out = cfg.sim_matrix_output or str(Path(cfg.workdir) / "sim_matrix.csv")
        cmd += ["--ifp-matrix", out]

    if cfg.out_pse:
        cmd.append("--pse")
        if cfg.pse_path:
            cmd += ["--psedir", cfg.pse_path]

    if cfg.filter_binding_modes and cfg.binding_modes_cfg:
        cmd += ["--bind", cfg.binding_modes_cfg]

    if cfg.fork_from:
        cmd += ["--fork", cfg.fork_from]

    if cfg.overwrite:
        cmd.append("--overwrite")

    # Always pass nproc explicitly — LUNA's default is -1 (all cores),
    # which breaks on Windows due to a multiprocessing pickling bug.
    cmd += ["--nproc", str(max(1, cfg.nproc or 1))]

    return cmd


def validate(cfg: ProjectConfig) -> list[str]:
    """Return a list of human-readable errors. Empty list = OK."""
    errs: list[str] = []
    protein_path = Path(cfg.protein_file) if cfg.protein_file else None
    if not cfg.protein_file or protein_path is None or not protein_path.exists():
        errs.append("Arquivo de proteína não encontrado.")
    elif protein_path.is_dir() and not any(protein_path.glob("*.pdb")):
        errs.append("A pasta de proteínas não contém arquivos .pdb.")
    ligand_path = Path(cfg.ligand_file) if cfg.ligand_file else None
    if not cfg.ligand_file or ligand_path is None or not ligand_path.exists():
        errs.append("Arquivo de ligantes não encontrado.")
    elif ligand_path.is_dir() and not any(
        item.is_file() and item.suffix.lower() in {".mol2", ".sdf", ".sd", ".mol", ".pdb", ".ent"}
        for item in ligand_path.iterdir()
    ):
        errs.append("A pasta de ligantes não contém arquivos .mol2, .sdf, .mol, .pdb ou .ent.")
    if not cfg.selected_ligands:
        errs.append("Nenhum ligante selecionado.")
    if not cfg.workdir:
        errs.append("Diretório de trabalho não definido.")
    if cfg.filter_binding_modes and not Path(cfg.binding_modes_cfg or "").exists():
        errs.append("Arquivo de binding modes (.cfg) não encontrado.")
    if cfg.fork_from and not Path(cfg.fork_from).exists():
        errs.append("Projeto de fork não encontrado.")
    return errs
