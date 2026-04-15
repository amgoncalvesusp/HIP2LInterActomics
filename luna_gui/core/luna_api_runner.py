"""Run LUNA via its Python API (inside luna-env) instead of the CLI.

Used when the user enables any advanced option (water handling, custom
InteractionCalculator flags, DefaultInteractionConfig overrides, per-type
PSE filtering). Writes a self-contained Python script to the workdir and
returns the argv needed to run it.
"""
from __future__ import annotations

import json
from pathlib import Path

from .project import ProjectConfig


API_RUNNER_SCRIPT = r'''
"""Auto-generated LUNA runner. Reads params from argv[1] (JSON file)."""
import sys, os, json, csv
from pathlib import Path

params_file = sys.argv[1]
with open(params_file, "r", encoding="utf-8") as fh:
    p = json.load(fh)

import luna
from luna.mol.entry import MolFileEntry
from luna.projects import LocalProject
from luna.interaction.config import DefaultInteractionConfig
from luna.interaction.filter import InteractionFilter
from luna.interaction.calc import InteractionCalculator
from luna.interaction.fp.type import IFPType

workdir   = p["workdir"]
pdb_id    = p["pdb_id"]
pdb_dir   = p["pdb_dir"]
lig_file  = p["lig_file"]
entries   = p["entries"]        # list of molecule names

# ----- Build MolFileEntry objects -----
entry_objs = []
for name in entries:
    try:
        e = MolFileEntry.from_mol_file(
            pdb_id=pdb_id, mol_id=name, mol_file=lig_file,
            is_multimol_file=True, mol_obj_type="openbabel",
        )
        entry_objs.append(e)
    except Exception as ex:
        print(f"[warn] pulando '{name}': {ex}", flush=True)
print(f"[luna-api] {len(entry_objs)} entries carregadas", flush=True)

# ----- DefaultInteractionConfig with overrides -----
inter_config = DefaultInteractionConfig()
for k, v in (p.get("inter_config_overrides") or {}).items():
    try:
        inter_config[k] = v
        print(f"[luna-api] inter_config[{k}] = {v}", flush=True)
    except Exception as ex:
        print(f"[warn] inter_config[{k}] inválido: {ex}", flush=True)

# ----- InteractionFilter + Calculator -----
inter_filter = InteractionFilter.new_pli_filter(
    ignore_self_inter=p.get("ic_ignore_self_inter", True),
    ignore_any_h2o=not p.get("include_waters", False),
)
# Build IC kwargs based on LUNA version compatibility
ic_kwargs = {
    "inter_filter": inter_filter,
    "inter_config": inter_config,
    "add_proximal": p.get("ic_add_proximal", False),
    "add_atom_atom": p.get("ic_add_atom_atom", False),
    "add_dependent_inter": p.get("ic_add_dependent_inter", True),
}
# Try to add h2o flag if supported (newer LUNA versions)
try:
    # Test if the flag is supported by checking the signature
    import inspect
    ic_sig = inspect.signature(InteractionCalculator.__init__)
    if "add_h2o_pairs_with_no_target" in ic_sig.parameters:
        ic_kwargs["add_h2o_pairs_with_no_target"] = p.get("ic_add_h2o_pairs_with_no_target", True)
except Exception:
    pass
inter_calc = InteractionCalculator(**ic_kwargs)

# ----- LocalProject -----
opts = dict(
    entries=entry_objs,
    working_path=workdir,
    pdb_path=pdb_dir,
    overwrite_path=p.get("overwrite", False),
    add_h=True, ph=7.4, amend_mol=True,
    calc_ifp=p.get("out_ifp", True),
    out_pse=p.get("out_pse", False),
    nproc=p.get("nproc", 1),
    inter_calc=inter_calc,
)
if p.get("out_ifp", True):
    opts["ifp_type"] = getattr(IFPType, p.get("ifp_type", "EIFP"))
    opts["ifp_num_levels"] = p.get("ifp_levels", 2)
    opts["ifp_radius_step"] = p.get("ifp_radius", 5.73171)
    opts["ifp_length"] = p.get("ifp_length", 4096)
    opts["ifp_count"] = not p.get("ifp_bit", False)

proj = LocalProject(**opts)
print("[luna-api] iniciando proj.run()...", flush=True)
proj.run()

# ----- Export IFP CSV (matches the CLI's output layout) -----
if p.get("out_ifp", True):
    ifp_out = p.get("ifp_output") or os.path.join(workdir, "results", "fingerprints", "ifp.csv")
    os.makedirs(os.path.dirname(ifp_out), exist_ok=True)
    with open(ifp_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand_id", "on_bits"])
        for entry in proj.entries:
            er = proj.get_entry_results(entry)
            fp = getattr(er, "ifp", None)
            if fp is None:
                continue
            on_bits = fp.get_on_bits() if hasattr(fp, "get_on_bits") else list(fp.indices)
            w.writerow([entry.to_string(), "\t".join(str(b) for b in on_bits)])
    print(f"[luna-api] IFP CSV salvo em {ifp_out}", flush=True)

# ----- Similarity matrix -----
if p.get("sim_matrix", False):
    sm_out = p.get("sim_matrix_output") or os.path.join(workdir, "sim_matrix.csv")
    ifps = list(proj.ifps)  # list of (entry, fp)
    labels = [e.to_string() for e, _ in ifps]
    with open(sm_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([""] + labels)
        for i, (_ei, fi) in enumerate(ifps):
            row = [labels[i]]
            for j, (_ej, fj) in enumerate(ifps):
                try:
                    row.append(f"{fi.calc_similarity(fj):.4f}")
                except Exception:
                    row.append("")
            w.writerow(row)
    print(f"[luna-api] Sim matrix salva em {sm_out}", flush=True)

# ----- PSE sessions (with optional per-type filter) -----
if p.get("out_pse", False):
    try:
        from luna.interaction.view import InteractionViewer
    except Exception as ex:
        print(f"[warn] InteractionViewer indisponível: {ex}", flush=True)
    else:
        pse_dir = p.get("pse_path") or os.path.join(workdir, "results", "pse")
        os.makedirs(pse_dir, exist_ok=True)
        types_filter = p.get("pse_interaction_types") or []
        for entry in proj.entries:
            er = proj.get_entry_results(entry)
            if types_filter:
                inter = list(er.interactions_mngr.filter_by_types(types_filter))
            else:
                inter = er.interactions_mngr
            safe_name = entry.to_string().replace(":", "_").replace("/", "_")
            pse_path = os.path.join(pse_dir, f"{safe_name}.pse")
            try:
                viewer = InteractionViewer(show_hydrop_surface=False)
                viewer.new_session([(entry, inter, proj.pdb_path)], pse_path)
            except Exception as ex:
                print(f"[warn] PSE para {safe_name} falhou: {ex}", flush=True)
        print(f"[luna-api] PSE salvos em {pse_dir}", flush=True)

print("[luna-api] concluído", flush=True)
'''


def write_runner(workdir: str | Path) -> str:
    """Write the API runner script to the workdir and return its path."""
    wd = Path(workdir); wd.mkdir(parents=True, exist_ok=True)
    path = wd / "_luna_api_runner.py"
    path.write_text(API_RUNNER_SCRIPT, encoding="utf-8")
    return str(path)


def write_params(workdir: str | Path, cfg: ProjectConfig, entries: list[str]) -> str:
    """Dump the relevant ProjectConfig subset as JSON next to the runner."""
    wd = Path(workdir)
    protein_path = Path(cfg.protein_file)
    params = {
        "workdir": cfg.workdir,
        "pdb_id": protein_path.stem,
        "pdb_dir": str(protein_path.parent),
        "lig_file": cfg.ligand_file,
        "entries": entries,
        "overwrite": cfg.overwrite,
        "nproc": max(1, cfg.nproc or 1),
        "include_waters": cfg.include_waters,
        "out_ifp": cfg.out_ifp,
        "ifp_type": cfg.ifp_type,
        "ifp_levels": cfg.ifp_levels,
        "ifp_radius": cfg.ifp_radius,
        "ifp_length": cfg.ifp_length,
        "ifp_bit": cfg.ifp_bit,
        "ifp_output": cfg.ifp_output,
        "sim_matrix": cfg.sim_matrix,
        "sim_matrix_output": cfg.sim_matrix_output,
        "out_pse": cfg.out_pse,
        "pse_path": cfg.pse_path,
        "pse_interaction_types": cfg.pse_interaction_types,
        "ic_add_proximal": cfg.ic_add_proximal,
        "ic_add_atom_atom": cfg.ic_add_atom_atom,
        "ic_add_dependent_inter": cfg.ic_add_dependent_inter,
        "ic_add_h2o_pairs_with_no_target": cfg.ic_add_h2o_pairs_with_no_target,
        "ic_ignore_self_inter": cfg.ic_ignore_self_inter,
        "inter_config_overrides": cfg.inter_config_overrides,
    }
    path = wd / "_luna_api_params.json"
    path.write_text(json.dumps(params, indent=2), encoding="utf-8")
    return str(path)


def build_api_command(py_exe: str, cfg: ProjectConfig, entries: list[str]) -> list[str]:
    """Return argv to run the Python-API runner."""
    runner = write_runner(cfg.workdir)
    params = write_params(cfg.workdir, cfg, entries)
    return [py_exe, runner, params]
