# Snap package

This directory contains the Linux Snapcraft recipe for HIP2LInterActomics.

Build on Linux:

```bash
sudo snap install snapcraft --classic
./dist/linux/build_snap.sh
```

Install the generated snap locally:

```bash
sudo snap install ./hip2linteractomics_0.1_amd64.snap --classic --dangerous
hip2linteractomics
```

The snap is Linux-only. Windows does not support `.snap` as a native installer
format; use `dist/windows/install_hip2linteractomics.ps1` on Windows, or run the
Linux installer inside WSL2 when parallel execution is required.

The snap packages the GUI and launchers. The LUNA runtime remains in the
separate `luna-env` conda environment, created either by the GUI setup tab or by
`dist/linux/install_hip2linteractomics.sh`.
