#!/usr/bin/env bash
# =============================================================================
#  feng GPU compute env — ONE-TIME SETUP.  *** NOT AUTO-RUN ***
# -----------------------------------------------------------------------------
#  DO NOT source or execute this from any pipeline. Review it, then run it by
#  hand exactly once to create ~/orcd/pool/feng-gpu-env. It downloads several GB
#  of CUDA wheels + conda packages.
#
#  BEFORE RUNNING:
#    1. Check quota:  cat ~/orcd/.quota   (POOL was 14% used / ~877 GB free on
#       2026-09-03; this env is ~8-15 GB. Stay well clear of the POOL limit.)
#    2. Run it on a COMPUTE node with internet (login node is fine for the
#       download; the heavy build is pip/conda only, no GPU needed to install).
#    3. Expect 20-40 min.
#
#  WHY A SEPARATE ENV: the CPU feng-env (ViennaRNA/IntaRNA/samtools) is kept
#  clean and torch-free. GPU libs (torch/SpliceAI/Pangolin/OpenMM) with their
#  CUDA stacks live only here and are gated at runtime by
#  feng.gpu.env.check_gpu_env(), so a missing GPU env degrades gracefully.
#
#  ⚠️  COMPATIBILITY CAVEATS — resolve these BEFORE trusting a single-env build:
#    * SpliceAI (Illumina) pins TensorFlow (historically TF 1.x / <=2.4-era) and
#      Keras — it can HARD-CONFLICT with a modern torch+CUDA stack in one env.
#      RECOMMENDED: give SpliceAI its OWN sub-env (feng-gpu-env-spliceai) or run
#      it CPU-only; do not assume torch and TF-SpliceAI coexist. Pangolin
#      (also splice effect) is torch-based and lighter — prefer Pangolin as the
#      primary deep-splice scorer and treat SpliceAI as an optional cross-check.
#    * torch CUDA build MUST match the cluster driver. cudnn/9.8.0.87-cuda12 +
#      cuda/12.9 → install the cu121/cu124 torch wheel line (see below); do NOT
#      pair a cu13x torch with the cuda-12 cudnn module.
#    * OpenMM: install the conda-forge CUDA build; verify it sees the toolkit
#      with `python -m openmm.testInstallation` on a GPU node after setup.
#      ⚠️ `conda install -c conda-forge openmm cudatoolkit` AFTER the pip cu124
#      torch install can OVERWRITE torch's bundled CUDA libs (a known single-env
#      footgun) — torch.cuda may then break. Prefer a SEPARATE openmm sub-env,
#      or install openmm FIRST and torch last, and re-check
#      `torch.cuda.is_available()` after every conda step.
#    * Pin exact versions once a working combo is found and record them here.
# =============================================================================
set -euo pipefail

echo "This script is a documented one-time installer and is intentionally guarded."
echo "Remove the following exit to run it, after reading the caveats above."
exit 0   # <-- SAFETY GUARD: delete this line to actually run the setup.

# ----------------------------------------------------------------------------
# Everything below is the reference recipe (guarded by the exit above).
# ----------------------------------------------------------------------------
GPU_ENV_PREFIX="${FENG_GPU_ENV_PREFIX:-${HOME}/orcd/pool/feng-gpu-env}"
PY_VERSION="${FENG_GPU_PY:-3.11}"

module load miniforge/25.11.0-0
module load cuda/12.9.1
module load cudnn/9.8.0.87-cuda12
CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# --- idempotent: skip if the env already imports torch -----------------------
if [[ -x "${GPU_ENV_PREFIX}/bin/python" ]] \
   && "${GPU_ENV_PREFIX}/bin/python" -c "import torch" 2>/dev/null; then
    echo "GPU env already present + torch imports: ${GPU_ENV_PREFIX}"
    exit 0
fi

echo "Creating GPU env at ${GPU_ENV_PREFIX} (python ${PY_VERSION})..."
conda create -y -p "${GPU_ENV_PREFIX}" "python=${PY_VERSION}"
conda activate "${GPU_ENV_PREFIX}"

# --- VALIDATED COMBO (installed + imported 2026-09-03, env 6.2 GB) -----------
#   python=3.11.16  openmm=8.2.0.dev  torch=2.6.0+cu124  pangolin=1.0.2
#   (from github.com/tkzeng/Pangolin). Ordering below is the FIXED order that
#   actually worked: OpenMM/cudatoolkit FIRST, torch cu124 LAST, so the conda
#   step cannot clobber torch's bundled CUDA libs (re-checked import after each).

# --- OpenMM FIRST (MD; conda-forge CUDA build) -------------------------------
conda install -y -c conda-forge openmm cudatoolkit
python -c "import openmm; print('openmm', openmm.version.version)"

# --- torch LAST (cu124 wheel; matches cuda/12.9 + cudnn9.8) -------------------
pip install --index-url https://download.pytorch.org/whl/cu124 \
    torch torchvision torchaudio
python -c "import torch, openmm; print('torch', torch.__version__)"  # openmm must still import

# --- Pangolin (primary deep-splice; torch-based) -----------------------------
# The pip name `pangolin-splice` may not match; the source of record is
# github.com/tkzeng/Pangolin (imports as `pangolin`; needs gffutils/pyfasta).
pip install pangolin-splice || \
    pip install "git+https://github.com/tkzeng/Pangolin.git" || \
    echo "NOTE: install Pangolin from github.com/tkzeng/Pangolin if both fail"

# --- SpliceAI: SEPARATE env recommended (TF conflict — see caveats) ----------
# Do NOT co-install into this torch env without resolving the TF/Keras pin.
#   conda create -y -p "${GPU_ENV_PREFIX}-spliceai" python=3.8
#   conda activate "${GPU_ENV_PREFIX}-spliceai"
#   pip install spliceai tensorflow==2.4.1

echo "GPU env build complete. Validate on a GPU node:"
echo "  ${GPU_ENV_PREFIX}/bin/python -c 'import torch; print(torch.cuda.is_available())'"
echo "  ${GPU_ENV_PREFIX}/bin/python -m openmm.testInstallation"
