#!/bin/bash
# Plan 96 phase 2: the torch-fdmt Crab recover-a-known + benchmark on an NVIDIA GPU, plus the
# torch-dsp cross-device kernel checks. Runs as root via SSM Run Command on the jansky-gpu
# launch template (NVIDIA-driver base AMI). Usage: cuda_validation.sh <commit-sha>
# Writes everything to /opt/job/out; nothing goes back into a results/ tree on the instance.
set -euxo pipefail
SHA="$1"
TORCH_VERSION="2.13.0"   # the version uv.lock pins (CPU wheel); same version, CUDA build
export HOME=/root PATH=/root/.local/bin:$PATH
mkdir -p /opt/job/out && cd /opt/job

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > out/nvidia-smi.csv

curl -LsSf https://astral.sh/uv/install.sh | sh
git clone -q https://github.com/joebarbere/jansky-research.git repo
cd repo && git checkout -q "$SHA"
uv sync -q --extra fdmt   # the dev group (pytest) is included by default
# uv.lock pins torch to the CPU wheel index; swap in the CUDA build of the same version, and
# call .venv/bin/python directly so `uv run` cannot re-sync it back to CPU.
uv pip install -q --python .venv/bin/python --reinstall "torch==${TORCH_VERSION}" \
    --index-url https://download.pytorch.org/whl/cu128
PY=.venv/bin/python
$PY -c "import torch; assert torch.cuda.is_available(), 'no CUDA'; print(torch.__version__, torch.version.cuda, torch.cuda.get_device_name(0))" | tee out/torch.txt

# 1. The repo's own tests for the three torch slices (CPU path; proves the environment).
$PY -m pytest -q tests/test_fdmt.py tests/test_singlepulse.py tests/test_torchdsp.py 2>&1 | tail -3 | tee out/pytest.txt

# 2. CPU vs CUDA kernel parity.
$PY /opt/job/repo/infra/jobs/fdmt_parity.py > out/fdmt_parity.json

# 3. The torchfdmt slice itself: science leg on CUDA, single-session CPU+CUDA benchmark.
#    --out is a scratch dir, never the repo root (the CLAUDE.md results-clobber lesson).
$PY -m jansky_research.singlepulse --out /opt/job/sp --device cuda --benchmark \
    --bench-devices cpu,cuda > /opt/job/sp.log 2>&1
cp /opt/job/sp/results/singlepulse_metrics.json out/singlepulse_cuda.json

# 4. torch-dsp single-session benchmark + cross-device kernel checks (no download).
$PY -m jansky_research.torchdsp --out /opt/job/dsp --benchmark-only --device cuda > /opt/job/dsp.log 2>&1
cp /opt/job/dsp/results/torchdsp_metrics.json out/torchdsp_cuda.json

echo "$SHA" > out/commit.txt
echo JOB-DONE
