#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

REPO_DIR="${HOME}/Repositories"
PYTHON_VERSION="3.9"

# -----------------------------
# Install uv if missing
# -----------------------------

if ! command -v uv >/dev/null 2>&1; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="${HOME}/.local/bin:${PATH}"
else
    echo "uv already installed."
fi

# -----------------------------
# Clone repositories
# -----------------------------

mkdir -p "${REPO_DIR}"
cd "${REPO_DIR}"

if [ ! -d "pyHiM" ]; then
    git clone https://github.com/pyHi-M/pyHiM.git
else
    echo "pyHiM repository already exists."
fi

if [ ! -d "traceratops" ]; then
    git clone https://github.com/pyHi-M/traceratops.git
else
    echo "traceratops repository already exists."
fi

# -----------------------------
# Create fresh uv environment
# -----------------------------

cd "${REPO_DIR}/pyHiM"

if [ -d ".venv" ]; then
    echo "Removing existing .venv..."
    rm -rf .venv
fi

uv venv .venv --python "${PYTHON_VERSION}"
source .venv/bin/activate

# -----------------------------
# Install pinned scientific stack
# -----------------------------

uv pip install --upgrade pip wheel
uv pip install "setuptools<70"

uv pip install \
    "numpy==1.26.4" \
    "scipy==1.13.1" \
    "astropy==6.0.1" \
    "photutils==1.11.0" \
    "scikit-image==0.19.2" \
    "scikit-learn==1.6.1" \
    "pandas==2.2.3" \
    "dask==2024.8.0" \
    "distributed==2024.8.0" \
    "numba==0.60.0" \
    "tifffile==2024.8.30" \
    "pillow==11.1.0" \
    "matplotlib" \
    "pytest==8.3.5"

uv pip install \
    "csbdeep==0.8.1" \
    "stardist==0.9.1" \
    "tensorflow==2.19.0" \
    "tensorflow-addons==0.23.0" \
    "tensorflow-io-gcs-filesystem==0.37.1"

uv pip install \
    "roipoly==0.5.3" \
    "rich==13.9.4" \
    "requests==2.32.4" \
    "pyyaml==6.0.2" \
    "tqdm==4.67.1" \
    "cloudpickle==3.1.1" \
    "pympler==1.1" \
    "dataclasses-json==0.6.7"

# -----------------------------
# Install traceratops and pyHiM
# -----------------------------

uv pip install -e "${REPO_DIR}/traceratops"
uv pip install -e "${REPO_DIR}/pyHiM"

# -----------------------------
# Clone apiFISH
# -----------------------------

cd "${REPO_DIR}"

if [ ! -d "apiFISH" ]; then
    git clone https://github.com/apiFISH/apiFISH.git
else
    echo "apiFISH repository already exists."
fi

cd "${REPO_DIR}/apiFISH"
git checkout development

# -----------------------------
# Optional PATH/PYTHONPATH setup
# -----------------------------

echo ""
echo "Add this to your ~/.bashrc if pyHiM.py is not found:"
echo ""
echo "export PATH=\"\$PATH:${REPO_DIR}/pyHiM/src\""
echo "export PATH=\"\$PATH:${REPO_DIR}/pyHiM/src/toolbox/file_handling\""
echo "export PATH=\"\$PATH:${REPO_DIR}/pyHiM/src/postProcessing\""
echo "export PYTHONPATH=\"\$PYTHONPATH:${REPO_DIR}/pyHiM/src\""
echo "export PYTHONPATH=\"\$PYTHONPATH:${REPO_DIR}/apiFISH\""
echo "export MPLBACKEND=agg"
echo ""

# -----------------------------
# Test installation
# -----------------------------

python -c "import numpy, scipy, astropy; print('numpy', numpy.__version__); print('scipy', scipy.__version__); print('astropy', astropy.__version__)"
python -c "import traceratops; print('traceratops OK')"
python -c "import pkg_resources; print('pkg_resources OK')"

echo ""
echo "Installation complete."
echo "Activate with:"
echo "source ${REPO_DIR}/pyHiM/.venv/bin/activate"
echo ""
echo "Test pyHiM with:"
echo "python ${REPO_DIR}/pyHiM/src/pyHiM.py --help"
