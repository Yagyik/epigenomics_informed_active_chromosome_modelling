#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv}"
INSTALL_SYSTEM_DEPS=1

if [[ "${1:-}" == "--python-only" ]]; then
  INSTALL_SYSTEM_DEPS=0
fi

run_with_sudo() {
  if command -v sudo >/dev/null 2>&1; then
    sudo "$@"
  else
    "$@"
  fi
}

install_linux_deps() {
  if command -v apt-get >/dev/null 2>&1; then
    run_with_sudo apt-get update
    run_with_sudo apt-get install -y \
      build-essential \
      bc \
      gcc \
      g++ \
      make \
      openmpi-bin \
      libopenmpi-dev \
      libgsl-dev \
      pkg-config \
      python3-venv
    return
  fi

  if command -v dnf >/dev/null 2>&1; then
    run_with_sudo dnf install -y \
      bc \
      gcc \
      gcc-c++ \
      make \
      openmpi \
      openmpi-devel \
      gsl \
      gsl-devel \
      pkgconf-pkg-config \
      python3
    return
  fi

  echo "Unsupported Linux package manager. Install gsl/openmpi/bc/build tools manually." >&2
}

install_macos_deps() {
  if ! command -v brew >/dev/null 2>&1; then
    echo "Homebrew is required on macOS: https://brew.sh" >&2
    exit 1
  fi
  brew update
  brew install bc gsl open-mpi pkg-config python@3.12
}

if [[ "$INSTALL_SYSTEM_DEPS" -eq 1 ]]; then
  OS_NAME="$(uname -s)"
  case "$OS_NAME" in
    Linux) install_linux_deps ;;
    Darwin) install_macos_deps ;;
    *)
      echo "Unsupported OS ($OS_NAME) for automatic system dependency install." >&2
      echo "Install gsl/openmpi/bc/build tools manually, then re-run with --python-only." >&2
      ;;
  esac
fi

"$PYTHON_BIN" -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip
pip install -r "$ROOT_DIR/requirements.txt"

echo "Environment ready."
echo "Activate with: source \"$VENV_DIR/bin/activate\""
