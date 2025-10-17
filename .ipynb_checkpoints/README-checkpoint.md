# NSRL
NSRL lattice file description and uniform beam optimization using SciBmad ecosystem.
Tutorial notebooks for IJulia and IPython usage.


## Installation of Xopt for PythonCall within a Conda environment

```
conda install -c conda-forge xopt

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"

cat > "$CONDA_PREFIX/etc/conda/activate.d/julia_python.sh" <<'EOF'
export JULIA_PYTHONCALL_EXE="$CONDA_PREFIX/bin/python"
EOF

echo 'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"' > "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh"
```
