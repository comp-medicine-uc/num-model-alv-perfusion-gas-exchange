# Tests

Examples and tests go here. For example, the `main` file or Jupyter notebooks.

## Files

- `main.py`: executable for simulations. To succesfully import from `src`, the file needs to be run with the `PYTHONPATH=$PYTHONPATH:.` setting (which adds the current directory to the Python PATH for module search). That is
`PYTHONPATH=$PYTHONPATH:. python3 tests/main.py`
or use another way to import from sibling directories.