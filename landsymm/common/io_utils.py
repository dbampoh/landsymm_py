"""I/O utilities (stub)."""
#
# TODO:
# - delete existing outputs with same semantics as MATLAB
# - preserve .mat cache behavior if needed
from __future__ import annotations

import os


def delete_existing_outputs(*args, **kwargs):
    """Delete existing outputs if needed (TODO)."""
    filename = args[0] if args else kwargs.get("filename")
    if not filename:
        raise RuntimeError("filename required")

    suffix_list = [".mat", ".maps.mat", ".garr.mat"]
    for suffix in suffix_list:
        mat_filename = f"{filename}{suffix}"
        if os.path.exists(mat_filename):
            try:
                print(f"Deleting existing {mat_filename}")
                os.remove(mat_filename)
            except Exception:
                print(
                    f"Warning: Unable to delete existing {mat_filename}. "
                    f"MATLAB may try to read that instead of the new {filename}."
                )
