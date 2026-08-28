"""Small, strict utilities for editing BMM Fortran namelists.

The previous iSKYLAB workflow replaced exact strings copied from a particular
namelist template.  Small changes in comments/spacing therefore caused silent
failures, including aerosol PSDs not being updated.  These helpers edit values
by variable name and fail loudly if an expected variable/block is absent.

This is intentionally not a general Fortran namelist parser.  It implements the
limited operations needed by the BMM batch scripts while preserving the rest
of the template verbatim.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable


def _format_scalar(value) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, Path):
        return repr(str(value))
    if isinstance(value, str):
        # Caller can pass an already quoted Fortran string if desired.
        if (value.startswith("'") and value.endswith("'")) or (
            value.startswith('"') and value.endswith('"')
        ):
            return value
        return repr(value)
    return f"{value}"


def set_value(text: str, name: str, value, *, required: bool = True) -> str:
    """Replace one scalar/array assignment by variable name.

    ``name`` may include an explicit Fortran slice, e.g. ``n_aer1(1:3,1:1)``.
    The replacement is restricted to a single physical line, which is how the
    BMM template stores the control/aerosol assignments edited by this repo.
    """
    # Stop at the next namelist assignment on the same line, at a newline, or
    # at the group terminator.  Some historical BMM templates place two array
    # assignments on one physical line (e.g. n_aer1 followed by d_aer1).
    pattern = re.compile(
        rf"(?m)(?P<indent>^[ \t]*|(?<=[,])){re.escape(name)}\s*=\s*"
        rf"[^\n]*?(?=(?:[A-Za-z_]\w*(?:\([^\n=]*?\))?\s*=)|$)"
    )
    matches = list(pattern.finditer(text))
    if not matches:
        if required:
            raise KeyError(f"Namelist variable not found: {name}")
        return text
    if len(matches) != 1:
        raise ValueError(f"Expected one assignment for {name}, found {len(matches)}")
    replacement = rf"\g<indent>{name} = {_format_scalar(value)},"
    return pattern.sub(replacement, text, count=1)


def set_array(text: str, name: str, values: Iterable[float]) -> str:
    """Replace one one-line numeric array assignment."""
    values = list(values)
    value_text = ", ".join(f"{float(v):.12g}" for v in values)
    pattern = re.compile(
        rf"(?m)(?P<indent>^[ \t]*|(?<=[,])){re.escape(name)}\s*=\s*"
        rf"[^\n]*?(?=(?:[A-Za-z_]\w*(?:\([^\n=]*?\))?\s*=)|$)"
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise KeyError(f"Expected exactly one namelist assignment for {name}; found {len(matches)}")
    return pattern.sub(rf"\g<indent>{name} = {value_text},", text, count=1)



def set_literal_array(text: str, name: str, values) -> str:
    """Replace a one-line array assignment with arbitrary Fortran literals.

    Useful for logical and character arrays, for which :func:`set_array`'s
    numeric formatting is inappropriate.
    """
    def fmt(v):
        if isinstance(v, bool):
            return ".true." if v else ".false."
        if isinstance(v, str):
            if (v.startswith("'") and v.endswith("'")) or (v.startswith('"') and v.endswith('"')):
                return v
            return repr(v)
        return f"{v}"
    value_text = ", ".join(fmt(v) for v in values)
    pattern = re.compile(
        rf"(?m)(?P<indent>^[ \t]*|(?<=[,])){re.escape(name)}\s*=\s*"
        rf"[^\n]*?(?=(?:[A-Za-z_]\w*(?:\([^\n=]*?\))?\s*=)|$)"
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise KeyError(f"Expected exactly one namelist assignment for {name}; found {len(matches)}")
    return pattern.sub(rf"\g<indent>{name} = {value_text},", text, count=1)

def replace_group(text: str, group: str, body: str) -> str:
    """Replace an entire ``&group ... /`` block.

    This is used for chamber time-series data, whose array lengths change from
    experiment to experiment and should never be matched against a hard-coded
    1853-element template string.
    """
    pattern = re.compile(
        rf"(?ims)^\s*&{re.escape(group)}\b.*?/"
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise KeyError(f"Expected exactly one &{group} group; found {len(matches)}")
    block = f"&{group}\n{body.rstrip()}\n/"
    return pattern.sub(block, text, count=1)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(path)
    return path.read_text()
