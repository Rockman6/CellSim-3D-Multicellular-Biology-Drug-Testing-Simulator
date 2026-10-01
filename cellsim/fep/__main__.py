"""Module entry point for `python -m cellsim.fep <args>`."""

from __future__ import annotations

import sys

from cellsim.fep import main


if __name__ == "__main__":
    sys.exit(main())
