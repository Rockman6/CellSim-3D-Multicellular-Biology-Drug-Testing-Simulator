"""Module entry point for `python -m src.fep <args>`."""

from __future__ import annotations

import sys

from src.fep import main


if __name__ == "__main__":
    sys.exit(main())
