"""Entry point for the ``allopep`` command and ``python -m allopep``."""

import sys

from .cli import main


def run():
    """Run the pipeline, reporting any failure as exit code 1."""
    try:
        return main()
    except Exception as exc:
        print(f'AlloPep failed: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(run())
