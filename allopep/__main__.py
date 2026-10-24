"""Entry point for ``python -m allopep``."""

import sys

from .cli import main

if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'AlloPep failed: {exc}', file=sys.stderr)
        sys.exit(1)
