#!/usr/bin/env python
from .classmodule import Probeit
import sys


def main():
    try:
        args = sys.argv[1:]
        probeit = Probeit(args)
        probeit.checkArgs()
        return 0
    except:
        return 1


if __name__ == '__main__':
    sys.exit(main())
