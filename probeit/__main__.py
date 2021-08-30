#!/usr/bin/env python
from .classmodule import Probeit
import sys


def main(argv):
    try:
        probeitProcess = Probeit(argv[1:])
        probeitProcess.checkArgs()
        return 0
    except:
        return 1


if __name__ == '__main__':
    main(sys.argv)
