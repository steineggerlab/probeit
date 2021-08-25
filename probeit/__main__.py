#!/usr/bin/env python
from .classmodule import Probeit
import sys

def main():
    args = sys.argv[1:]
    probeit = Probeit(args)
    probeit.checkArgs()
    quit()


if __name__ == '__main__':
    main()
