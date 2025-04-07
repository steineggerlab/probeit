#!/usr/bin/env python
from .Probeit import Probeit
import sys

def main():
    args = sys.argv[1:]
    probeitProcess = Probeit(args)
    probeitProcess.do()
    return 0
    

if __name__ == '__main__':
    main()
    
    
    
    
