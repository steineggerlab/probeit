#!/usr/bin/env python
from .classmodule import Probeit
import sys

def main():
    args = sys.argv[1:]
    probeitProcess = Probeit(args)
    probeitProcess.checkArgs()
    

if __name__ == '__main__':
    main()
    
    
    
    
