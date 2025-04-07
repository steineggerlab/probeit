#!/usr/bin/env python
import os
from .Posnegset import PosNegSet
from .SNP import SNP
from .Primer import Primer

class Probeit:
    args = []
    subWork = None

    def __init__(self, args):
        self.args = args

    def do(self):
        if self.args == []:
            self.printUsage()

        workflow = self.args[0]
        match workflow:
            case '-h' | '--help':
                self.printUsage()

            case 'posnegset':
                print('CURRENT: ', os.getcwd())
                subWork = PosNegSet(self.args[1:])
                subWork.run()
                return

            case 'snp':
                print('CURRENT: ', os.getcwd())
                subWork = SNP(self.args[1:])
                subWork.run()
                return

            # case 'primer':
            #     print('CURRENT: ', os.getcwd())
            #     subWork = Primer(self.args[1:])
            #     subWork.run()
            #     return

            case _:
                self.printUsage()

    @staticmethod
    def printUsage():
        print("PROBEIT V2.2")
        print('probeit <workflow> [<args>]')
        print("WORKFLOWS")
        print("posnegset: make two-sets probes with positive and negative sets")
        print("snp: make two-sets probes with wild type genome, strain genome and SNPs")
        quit()