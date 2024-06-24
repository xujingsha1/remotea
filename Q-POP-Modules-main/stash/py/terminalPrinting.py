# Import Modules/Libraries
import numpy as np
import itertools as itools
import textwrap
import sys

# Build Terminal Printing Code

class terminalPrint():

    def printLine(symbol='*'):
        print(symbol*(int(70/len(symbol))))

    def printErrorLine(msg):
        printLine(symbol='ERROR ')
        print('')
        terminalPrint.printNamedLine(msg, symbol=' ')
        print('')
        printLine(symbol='Error')
        sys.exit(1)

    def printSolveUpdate(msg, sec):
        colLength = 70
        msgLength = len(msg)
        str3 = ': {} seconds'.format(round(sec,4))
        str2 = ' '*(colLength - msgLength - len(str3))
        print('{}{} {}'.format(msg, str2, str3))

    def printNamedLine(name, symbol='*'):
        colLength = 70
        nameLength = len(name)
        print('{}{}{}'.format(symbol*int((colLength-nameLength)/2), name, symbol*int((colLength-nameLength)/2)))

    def printIntroMessage():

        terminalPrint.printLine(symbol='-')
        print("       ________                      ________    ________      ________   ")
        print("      |\   __  \                    |\   __  \  |\   __  \    |\   __  \  ")
        print("      \ \  \|\  \     ____________  \ \  \|\  \ \ \  \|\  \   \ \  \|\  \ ")
        print("       \ \  \\\   \   |\____________\ \ \   ____\ \ \  \\\   \   \ \   ____\ ")
        print("        \ \  \\\   \  \|____________|  \ \  \___|  \ \  \\\   \   \ \  \___|")
        print("         \ \_____  \                   \ \__\      \ \_______\   \ \__\   ")
        print("          \|___| \__\                   \|__|       \|_______|    \|__|   ")
        print("                \|__|                                                     ")
        print('')
        terminalPrint.printNamedLine('Welcome to Q-POP',' ')
        terminalPrint.printNamedLine('The Quantum Phase-field Open-source Package (Q-POP)', ' ')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        terminalPrint.printNamedLine('This Program is developed for the phase-field simulation',' ')
        terminalPrint.printNamedLine('of quantum materials and their electronic and structural transitions', ' ')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        terminalPrint.printNamedLine('Q-POP',' ')
        terminalPrint.printNamedLine('Copyright Notice (c) 2022 COMMS',' ')
        terminalPrint.printNamedLine('Created at Pennsylvania State University',' ')
        terminalPrint.printNamedLine('Version: 0.1.0',' ')
        terminalPrint.printNamedLine('Last Updated: 7/8/2022',' ')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        terminalPrint.printNamedLine('Help and User Information: Can be found in the accompanying',' ')
        terminalPrint.printNamedLine('README file and Q-POP manual.',' ')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        terminalPrint.printNamedLine('This Software is Licensed via the XXX License.', ' ')
        terminalPrint.printNamedLine('Information regarding this license can be found in ', ' ')
        terminalPrint.printNamedLine('the provided License.txt software provided with this', ' ')
        terminalPrint.printNamedLine('software or at www.gitlab.com/XXXX', ' ')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')

    def printOutroMessage():

        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        terminalPrint.printNamedLine('Copyright Notice (c) 2022 COMMS',' ')
        terminalPrint.printNamedLine('Created at Pennsylvania State University',' ')
        terminalPrint.printNamedLine('Version: 0.1.0',' ')
        terminalPrint.printNamedLine('Last Updated: 7/8/2022',' ')
        print('')
        terminalPrint.printLine(symbol='-')
        terminalPrint.printNamedLine(" ________                      ________    ________      ________   ",' ')
        terminalPrint.printNamedLine("|\   __  \                    |\   __  \  |\   __  \    |\   __  \  ",' ')
        terminalPrint.printNamedLine("\ \  \|\  \     ____________  \ \  \|\  \ \ \  \|\  \   \ \  \|\  \ ",' ')
        terminalPrint.printNamedLine("  \ \  \\\  \   |\____________\ \ \   ____\ \ \  \\\  \   \ \   ____\ ", ' ')
        terminalPrint.printNamedLine("  \ \  \\\  \  \|____________|  \ \  \___|  \ \  \\\  \    \ \  \___|",' ')
        terminalPrint.printNamedLine("   \ \_____  \                   \ \__\      \ \_______\    \ \__\   ",' ')
        terminalPrint.printNamedLine("    \|___| \__\                   \|__|       \|_______|     \|__|   ",' ')
        terminalPrint.printNamedLine("          \|__|                                                  ",' ')
        print('')
        print('')
