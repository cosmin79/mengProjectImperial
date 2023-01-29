#!/usr/bin/python
from subprocess import call
from os.path import expanduser
from os import chdir
from optparse import OptionParser
import re

SMPF_DIR = "inputSmpf"
INPUT_DIR = "input"

def declareGlobals():
  global SMPF_DIR


def command_line_arguments(parser):
  parser.add_option("--file", dest="file", action="store",
    type="string")

  return parser.parse_args()

def main():
  declareGlobals()
  parser = OptionParser(usage="")
  (options, args) = command_line_arguments(parser)

  if not options.file:
    parser.error("You are required to provide a password field")
  file = options.file

  inputPath = "%s/%s" % (SMPF_DIR, file)
  outputPath = "%s/%s" % (INPUT_DIR, file)
  with open(inputPath, "r") as fileR:
    lines = fileR.readlines()
    with open(outputPath, "w") as fileW:
      for idx, line in enumerate(lines):
        array = [int(x) for x in line.split()]
        array = filter(lambda x: x > 0, array)
        if idx == 0:
          fileW.write("%d %d\n" % (len(lines), len(array)))

        array = [str(idx + 1)] + map(lambda x: str(x), array)
        #print array
        fileW.write(" ".join(array) + "\n")

if __name__ == "__main__":
  main()
