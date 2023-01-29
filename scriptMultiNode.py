#!/usr/bin/python
from subprocess import call
from os.path import expanduser
from os import chdir
from optparse import OptionParser
import subprocess
import re
import json
import time

# Some default values
defaultEps = 5000.0
defaultMinPts = 5
defaultDepth = 2
defaultInput = "~/homedir/cac112/input/10D_500K.ds"
defaultOutDir = "~/homedir/cac112/inputSplit"
defaultStrategy = 0

# Files in which we write measurements
PATH_MEASUREMENTS = "measurements"
PATH_INDIVIDUAL_MEASUREMENTS = "individMeasurements_v4"
PATH_PARTITION_MEASUREMENTS = "partitionMeasurements"

class BruteDBSCAN(object):
  CLUSTER_REGEX = re.compile('^\(dbScanDSU\).*?clusters.*?(\d+)$')
  RUNTIME_REGEX = re.compile('\(dbScanDSU\).*?running.*?(\d+)$')

  def __init__(self, output):
    for line in output:
      mCluster = self.CLUSTER_REGEX.match(line)
      mRuntime = self.RUNTIME_REGEX.match(line)

      if mCluster:
        self.clusters = int(mCluster.group(1))
      if mRuntime:
        self.runtime = int(mRuntime.group(1))

    if self.clusters is None:
      raise ValueError("The cluster field is not set")
    if self.runtime is None:
      raise ValueError("The runtime field is not set")

class InputSplit(object):
  READ_TIME_REGEX = re.compile('^\(inputGen\).*?reading.*?(\d+)$')
  GEN_TIME_REGEX = re.compile('^\(inputGen\).*?generating\spartition.*?(\d+)$')
  PARTITION_REGEX = re.compile('\(inputGen\).*?Partition.*?(\d+).*?noElems.*?(\d+).*?total.*?(\d+)$')
  WRITE_TIME_REGEX = re.compile('^\(inputGen\).*?writing.*?(\d+)$')

  def __init__(self, output):
    self.partitions = []
    for line in output:
      mRead = self.READ_TIME_REGEX.match(line)
      mGen = self.GEN_TIME_REGEX.match(line)
      mPart = self.PARTITION_REGEX.match(line)
      mWrite = self.WRITE_TIME_REGEX.match(line)

      if mRead:
        self.read = int(mRead.group(1))
      if mGen:
        self.gen = int(mGen.group(1))
      if mPart:
        self.partitions.append({
          "partitionIdx": int(mPart.group(1)),
          "noPoints": int(mPart.group(2)),
          "totalPoints": int(mPart.group(3))})
      if mWrite:
        self.write = int(mWrite.group(1))

    if self.read is None:
      raise ValueError("The read time is not set")
    if self.gen is None:
      raise ValueError("The generation time is not set")
    if not self.partitions:
      raise ValueError("The partitions data is missing")
    if self.write is None:
      raise ValueError("The write time parameter is not set")

    self.write = self.write - self.gen
    self.gen = self.gen - self.read

class MultiNode(object):
  LOCAL_READ_REGEX = re.compile('^\(dbScanMulti\).*?local\sread.*?(\d+)$')
  LOCAL_INDEX_REGEX = re.compile('^\(dbScanMulti\).*?local\sindex.*?(\d+)$')
  LOCAL_REGEX = re.compile('^\(dbScanMulti\).*?local\sstep.*?(\d+)$')
  GET_TREES_REGEX = re.compile('^\(dbScanMulti\).*?local\sclusters.*?(\d+)$')
  GET_FOREIGN_REGEX = re.compile('^\(dbScanMulti\).*?foreign\sedges.*?(\d+)$')
  FOREIGN_CNT_REGEX = re.compile('^\(dbScanMulti\).*?FOREIGN.*?NO.*?(\d+)$')
  RUNTIME_REGEX = re.compile('^\(dbScanMulti\).*?runtime.*?(\d+)$')

  def __init__(self, output):
    for line in output:
      mLocalRead = self.LOCAL_READ_REGEX.match(line)
      mLocalIndex = self.LOCAL_INDEX_REGEX.match(line)
      mLocalTime = self.LOCAL_REGEX.match(line)
      mTreesTime = self.GET_TREES_REGEX.match(line)
      mForeignTime = self.GET_FOREIGN_REGEX.match(line)
      mForeignCnt = self.FOREIGN_CNT_REGEX.match(line)
      mRuntime = self.RUNTIME_REGEX.match(line)

      if mLocalRead:
        self.localRead = int(mLocalRead.group(1))
      if mLocalIndex:
        self.localIndex = int(mLocalIndex.group(1))
      if mLocalTime:
        self.localTime = int(mLocalTime.group(1))
      if mTreesTime:
        self.sendTrees = int(mTreesTime.group(1))
      if mForeignTime:
        self.sendForeign = int(mForeignTime.group(1))
      if mForeignCnt:
        self.foreignEdges = int(mForeignCnt.group(1))
      if mRuntime:
        self.parallelRuntime = int(mRuntime.group(1))

    if self.localRead is None:
      raise ValueError("The local read field is not set")
    if self.localIndex is None:
      raise ValueError("The local index field is not set")
    if self.localTime is None:
      raise ValueError("The local time field is not set")
    if self.sendTrees is None:
      raise ValueError("The send tree time is not set")
    if self.sendForeign is None:
      raise ValueError("The send foreign time is not set")
    if self.foreignEdges is None:
      raise ValueError("The foreignEdges field is not set")
    if self.parallelRuntime is None:
      raise ValueError("The parallelRuntime field is not set")

class ExecutionResult(object):

  def __init__(self, eps, minPts, depth, brute, inputSplit, multiNode):
    self.eps = eps
    self.minPts = minPts
    self.noMachines = 2**depth
    self.brute = brute.__dict__
    self.inputSplit = inputSplit.__dict__
    self.multiNode = multiNode.__dict__

class PartitionExperiment(object):

  def __init__(self, speed, accuracy, eps, depth):
    self.speed = speed
    self.accuracy = accuracy
    self.eps = epsValues
    self.noMachines = 2**depth

def declareGlobals():
  global defaultEps
  global defaultMinPts
  global defaultDepth
  global defaultInput
  global defaultOutDir
  global defaultStrategy

  global PATH_MEASUREMENTS
  global PATH_INDIVIDUAL_MEASUREMENTS
  global PATH_PARTITION_MEASUREMENTS

def setupInput(inputFile, outputDir, depth, eps, strategy):
  # Compile generator
  command = ["make", "compile_gen"]
  call(command)

  # Run the generator
  command = ["make", "exec_gen"]
  command.extend(["FILE=%s" % (inputFile), \
                  "OUT_DIR=%s" % (outputDir), \
                  "DEPTH=%d" % (depth), \
                  "EPS=%lf" % (eps), \
                  "STRATEGY=%d" % (strategy)
                ])
  cmd = subprocess.Popen(" ".join(command), shell=True, stdout=subprocess.PIPE)

  return InputSplit(cmd.stdout)

def runBruteDSU(inputFile, eps, minPts):
  # Compile
  command = ["make", "compile_dsu"]
  call(command)  

  # Run
  command = ["time", "make", "exec_dsu"]
  command.extend(["FILE=%s" % (inputFile), \
                  "EPS=%lf" % (eps), \
                  "MIN_PTS=%d" % (minPts), \
                  "INDEX=1"
                ])
  cmd = subprocess.Popen(" ".join(command), shell=True, stdout=subprocess.PIPE)

  return BruteDBSCAN(cmd.stdout)

def runMultiNode(inputFile, outputDir, depth, eps, minPts, optimizedComm=0, verbose=0):
  # Compile
  command = ["make", "compile_multinode"]
  call(command)

  # Run
  command = ["make", "exec_multinode"]
  command.extend(["OUT_DIR=%s" % (outputDir), \
                  # i.e. whatever is after the last / 
                  "ORIG_FILE=%s" % (inputFile.rsplit('/', 1)[-1]), \
                  "NO_MACHINES=%d" % (2**depth), \
                  "EPS=%lf" % (eps), \
                  "MIN_PTS=%d" % (minPts), \
                  "OPTIMIZE_COMM=%d" % (optimizedComm), \
                  "VERBOSE=%d" % (verbose)
                  ])
  cmd = subprocess.Popen(" ".join(command), shell=True, stdout=subprocess.PIPE)
  
  return MultiNode(cmd.stdout)

def command_line_arguments():
  parser = OptionParser(usage="")
  parser.add_option("--EXP", dest="experiment", action="store_true",
    default=False, help="Set this option if you want to set up an experiment in which you\
    run all the different algorithms with varying parameters")

  parser.add_option("--PEXP", dest="partitionExperiment", action="store_true",
    default=False, help="Set this option if you want to set up an experiment for partition generator")

  parser.add_option("--BULK", dest="writeBulk", action="store_true",
    default=False, help="Set this if you want to write all the experiment results in a\
    single file")

  parser.add_option("--EPS", dest="eps", action="store",
    type="string", default=defaultEps)

  parser.add_option("--MP", dest="minPts", action="store", type="int",
    default=defaultMinPts, help="minPts in DBSCAN")

  parser.add_option("--D", dest="depth", action="store", type="int",
    default=defaultDepth, help="if this is d, the no. of partitions is 2 ^ d")

  parser.add_option("--F", dest="inputFile", action="store", type="string",
    default=defaultInput, help="location of the input file")

  parser.add_option("--DIR", dest="outputDir", action="store", type="string",
    default=defaultOutDir, help="location of the output directory.\
    This is where the partitions will be written to")

  parser.add_option("--GEN_S", dest="strategy", action="store", type="int",
    default=defaultStrategy, help="strategy type for input generation.\
    Choose 0 for speed, 1 for accuracy")

  return parser.parse_args()

def setupWorkspace():
  command = ["mkdir", "-p", PATH_MEASUREMENTS]
  call(command)

  command = ["mkdir", "-p", PATH_INDIVIDUAL_MEASUREMENTS]
  call(command)

def currTime():
  return int(round(time.time() * 1000))

def writeBulk(measurements, timestamp, fileName):
  writeFolder = "%s/%d" % (PATH_MEASUREMENTS, timestamp)
  command = ["mkdir", "-p", writeFolder]
  call(command)
  jsonFileResult = "%s/%s.result" % (writeFolder, fileName)

  with open(jsonFileResult, "wb") as outputFile:
    results = [measurement.__dict__ for measurement in measurements]
    jsonResults = json.dumps(results, sort_keys = True, indent = 4, separators=(',', ': '))
    outputFile.write(jsonResults)

def writeSingle(measurement, timestamp, fileName):
  writeFolder = "%s/%d" % (PATH_INDIVIDUAL_MEASUREMENTS, timestamp)
  command = ["mkdir", "-p", writeFolder]
  call(command)
  jsonFileResult = "%s/FILE_%s_EPS_%d_MP_%d_MACHINES_%d.result" % \
    (writeFolder, fileName, int(measurement.eps), measurement.minPts, measurement.noMachines)
  print "Writing in ", jsonFileResult
  with open(jsonFileResult, "wb") as outputFile:
    jsonResult = json.dumps(measurement.__dict__, sort_keys = True, indent = 4, separators=(',', ': '))
    outputFile.write(jsonResult)

def writeSinglePartition(measurement, timestamp, fileName):
  writeFolder = "%s/%d" % (PATH_PARTITION_MEASUREMENTS, timestamp)
  command = ["mkdir", "-p", writeFolder]
  call(command)
  # MP_0 is a hack for the parsing in plots to work
  jsonFileResult = "%s/FILE_%s_EPS_%d_MP_0_MACHINES_%d.result" % \
    (writeFolder, fileName, int(measurement.eps), measurement.noMachines)
  print "Writing partition in ", jsonFileResult
  with open(jsonFileResult, "wb") as outputFile:
    jsonResult = json.dumps(measurement.__dict__, sort_keys = True, indent = 4, separators=(',', ': '))
    outputFile.write(jsonResult)

def main():
  setupWorkspace()
  declareGlobals()
  (options, args) = command_line_arguments()

  # get parameters needed by other methods from looking at the options value
  inputFile = options.inputFile
  outputDir = options.outputDir
  strategy = int(options.strategy)

  if options.experiment:
    epsValues = [5000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0, 400000.0, 600000.0]
    minPtsValues = [3, 5, 8, 10]
    depthValues = [1, 2, 3, 4]

    timestamp = currTime()
    fileName = inputFile.rsplit('/', 1)[-1]

    if options.partitionExperiment:
      for eps in epsValues:
        for depth in depthValues:
          fastInput = setupInput(inputFile, outputDir, depth, eps, 0)
          accurateInput = setupInput(inputFile, outputDir, depth, eps, 1)

          partitionMeasurement = PartitionExperiment(fastInput, accurateInput, eps, depth)
          writeSinglePartition(partitionMeasurement, timestamp, fileName)
    else:
      bruteDict = {}
      inputDict = {}
      measurements = []
      # inputSplit only depends on eps and depth
      # brute only depends on eps and minPts
      # Note that it is crucial the for loops are processed in this order in order for the partitions to exist
      # when running the multinode algorithm
      for epsIdx, eps in enumerate(epsValues):
        for depthIdx, depth in enumerate(depthValues):
          for minPtsIdx, minPts in enumerate(minPtsValues):
            if (epsIdx, depthIdx) in inputDict:
              inputSplit = inputDict[(epsIdx, depthIdx)]
            else:
              inputSplit = setupInput(inputFile, outputDir, depth, eps, strategy)
              inputDict[(epsIdx, depthIdx)] = inputSplit

            if (epsIdx, minPtsIdx) in bruteDict:
              brute = bruteDict[(epsIdx, minPtsIdx)]
            else:
              brute = runBruteDSU(inputFile, eps, minPts)
              bruteDict[(epsIdx, minPtsIdx)] = brute

            multiNode = runMultiNode(inputFile, outputDir, depth, eps, minPts, verbose=1, optimizedComm=1)
            measurement = ExecutionResult(eps, minPts, depth, brute, inputSplit, multiNode)

            if options.writeBulk:
              measurements.append(measurement)
            else:
              writeSingle(measurement, timestamp, fileName)

      if options.writeBulk:
        writeBulk(measurements, timestamp, fileName)

  else:
    eps = float(options.eps)
    minPts = int(options.minPts)
    depth = int(options.depth)
    #inputResult = setupInput(inputFile, outputDir, depth, eps, strategy)
    #print inputResult.__dict__
    #bruteResult = runBruteDSU(inputFile, eps, minPts)
    #print bruteResult.__dict__
    multiNodeResult = runMultiNode(inputFile, outputDir, depth, eps, minPts, verbose=1, optimizedComm=1)
    print multiNodeResult.__dict__

if __name__ == "__main__":
  main()
