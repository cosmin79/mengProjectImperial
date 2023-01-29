#!/usr/bin/python
from subprocess import call
from optparse import OptionParser
import subprocess
import re
import json
import os
from sets import Set
import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, axes, pie, title, show

PATH_MEASUREMENTS = "measurements"
PATH_INDIVIDUAL_MEASUREMENTS = "individMeasurements_v3"
PATH_PARTITIONING_MEASUREMENTS = "partitionMeasurements"
PATH_PLOTS = "figures/v3"
PATH_EPS_PLOTS = PATH_PLOTS + "/eps"
PATH_MIN_PTS_PLOTS = PATH_PLOTS + "/minPts"
PATH_INPUT_PLOTS = PATH_PLOTS + "/inputSplit"
PATH_INPUT_ANALYSIS_PLOTS = PATH_PLOTS + "/inputSplitAnalysis"
PATH_CLUSTER_ANALYSIS = PATH_PLOTS + "/clusterAnalysis"
PATH_SCALABILITY_ANALYSIS = PATH_PLOTS + "/scalability"
PATH_DISTRIBUTED_ANALYSIS = PATH_PLOTS + "/distributedAlgo"
PATH_DIMENSIONALITY_ANALYSIS = PATH_PLOTS + "/dimensionality"
PATH_PARTITIONING_ANALYSIS = PATH_PLOTS + "/partitioning"
colors = ['b', 'y', 'r', 'g', 'c', 'm', 'k']

class Measurement(object):

  def __init__(self, file, eps, minPts, noMachines, measurement):
    self.file = file
    self.eps = eps
    self.minPts = minPts
    self.noMachines = noMachines
    # Note that the measurement object corresponds to a json object (see scriptMultiNode.py for more details)
    self.measurement = measurement

  def __eq__(self, other):
    return self.file == other.file and self.eps == other.eps and \
      self.minPts == other.minPts and self.noMachines == other.noMachines

  def __hash__(self):
    return hash((self.file, self.eps, self.minPts, self.noMachines))

def declareGlobals():
  global PATH_MEASUREMENTS
  global PATH_INDIVIDUAL_MEASUREMENTS
  global PATH_PARTITIONING_MEASUREMENTS
  global PATH_PLOTS
  global PATH_EPS_PLOTS
  global PATH_MIN_PTS_PLOTS
  global PATH_INPUT_ANALYSIS_PLOTS
  global PATH_CLUSTER_ANALYSIS
  global PATH_SCALABILITY_ANALYSIS
  global PATH_DISTRIBUTED_ANALYSIS
  global PATH_DIMENSIONALITY_ANALYSIS
  global PATH_PARTITIONING_ANALYSIS

def setupWorkspace():
  for folder in [PATH_MEASUREMENTS, PATH_INDIVIDUAL_MEASUREMENTS, PATH_PARTITIONING_MEASUREMENTS, PATH_EPS_PLOTS,\
      PATH_MIN_PTS_PLOTS, PATH_INPUT_PLOTS, PATH_INPUT_ANALYSIS_PLOTS, PATH_CLUSTER_ANALYSIS,\
      PATH_SCALABILITY_ANALYSIS, PATH_DISTRIBUTED_ANALYSIS, PATH_DIMENSIONALITY_ANALYSIS, PATH_PARTITIONING_ANALYSIS]:
    command = ["mkdir", "-p", folder]
    call(command)

# Result is a list of tuples (path, filename) where path is the absolute path of the file in question
def getFilesInDir(dirPath):
  absoluteDirPath = os.path.join(os.getcwd(), dirPath)

  result = []
  for root, dirs, files in os.walk(absoluteDirPath):
    for file in files:
      absoluteFilePath = os.path.join(root, file)
      if not os.path.isfile(absoluteFilePath) or absoluteFilePath.rsplit('.')[-1] != "result":
        continue
      fileName = re.match("^(.*)\.result$", file).group(1)
      result.append((absoluteFilePath, fileName))

  return result

def indexIndividualResults(jsonData, fileName):
  fileEncoding = "FILE_%s_EPS_%d_MP_%d_MACHINES_%d.result"
  # jsonData is an array of measurements
  for measurement in jsonData:
    newFile = os.path.join(PATH_INDIVIDUAL_MEASUREMENTS, fileEncoding % 
      (fileName, int(measurement["eps"]), measurement["minPts"], measurement["noMachines"]))
    print newFile
    with open(newFile, "wb") as outputFile:
      jsonMeasurement = json.dumps(measurement, sort_keys = True, indent = 4, separators=(',', ': '))
      outputFile.write(jsonMeasurement)

def createIndividualResults(files):
  for path, fileName in files:
    with open(path) as file:
      data = json.load(file)
      indexIndividualResults(data, fileName)

def updateIndividualResults():
  createIndividualResults(getFilesInDir(PATH_MEASUREMENTS))

def parseFromNormalMeasurements():
  measurements = []
  for path, fileName in getFilesInDir(PATH_MEASUREMENTS):
    with open(path) as file:
      data = json.load(file)
      for measurement in data:
        measurements.append(Measurement(fileName, 
          int(measurement["eps"]), measurement["minPts"], measurement["noMachines"], measurement))

  return measurements

def parseFromIndividualMeasurements(dirPath):
  measurements = []
  for path, fileName in getFilesInDir(dirPath):
    m = re.match("^FILE_(.*?)_EPS_(\d+)_MP_(\d+)_MACHINES_(\d+)$", fileName)
    (name, eps, minPts, noMachines) = (m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)))
    with open(path) as file:
      data = json.load(file)
      measurements.append(Measurement(name, eps, minPts, noMachines, data))

  return measurements

def parseFromIndividualOrNormal():
  return list(Set(parseFromNormalMeasurements() + parseFromIndividualMeasurements()))

def getDistinctFiles(measurements):
  files = Set()
  for measurement in measurements:
    files.add(measurement.file)

  return list(sorted(files))

def getEpsValues(measurements, file):
  epsValues = Set()
  for measurement in measurements:
    if measurement.file == file:
      epsValues.add(measurement.eps)

  return list(sorted(epsValues))

def getMinPtsValues(measurements, file):
  minPtsValues = Set()
  for measurement in measurements:
    if measurement.file == file:
      minPtsValues.add(measurement.minPts)

  return list(sorted(minPtsValues))

def getMeasurement(measurements, file, eps, minPts, noMachines):
  for measurement in measurements:
    if measurement.file == file and measurement.eps == eps and\
      measurement.minPts == minPts and measurement.noMachines == noMachines:
      return measurement
  raise ValueError("Corresponding entry not found")

def setAxisAndLabels(xValues, xLabel, yLabel, font=20, legend=True):
  plt.gca().set_xlim([0, max(xValues) * 1.1])
  if legend:
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', numpoints=1, fancybox=True, shadow=True)
  plt.xlabel(xLabel, fontweight='bold', fontsize=font)
  plt.ylabel(yLabel, fontweight='bold', fontsize=font)

def averageParams(file, measurements, func, eps=None, minPts=None, depth=None):
  epsValues = getEpsValues(measurements, file) if eps is None else [eps]
  minPtsValues = getMinPtsValues(measurements, file) if minPts is None else [minPts]
  depthValues = range(1, 5) if depth is None else [depth]

  sum = 0
  for eps in epsValues:
    for minPts in minPtsValues:
      for depth in depthValues:
        sum += func(getMeasurement(measurements, file, eps, minPts, 2**depth))
  sum /= (len(epsValues) * len(minPtsValues) * len(depthValues))

  return sum

def plotEpsGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  for minPts in minPtsValues:
    outputFile = "%s/FILE_%s_MIN_PTS_%d.png" % (PATH_EPS_PLOTS, file, minPts)

    plt.suptitle("Values for input %s with minPts param: %d" % (file, minPts), fontsize=25, fontweight='bold')
    # The first subfigure corresponds to the runtime using the brute or 2, 4, 8, 16 machines
    plt.subplot(311)
    runtimeBrute = [0] * len(epsValues)
    sqrtEpsValues = map(lambda x: x**0.5, epsValues)
    for depth in range(1, 5):
      runtime = []
      for index, eps in enumerate(epsValues):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        runtime.append(float(entry.measurement["multiNode"]["parallelRuntime"]) / 1000)
        runtimeBrute[index] += entry.measurement["brute"]["runtime"]

      plt.plot(sqrtEpsValues, runtime, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))
    
    # Add the graph for Brute DBSCAN
    # Note that we divide by 1000 to get the time in millis + divide by 4 (number of depth params) to average the time
    # across those runs
    runtimeBrute = map(lambda x: float(x) / (1000 * 4), runtimeBrute)
    plt.plot(sqrtEpsValues, runtimeBrute, marker='o', linestyle='-', color=colors[4], label="Brute DBScan")
    setAxisAndLabels(sqrtEpsValues, "EPS", "Runtime(s)")

    # The second subfigure corresponds to the communication volume
    plt.subplot(312)
    for depth in range(1, 5):
      commVolume = []
      for eps in epsValues:
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        commVolume.append(entry.measurement["multiNode"]["foreignEdges"])

      plt.plot(sqrtEpsValues, commVolume, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))

    plt.yscale('log')
    setAxisAndLabels(sqrtEpsValues, "EPS", "Communication volume")

    # The third subfigure corresponds to the parallel overhead (inputGen)
    plt.subplot(313)
    for depth in range(1, 5):
      runtime = []
      for eps in epsValues:
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        inputGen = entry.measurement["inputSplit"]
        runtime.append(float(inputGen["read"] + inputGen["gen"] + inputGen["write"]) / 1000)

      plt.plot(sqrtEpsValues, runtime, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))

    setAxisAndLabels(sqrtEpsValues, "EPS", "Generation overhead(s)")

    # Final touches
    plt.tight_layout()
    plt.gcf().set_size_inches(18, 18)
    plt.subplots_adjust(top=0.95)
    plt.savefig(outputFile, bbox_inches='tight', dpi=100)
    plt.clf()

def plotMinPtsGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  for eps in epsValues:
    outputFile = "%s/FILE_%s_EPS_%d.png" % (PATH_MIN_PTS_PLOTS, file, eps)

    plt.suptitle("Values for input %s with EPS param(squared): %d" % (file, eps), fontsize=25, fontweight='bold')
    # The first subfigure corresponds to the runtime using the brute or 2, 4, 8, 16 machines
    plt.subplot(311)
    runtimeBrute = [0] * len(minPtsValues)
    for depth in range(1, 5):
      runtime = []
      for index, minPts in enumerate(minPtsValues):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        runtime.append(float(entry.measurement["multiNode"]["parallelRuntime"]) / 1000)
        runtimeBrute[index] += entry.measurement["brute"]["runtime"]

      plt.plot(minPtsValues, runtime, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))
    
    # Add the graph for Brute DBSCAN
    # Note that we divide by 1000 to get the time in millis + divide by 4 (number of depth params) to average the time
    # across those runs
    runtimeBrute = map(lambda x: float(x) / (1000 * 4), runtimeBrute)
    plt.plot(minPtsValues, runtimeBrute, marker='o', linestyle='-', color=colors[4], label="Brute DBScan")
    setAxisAndLabels(minPtsValues, "minPts", "Runtime(s)")

    # The second subfigure corresponds to the communication volume
    plt.subplot(312)
    for depth in range(1, 5):
      commVolume = []
      for minPts in minPtsValues:
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        commVolume.append(entry.measurement["multiNode"]["foreignEdges"])

      plt.plot(minPtsValues, commVolume, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))

    setAxisAndLabels(minPtsValues, "minPts", "Communication volume")

    # The third subfigure corresponds to the parallel overhead (inputGen)
    plt.subplot(313)
    for depth in range(1, 5):
      runtime = []
      for minPts in minPtsValues:
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        inputGen = entry.measurement["inputSplit"]
        runtime.append(float(inputGen["read"] + inputGen["gen"] + inputGen["write"]) / 1000)

      plt.plot(minPtsValues, runtime, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))

    setAxisAndLabels(minPtsValues, "minPts", "Generation overhead(s)")

    # Final touches
    plt.tight_layout()
    plt.gcf().set_size_inches(18, 18)
    plt.subplots_adjust(top=0.95)
    plt.savefig(outputFile, bbox_inches='tight', dpi=100)
    plt.clf()

def plotInputGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  labels = ("read", "generate", "write")
  for eps in epsValues:
    outputFile = "%s/FILE_%s_EPS_%d.png" % (PATH_INPUT_PLOTS, file, eps) 

    plt.suptitle("Partitioning components for file: %s, EPS(sq): %d" % (file, eps), size=14, fontweight='bold')
    read = [0.0] * 4; gen = [0.0] * 4; write = [0.0] * 4
    for depth in range(1, 5):
      for minPts in minPtsValues:
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth).measurement["inputSplit"]
        read[depth - 1] += entry["read"]; gen[depth - 1] += entry["gen"]; write[depth - 1] += entry["write"]
      read[depth - 1] /= len(minPtsValues) * 1000; gen[depth - 1] /= len(minPtsValues) * 1000; write[depth - 1] /= len(minPtsValues) * 1000

    objects = ("2 machines", "4 machines", "8 machines", "16 machines")
    y_pos = np.arange(len(objects))
    width = 0.35

    p1 = plt.bar(y_pos, read, width, color=colors[0])
    p2 = plt.bar(y_pos, gen, width, color=colors[1], bottom=read)
    p3 = plt.bar(y_pos, write, width, color=colors[2], bottom=[r + g for r, g in zip(read, gen)])
    plt.ylabel('Partitioning time', fontsize=15, fontweight='bold')
    plt.xticks(y_pos + width / 2., objects)
    plt.legend((p1[0], p2[0], p3[0]), ('read', 'gen', 'write'),\
      bbox_to_anchor=(1, 0.5), loc='center left', numpoints=1, fancybox=True, shadow=True)
    #plt.legend((p1[0], p2[0], p3[0]), ('read', 'gen', 'write'))

    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.savefig(outputFile, bbox_inches='tight', dpi=100)
    plt.clf()

def plotInputAnalysisGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  sqrtEpsValues = map(lambda x: x**0.5, epsValues)

  clusterAnalysisFile = "%s/FILE_%s_clustering.png" % (PATH_INPUT_ANALYSIS_PLOTS, file)
  phantomAnalaysisFile = "%s/FILE_%s_phantom.png" % (PATH_INPUT_ANALYSIS_PLOTS, file)

  # The first subplot corresponds to the stdev of the partition sizes
  variance = [0] * 4
  for depth in range(1, 5):
    for eps in epsValues:
      varForEps = 0.0
      for minPts in minPtsValues:
        entries = getMeasurement(measurements, file, eps, minPts, 2**depth).measurement["inputSplit"]["partitions"]
        meanPointsPartition = float(sum(map(lambda entry: entry["noPoints"], entries))) / len(entries)
        varForEps += sum(map(lambda entry: (float(entry["noPoints"]) - meanPointsPartition)**2, entries))**0.5
      variance[depth - 1] += varForEps / len(minPtsValues)
    variance[depth - 1] /= len(epsValues)

  objects = ("2 machines", "4 machines", "8 machines", "16 machines")
  y_pos = np.arange(len(objects))
  plt.bar(y_pos, variance, align = "center", alpha=0.5)
  plt.xticks(y_pos, objects, fontweight='bold')
  plt.ylabel('Stdev for partitioning', fontsize=15, fontweight='bold')
  plt.title("Analysis of the partitioning quality for file: %s" % file, fontsize=20, fontweight='bold', y=1.08)
  # Final touches
  plt.tight_layout()
  plt.savefig(clusterAnalysisFile, bbox_inches='tight', dpi=100)
  plt.clf()

  # The second plot corresponds to the total no. of phantom points
  plt.title("Phantom points for file: %s" % file, fontsize=20, fontweight='bold', y=1.08)
  for depth in range(1, 5):
    phantomPoints = [0] * len(epsValues)
    for epsIdx, eps in enumerate(epsValues):
      for minPts in minPtsValues:
        entries = getMeasurement(measurements, file, eps, minPts, 2**depth).measurement["inputSplit"]["partitions"]
        phantomPoints[epsIdx] += sum(map(lambda entry: entry["totalPoints"] - entry["noPoints"], entries))
      phantomPoints[epsIdx] /= len(minPtsValues)

    plt.plot(sqrtEpsValues, phantomPoints, marker='o', linestyle='-', color=colors[depth - 1], label="%d machines" % (2**depth))
    setAxisAndLabels(sqrtEpsValues, "EPS", "Number of phantom points", font=15)
  plt.yscale('log')

  # Final touches
  plt.tight_layout()
  plt.savefig(phantomAnalaysisFile, bbox_inches='tight', dpi=100)
  plt.clf()

def plotClusterAnalysisGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  sqrtEpsValues = map(lambda x: x**0.5, epsValues)
  # The first type of graphs plot EPS against the no of clusters
  for minPts in minPtsValues:
    outputFile = "%s/FILE_%s_MIN_PTS_%d.png" % (PATH_CLUSTER_ANALYSIS, file, minPts)

    plt.title("Clusters for file %s, minPts: %d" % (file, minPts), fontsize=20, fontweight='bold', y=1.08)
    clusters = [0] * len(epsValues)
    for depth in range(1, 5):
      for index, eps in enumerate(epsValues):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        clusters[index] += entry.measurement["brute"]["clusters"]
    clusters = map(lambda x: x / 4, clusters)

    plt.plot(sqrtEpsValues, clusters, marker='o', linestyle='-', color=colors[0])
    setAxisAndLabels(sqrtEpsValues, "EPS", "Number of clusters", legend=False)

    # Final touches
    plt.tight_layout()
    plt.savefig(outputFile, bbox_inches='tight', dpi=100)
    plt.clf()

  # The second type of graphs plot MIN_PTS against the no of clusters for ct. EPS
  for eps in epsValues:
    outputFile = "%s/FILE_%s_EPS_%d.png" % (PATH_CLUSTER_ANALYSIS, file, eps)

    plt.title("Clusters for file %s, EPS(sq): %d" % (file, eps), fontsize=20, fontweight='bold', y=1.08)
    clusters = [0] * len(minPtsValues)
    for depth in range(1, 5):
      for index, minPts in enumerate(minPtsValues):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        clusters[index] += entry.measurement["brute"]["clusters"]
    clusters = map(lambda x: x / 4, clusters)

    plt.plot(minPtsValues, clusters, marker='o', linestyle='-', color=colors[0])
    setAxisAndLabels(minPtsValues, "minPts", "Number of clusters", legend=False)

    # Final touches
    plt.tight_layout()
    plt.savefig(outputFile, bbox_inches='tight', dpi=100)
    plt.clf()

def formatFile(file):
  regex = re.compile('^(.*?)(\.ds)?$')
  result = regex.match(file).group(1).replace("_", "-")

  return result

def formatTitle(file):
  regex = re.compile('^(.*?)(_synthetic)?(\.ds)?$')
  result = (regex.match(file).group(1) + ".ds")

  return result

def createDirs(listDirs):
  for folder in listDirs:
    command = ["mkdir", "-p", folder]
    call(command)

def copyFileFolder(filePath, destFolder, fileName):
  command = ["cp", filePath, "%s/%s" % (destFolder, fileName)]
  call(command)

def plotScalabilityGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)

  for eps in epsValues:
    for minPts in minPtsValues:
      fileName = "FILE-%s-EPS-%d-MINPTS-%d.png" % (formatFile(file), eps, minPts)

      brute = 0.0
      runtime = [0.0] * 4
      for depth in range(1, 5):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth)
        brute += entry.measurement["brute"]["runtime"]
        runtime[depth - 1] = entry.measurement["multiNode"]["parallelRuntime"]
      brute /= 4
      speedup = map(lambda x: float(brute) / x, runtime)

      objects = ("2 machines", "4 machines", "8 machines", "16 machines")
      y_pos = np.arange(len(objects))
      plt.bar(y_pos, speedup, align = "center", alpha=0.5)
      plt.xticks(y_pos, objects, fontweight='bold', fontsize=35)
      plt.ylabel('Speedup', fontweight='bold', fontsize=35)
      plt.title("Speedup for file: %s, EPS(sq): %d, MIN PTS: %d" % (formatTitle(file), eps, minPts), fontsize=35, fontweight='bold', y=1.08)

      # Final touches
      plt.gcf().set_size_inches(18, 12)
      plt.rcParams['ytick.labelsize'] = 40

      minPtsDirPath = "%s/minPts-%d" % (PATH_SCALABILITY_ANALYSIS, minPts)
      epsDirPath = "%s/eps-%d" % (PATH_SCALABILITY_ANALYSIS, eps)
      createDirs([minPtsDirPath, epsDirPath])

      # save in minPts dir
      minPtsPath = "%s/%s" % (minPtsDirPath, fileName)

      plt.savefig(minPtsPath, bbox_inches='tight', dpi=100)
      plt.clf()
      # copy it to the corresponding eps dir
      copyFileFolder(minPtsPath, epsDirPath, fileName)
      
def sumList(listA, listB):
  return [a + b for a, b in zip(listA, listB)]

def plotDistributedGraphs(file, measurements):
  epsValues = getEpsValues(measurements, file)
  minPtsValues = getMinPtsValues(measurements, file)
  for eps in epsValues:
    for minPts in minPtsValues:
      fileName = "FILE-%s-EPS-%d-MINPTS-%d.png" % (formatFile(file), eps, minPts)

      plt.suptitle("Breakdown for file: %s, EPS(sq): %d, MIN_PTS: %d" % (file, eps, minPts), size=35, fontweight='bold')
      localRead = [0.0] * 4; localIndex = [0.0] * 4; algorithm = [0.0] * 4; send = [0.0] * 4; merge = [0.0] * 4
      for depth in range(1, 5):
        entry = getMeasurement(measurements, file, eps, minPts, 2**depth).measurement["multiNode"]
        localRead[depth - 1] = entry["localRead"]
        localIndex[depth - 1] = entry["localIndex"] - entry["localRead"]
        algorithm[depth - 1] = entry["localTime"] - entry["localIndex"]
        send[depth - 1] = (entry["sendForeign"] - entry["localTime"])
        merge[depth - 1] = (entry["parallelRuntime"] - entry["sendForeign"])

      objects = ("2 machines", "4 machines", "8 machines", "16 machines")
      y_pos = np.arange(len(objects))
      #print merge

      sumSoFar = [0.0] * 4
      p1 = plt.bar(y_pos, localRead, align = "center", color=colors[0], bottom=sumSoFar)
      sumSoFar = sumList(sumSoFar, localRead)

      p2 = plt.bar(y_pos, localIndex, align = "center", color=colors[1], bottom=sumSoFar)
      sumSoFar = sumList(sumSoFar, localIndex)

      p3 = plt.bar(y_pos, algorithm, align = "center", color=colors[2], bottom=sumSoFar)
      sumSoFar = sumList(sumSoFar, algorithm)

      p4 = plt.bar(y_pos, send, align = "center", color=colors[3], bottom=sumSoFar)
      sumSoFar = sumList(sumSoFar, send)

      p5 = plt.bar(y_pos, merge, align = "center", color=colors[4], bottom=sumSoFar)

      plt.ylabel('Time (ms)', fontsize=40, fontweight='bold')
      plt.xticks(y_pos, objects, fontweight='bold', fontsize=35)
      plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('read local', 'construct index', 'algorithm', 'communication', 'merge'), fontsize=35)

      minPtsDirPath = "%s/minPts-%d" % (PATH_DISTRIBUTED_ANALYSIS, minPts)
      epsDirPath = "%s/eps-%d" % (PATH_DISTRIBUTED_ANALYSIS, eps)
      createDirs([minPtsDirPath, epsDirPath])

      # Final touches
      plt.tight_layout()
      plt.subplots_adjust(top=0.88)
      plt.gcf().set_size_inches(18, 12)
      plt.rcParams['ytick.labelsize'] = 40

      # save in minPts dir
      minPtsPath = "%s/%s" % (minPtsDirPath, fileName)
      plt.savefig(minPtsPath, bbox_inches='tight', dpi=100)
      plt.clf()
      # copy it to the corresponding eps dir
      copyFileFolder(minPtsPath, epsDirPath, fileName)

def plotDimensionalityCorrelation(file10, file100, measurements):
  epsValues = getEpsValues(measurements, file10)
  minPts = 5

  plt.rcParams['ytick.labelsize'] = 40
  plt.rcParams['xtick.labelsize'] = 40
  for depth in range(1, 5):
    fileName = "%s/DEPTH-%d-minPts-5.png" % (PATH_DIMENSIONALITY_ANALYSIS, depth)

    map10 = {}
    map100 = {}
    xLabelMap = {}
    for idx, eps in enumerate(epsValues):
      map10[eps] = getMeasurement(measurements, file10, eps, minPts, 2**depth).measurement["multiNode"]["parallelRuntime"] / 10
      map100[eps] = getMeasurement(measurements, file100, eps, minPts, 2**depth).measurement["multiNode"]["parallelRuntime"] / 100
      xLabelMap[idx] = "%dK" % (eps / 1000)
    
    p1 = plt.scatter(xLabelMap.keys(), map10.values(), color=colors[0], s=[20 * 4 * 3 for n in epsValues])
    p2 = plt.scatter(xLabelMap.keys(), map100.values(), color=colors[1], s=[20 * 4 * 3 for n in epsValues])
    plt.xticks(xLabelMap.keys(), xLabelMap.values())

    print "DEPTH: ", depth
    best = 0.0
    for eps in epsValues:
      best = max(best, 1.0 * map10[eps] / map100[eps])
    print best

    plt.title("Correlation between files 10D_500K and 100D-500K", fontsize=35, fontweight='bold', y=1.08)
    plt.ylabel('Normalized runtime(ms)', fontweight='bold', fontsize=35)
    plt.xlabel('EPS', fontweight='bold', fontsize=35)
    plt.legend((p1, p2), ("10D_500K", "100D-500K"), fontsize=35, scatterpoints=1, loc='lower left', ncol=5)

    # Final touches
    plt.gcf().set_size_inches(24, 18)
    plt.savefig(fileName, bbox_inches='tight', dpi=100)
    plt.clf()

def plotGraphs(measurements):
  # special type of graph: correlation between 10D - 500K and 100D - 500K
  #plotDimensionalityCorrelation("10D_500K_synthetic.ds", "100D_500K.ds", measurements)

  for file in getDistinctFiles(measurements):
    print file
    #plotScalabilityGraphs(file, measurements)
    plotDistributedGraphs(file, measurements)
    #plotClusterAnalysisGraphs(file, measurements)
    #plotInputAnalysisGraphs(file, measurements)
    #plotInputGraphs(file, measurements)
    #plotEpsGraphs(file, measurements)
    #plotMinPtsGraphs(file, measurements)

def plotPartitions(files, measurements):
  destDir = PATH_PARTITIONING_ANALYSIS + "/simplePartitions"
  createDirs([destDir])
  minPts = 0
  # Same eps values are used for all the files
  epsValues = getEpsValues(measurements, files[0])

  plt.rcParams['xtick.labelsize'] = 40
  plt.rcParams['ytick.labelsize'] = 40
  for depth in range(1, 5):
    fileName = "%s/DEPTH-%d.png" % (destDir, depth)
    noPartitions = 2**depth

    xLabelMap = {}
    for partition in range(noPartitions):
      xLabelMap[partition] = partition

    plotRefs = []
    for fileIdx, file in enumerate(files):
      val = [0] * noPartitions
      for eps in epsValues:
        partitions = getMeasurement(measurements, file, eps, minPts, noPartitions).measurement["accuracy"]["partitions"]
        for partition in partitions:
          val[partition["partitionIdx"]] += partition["noPoints"]
      val = map(lambda x: x / len(epsValues), val)

      p = plt.scatter(xLabelMap.values(), val, color=colors[fileIdx], s=[20 * 4 * 3 for n in range(noPartitions)])
      plotRefs.append(p)
    
    plt.xticks(xLabelMap.values(), xLabelMap.keys())
    plt.legend(tuple(plotRefs), tuple(files), fontsize=30, scatterpoints=1, loc='lower left', ncol=5)
    plt.title("Partition sizes with no. partitions: %d" % noPartitions, fontsize=35, fontweight='bold', y=1.08)
    plt.ylabel('No. elems / partition', fontweight='bold', fontsize=35)
    plt.xlabel('Partition idx', fontweight='bold', fontsize=35)
    plt.yscale('log')

    # Final touches
    plt.gcf().set_size_inches(18, 16)
    plt.savefig(fileName, bbox_inches='tight', dpi=100)
    plt.clf()
      
def plotOverallPartitions(files, measurements):
  destDir = PATH_PARTITIONING_ANALYSIS + "/partitions"
  createDirs([destDir])
  minPts = 0
  # Same eps values are used for all the files
  epsValues = getEpsValues(measurements, files[0])

  plt.rcParams['xtick.labelsize'] = 40
  plt.rcParams['ytick.labelsize'] = 40
  for file in files:
    for depth in range(1, 5):
      fileName = "%s/FILE-%s-DEPTH-%d.png" % (destDir, formatFile(file), depth)
      noPartitions = 2**depth

      xLabelMap = {}
      avgPartPoints = [0] * len(epsValues); avgTotalPoints = [0] * len(epsValues)
      for epsIdx, eps in enumerate(epsValues):
        partitions = getMeasurement(measurements, file, eps, minPts, noPartitions).measurement["speed"]["partitions"]
        sumPart = 0; sumTotal = 0
        for partition in partitions:
          sumPart += partition["noPoints"]; sumTotal += partition["totalPoints"]
        avgPartPoints[epsIdx] = sumPart / noPartitions; avgTotalPoints[epsIdx] = sumTotal / noPartitions
        xLabelMap[epsIdx] = "%dK" % (eps / 1000)

      p1 = plt.scatter(xLabelMap.keys(), avgPartPoints, color=colors[0], s=[20 * 4 * 3 for eps in epsValues])
      p2 = plt.scatter(xLabelMap.keys(), avgTotalPoints, color=colors[1], s=[20 * 4 * 3 for eps in epsValues])

      plt.xticks(xLabelMap.keys(), xLabelMap.values())
      plt.legend((p1, p2), ("partition points", "total points"), fontsize=30, scatterpoints=1, loc='lower left', ncol=5)
      plt.title("Partition analysis for file %s" % file, fontsize=35, fontweight='bold', y=1.08)
      plt.ylabel('Avg. no. elems / partition', fontweight='bold', fontsize=35)
      plt.xlabel('EPS', fontweight='bold', fontsize=35)
      plt.yscale('log')

      # Final touches
      plt.gcf().set_size_inches(18, 16)
      plt.savefig(fileName, bbox_inches='tight', dpi=100)
      plt.clf()

def plotPartitioningGraphs(measurements):
  filesOfInterest = ["100D_300K.ds", "100D_500K.ds", "15D_2M.ds"]
  plotPartitions(filesOfInterest, measurements)
  plotOverallPartitions(filesOfInterest, measurements)

def main():
  declareGlobals()
  setupWorkspace()
  #updateIndividualResults()

  plotGraphs(parseFromIndividualMeasurements(PATH_INDIVIDUAL_MEASUREMENTS))
  #plotPartitioningGraphs(parseFromIndividualMeasurements(PATH_PARTITIONING_MEASUREMENTS))
  
if __name__ == "__main__":
  main()
