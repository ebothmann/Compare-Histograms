#! /usr/bin/env python

from __future__ import print_function

import yoda
import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
import yaml
import sys
import os.path

# TODO: Include parsing of arguments
# These should have a higher priority compared to the configuration file

# Obtain a configuration file if available
try:
    configurationFile = open("Plot.yaml")
except IOError:
    print("You must provide a configuration file named `Plot.yaml`.",
          file=sys.stderr)
    sys.exit()

# Parse the configuration file
try:
    configuration = yaml.safe_load(configurationFile)
except yaml.YAMLError:
    print("The configuration file could not be parsed. "
          "Please provide a valid YAML-formatted file.",
          file=sys.stderr)
    sys.exit()

# Obtain the number of bins
try:
    binRange = configuration["binRange"]
    if binRange[1] < 1 or binRange[0] < 0:
        raise Exception("The `binRange` in your configuration file "
                        "is not valid")
except KeyError:
    binRange = [0, 0]

leftBinEdges = []
binWidths = []
histograms = []
yodaDistributions = []
normalizedHistogram = None

for distribution in configuration["distributions"]:

    try:
        if distribution["hidden"]:
            continue
    except KeyError:
        pass

    fileName = distribution["file"]
    extension = os.path.splitext(fileName)[-1].lower()

    if extension == '.yoda':
        yodaDistributions.append(distribution)

    elif extension == '.dat':
        data = numpy.loadtxt(fileName)
        try:
            column = distribution["column"]
        except KeyError:
            # The column with the bin heights is not explicitly given
            # Assume that the last column provides the bin heights
            column = data.shape[1] - 1
        if binRange[1] == 0:
            # Assume we want all bins
            binRange[1] = data.shape[0]
        try:
            isNormalized = distribution["isNormalized"]
        except KeyError:
            isNormalized = False
        binHeights = data[binRange[0]:binRange[0]+binRange[1], column]
        if isNormalized:
            if (normalizedHistogram):
                raise Exception("There may only be one normalized histogram")
            normalizedHistogram = binHeights
        else:
            histograms.append(binHeights)
        if len(leftBinEdges) == 0 and column > 1 and data.shape[1] > 2:
            # Assume that the first two columns give us the bin edges
            for edges in data[binRange[0]:binRange[0]+binRange[1], :2]:
                leftBinEdges.append(edges[0])
                binWidths.append(edges[1] - edges[0])

    else:
        raise Exception("Unknown file extension detected while parsing"
                        "distribution file names.")

for yodaDistribution in yodaDistributions:
    fileName = distribution["file"]
    yodaHistos = yoda.readYODA(fileName)
    try:
        histogramIndex = distribution["histogramIndex"]
    except KeyError:
        # The histogram index is not explicitly given
        histogramIndex = 0
        if binRange[1] > 0:
            # We already have a valid bin range
            # Try to find a matching histogram
            i = 0
            for yodaHisto in yodaHistos:
                if (isinstance(yodaHisto, yoda.core.Histo1D)):
                    if (len(yodaHisto.bins) >= binRange[0] + binRange[1]):
                        histogramIndex = i
                        continue
                i += 1
            if i == len(yodaHistos):
                raise Exception("There is no YODA histogram in " + fileName +
                                " which has enough bins to accommodate the bin"
                                " range")
    yodaHistoBins = yodaHistos[histogramIndex].bins
    if binRange[1] == 0:
        binRange[1] = yodaHistoBins.len()
    try:
        isNormalized = distribution["isNormalized"]
    except KeyError:
        isNormalized = False
    yodaValues = [bin.height for bin in yodaHistoBins]
    binHeights = yodaValues[binRange[0]:binRange[0]+binRange[1]]
    if isNormalized:
        if (normalizedHistogram):
            raise Exception("There may only be one normalized histogram")
        normalizedHistogram = binHeights
    else:
        histograms.append()
    if len(leftBinEdges) == 0:
        leftBinEdges = [bin.edges[0] for bin in yodaHistoBins]
        binWidths = [bin.edges[1] - bin.edges[0] for bin in yodaHistoBins]

assert(len(leftBinEdges))
assert(len(leftBinEdges) == len(binWidths))
assert(len(histograms))

if (normalizedHistogram):
    # Normalize histograms
    for histogram in histograms:
        for i in range(len(histogram)):
            histogram[i] = histogram[i] - normalizedHistogram[i]
            histogram[i] /= normalizedHistogram[i]

fig, ax = plt.subplots()

subBinWidths = [binWidth/(len(histograms)) for binWidth in binWidths]

colorCycle = mpl.rcParams['axes.color_cycle']

style = ""
try:
    style = configuration["histogramStyle"]
except KeyError:
    pass

if style == 'step':
    rightBinEdges = [leftBinEdge + binWidth
                     for leftBinEdge, binWidth
                     in zip(leftBinEdges, binWidths)]
    x = numpy.ravel(zip(leftBinEdges, rightBinEdges))
    if style == 'step':
        ax.plot(x, [0] * len(x), color='grey')

logarithmic = False
try:
    logarithmic = configuration["logarithmic"]
except KeyError:
    pass

for i in range(len(histograms)):
    leftSubBinEdges = list(leftBinEdges)
    if logarithmic:
        histograms[i] = abs(histograms[i])
    if style == 'step':
        # y = numpy.ravel(zip(histograms[i], histograms[i]))
        # ax.plot(x, y)
        ax.step(leftSubBinEdges, histograms[i])
    else:
        for j in range(len(leftBinEdges)):
            leftSubBinEdges[j] += binWidths[j] / len(histograms) * i
        ax.bar(leftSubBinEdges,
               histograms[i],
               subBinWidths,
               color=colorCycle[i],
               linewidth=0,
               log=logarithmic)

if logarithmic and style == 'step':
    plt.yscale('log')
plt.show()
