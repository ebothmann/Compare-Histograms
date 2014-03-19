#! /usr/bin/env python

from __future__ import print_function

import yoda
import numpy
import matplotlib.pyplot as plt
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
    for distribution in configuration["distributions"]:
        fileName = distribution["file"]
        extension = os.path.splitext(fileName)[-1].lower()
        if extension == '.yoda':
            yodaHistos = yoda.readYODA(fileName)
            if yodaHistos.len() == 0:
                bins = yodaHistos[0].bins
                binRange[1] = bins.len()
        if extension == '.dat':
            dataRows = numpy.loadtxt(fileName)
            binRange[1] = dataRows.shape[0]
        else:
            raise Exception("Unknown file extension detected while parsing"
                            "distribution file names.")
        if binRange[1] != 0:
            break
    if binRange[1] < 1 or binRange[0] < 0:
        raise Exception("The bin range could not be determined by looking "
                        "at the distribution files. Consider to explicitly "
                        "provide a range.")

sys.exit()

# Obtain bins from YODA histogram
yodaHistos = yoda.readYODA("Rivet.yoda")
yodaHisto = yodaHistos[0]
yodaBins = yodaHisto.bins

# Obtain YODA values
yodaValues = [bin.height for bin in yodaBins]

# Obtain bar edges and widths
# (assumed to be the same as in the convolution data)
leftEdges = [bin.edges[0] for bin in yodaBins]
widths = [bin.edges[1] - bin.edges[0] for bin in yodaBins]

# Obtain convolution data values
convolutionValues = numpy.loadtxt('Convolute.dat', usecols=(5,))
print(convolutionValues)

# Calculate normalized differences
differences = []
for i in range(len(yodaValues)):
    differences.append((yodaValues[i] - convolutionValues[i])/yodaValues[i])

plt.bar(leftEdges, differences, widths)

plt.show()
