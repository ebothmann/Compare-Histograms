#! /usr/bin/env python
# coding=utf-8

from __future__ import print_function

# System modules
import os.path
import argparse
import sys
import math

# Third party modules
import yoda
import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
import yaml
import grapefruit
import seaborn as sns

def openConfiguration(configurationFileName):
    try:
        configurationFile = open(configurationFileName)
    except IOError:
        print(configurationFileName + " could not be opened.",
              file=sys.stderr)
        sys.exit()
    return configurationFile


def parseConfiguration(configurationFileName):
    configurationFile = openConfiguration(configurationFileName)
    try:
        configuration = yaml.safe_load(configurationFile)
    except yaml.YAMLError:
        print("The configuration file could not be parsed. "
              "Please provide a valid YAML-formatted file.",
              file=sys.stderr)
        sys.exit()
    return configuration


def readConfiguration(configuration, key=None, default=None):
    try:
        value = configuration[key]
    except KeyError:
        value = default
    return value


def hasFileExtension(distribution, extension=''):
    fileName = distribution["file"]
    fileExtension = os.path.splitext(fileName)[-1].lower()
    return fileExtension == '.' + extension


def isDat(distribution):
    return hasFileExtension(distribution, 'dat')


def isYODA(distribution):
    return hasFileExtension(distribution, 'yoda')


def readHistogram(distribution, binRange, leftBinEdges, binWidths, scalingReversed, overallScaling):
    fileName = distribution["file"]
    print("Read " + fileName)
    if isDat(distribution):
        data = numpy.loadtxt(fileName)
        column = readConfiguration(distribution,
                                   key='column',
                                   default=2)
        if binRange[1] == 0:
            # Assume we want all bins
            binRange[1] = data.shape[0]
        binHeights = data[binRange[0]:binRange[0]+binRange[1], column]
        if len(leftBinEdges) == 0 and column > 1 and data.shape[1] > 2:
            # Assume that the first two columns give us the bin edges
            for edges in data[binRange[0]:binRange[0]+binRange[1], :2]:
                leftBinEdges.append(edges[0])
                binWidths.append(edges[1] - edges[0])
    elif isYODA(distribution):

        yodaHistos = yoda.readYODA(fileName)
        try:
            histogramName = distribution["histogramName"]
        except KeyError:
            # The histogram index is not explicitly given
            histogramName = ""
            if binRange[1] > 0:
                # We already have a valid bin range
                # Try to find a matching histogram
                for yodaHistoName in yodaHistos:
                    yodaHisto = yodaHistos[yodaHistoName]
                    if (isinstance(yodaHisto, yoda.core.Histo1D)):
                        lastBin = binRange[0] + binRange[1]
                        if (len(yodaHisto.bins) >= lastBin):
                            histogramName = yodaHistoName
                            continue
                if len(histogramName) == 0:
                    raise Exception("There is no YODA histogram in "
                                    + fileName +
                                    " which has enough bins "
                                    + "to accommodate the bin range")

        yodaHistoBins = yodaHistos[histogramName].bins
        if binRange[1] == 0:
            binRange[1] = len(yodaHistoBins)
        yodaValues = [bin.height for bin in yodaHistoBins]
        binHeights = yodaValues[binRange[0]:binRange[0]+binRange[1]]
        if scalingReversed:
            # Manually read out scale_factor (there seems to be no interface in yoda)
            print(histogramName)
            scale_factor = None
            with open(fileName) as f:
                is_correct_histogram = False
                for line in f:
                    if is_correct_histogram and line[:9] == 'ScaledBy=':
                        scale_factor = float(line[9:])
                        break
                    if line[-29:-1] == 'D0_2008_S7554427/' + 'd01-x01-y01':
                        is_correct_histogram = True
            assert(scale_factor is not None)
            print(yodaHistos[histogramName].numEntries())
            scale_factor = scale_factor * yodaHistos[histogramName].numEntries() / 100000000.0
            print('Division by ', scale_factor)
            binHeights = [h / scale_factor for h in binHeights]

        print('Multiplication by ', overallScaling)
        binHeights = [h * overallScaling for h in binHeights]
        if len(leftBinEdges) == 0:
            leftBinEdges = [bin.xEdges[0] for bin in yodaHistoBins]
            binWidths = [bin.xEdges[1] - bin.xEdges[0]
                         for bin in yodaHistoBins]
    else:
        raise Exception("Unknown file extension detected while parsing"
                        " distribution file names.")

    scaledBy = readConfiguration(distribution,
                                 key='scaledBy',
                                 default=None)

    if args.normalize:
        normalizationFactor = sum(binHeights)
        binHeights = [h / normalizationFactor for h in binHeights]
    elif scaledBy is not None:
        binHeights = [h / scaledBy for h in binHeights]

    return binRange, leftBinEdges, binWidths, binHeights


def plot(configurationFileName):
    configuration = parseConfiguration(configurationFileName)

    # Obtain the number of bins
    try:
        binRange = configuration["binRange"]
        if binRange[1] < 1 or binRange[0] < 0:
            raise Exception("The `binRange` in your configuration file "
                            "is not valid")
    except KeyError:
        binRange = [0, 0]

    yRange = readConfiguration(configuration,
                               key='yRange',
                               default=None)
    xRange = readConfiguration(configuration,
                               key='xRange',
                               default=None)


    leftBinEdges = []
    binWidths = []
    histograms = []
    labels = []
    normalizedHistogram = None

    distributions = sorted(configuration["distributions"],
                           key=lambda distribution: isYODA(distribution))

    scalingReversed = readConfiguration(configuration,
                                        key='scalingReversed',
                                        default=False)

    overallScaling = readConfiguration(configuration,
                                       key='overallScaling',
                                       default=1.0)

    for distribution in distributions:
        try:
            if distribution["hidden"]:
                continue
        except KeyError:
            pass

        fileName = distribution["file"]

        binRange, leftBinEdges, binWidths, binHeights = readHistogram(
            distribution,
            binRange,
            leftBinEdges,
            binWidths,
            scalingReversed,
            overallScaling)

        try:
            normalizeBy = distribution["normalizeBy"]
        except KeyError:
            normalizeBy = None

        if normalizeBy:
            distributionToNormalizeBy = None
            for otherDistribution in distributions:
                if otherDistribution["file"] == normalizeBy:
                    distributionToNormalizeBy = otherDistribution
                    break
            if (distributionToNormalizeBy is None):
                raise Exception("The distribution to normalize by "
                                + "with the file name "
                                + normalizeBy
                                + " does not exist.")
            binHeightsToNormalizeBy = readHistogram(
                distributionToNormalizeBy,
                binRange,
                leftBinEdges,
                binWidths,
                scalingReversed,
                overallScaling)[3]
            binHeights = [binHeight / otherBinHeight # / 93710965.0 * 100000000 # (93710965.0 / 106979031.0) # - 1.0
                          for binHeight, otherBinHeight
                          in zip(binHeights, binHeightsToNormalizeBy)]

        isNormalized = readConfiguration(distribution,
                                         key='isNormalized',
                                         default=False)

        if isNormalized:
            print("Normalized dat distribution detected")
            if (normalizedHistogram is not None):
                raise Exception("There may only be "
                                "one normalized histogram")
            normalizedHistogram = binHeights
        else:
            histograms.append(binHeights)
            label = readConfiguration(distribution,
                                      key='label',
                                      default=fileName)
            labels.append(label)

        

    assert(len(leftBinEdges))
    assert(len(leftBinEdges) == len(binWidths))
    assert(len(histograms))

    if (normalizedHistogram is not None):
        # Normalize histograms
        for histogram in histograms:
            for i in range(len(histogram)):
                if normalizedHistogram[i] == 0.0:
                    print("Normalized bin height is 0")
                    histogram[i] = 1.0
                else:
                    histogram[i] /= normalizedHistogram[i]

    fig, ax = plt.subplots()

    subBinWidths = [binWidth/(len(histograms)) for binWidth in binWidths]

    colorCycle = mpl.rcParams['axes.color_cycle']

    style = readConfiguration(configuration,
                              key='histogramStyle',
                              default='')

    if style == 'step':
        rightBinEdges = [leftBinEdge + binWidth
                         for leftBinEdge, binWidth
                         in zip(leftBinEdges, binWidths)]
    #     x = numpy.ravel(zip(leftBinEdges, rightBinEdges))
    #     ax.plot(x, [0] * len(x), color='grey', zorder=-2)

    logarithmic = readConfiguration(configuration,
                                    key='logarithmic',
                                    default=False)


    plotMode = readConfiguration(configuration,
                                 key='plotMode',
                                 default=None)

    # Calculate average error bands for replica plots
    if plotMode == 'replica':
        averages = []
        stdDevs = []
        numberOfBins = len(histograms[0])
        for bin in range(numberOfBins):
            average = 0.0
            stdDev = 0.0
            for histogram in histograms:
                average += histogram[bin]
            average /= numberOfBins
            for histogram in histograms:
                stdDev += (histogram[bin] - average)**2
            stdDev /= numberOfBins
            stdDev = math.sqrt(stdDev)
            averages.append(average)
            stdDevs.append(stdDev)
        histograms.append(averages)
        histograms.append([average + stdDev
                           for average, stdDev
                           in zip(averages, stdDevs)])
        histograms.append([average - stdDev
                           for average, stdDev
                           in zip(averages, stdDevs)])
        labels.append('µ')
        labels.append('µ + σ')
        labels.append('µ - σ')

    for i in range(len(histograms)):
        leftSubBinEdges = list(leftBinEdges)
        if logarithmic:
            histograms[i] = [abs(binHeight) for binHeight in histograms[i]]
        if style == 'step':
            # y = numpy.ravel(zip(histograms[i], histograms[i]))
            # ax.plot(x, y)
            leftSubBinEdges.append(rightBinEdges[-1])
            try:
                histogramList = histograms[i].tolist()
            except AttributeError:
                histogramList = histograms[i]
            histogramList.append(histogramList[-1])
            if plotMode == 'replica':
                # color = baseColorNames[0]
                if i >= len(histograms) - 3:
                    if i == len(histograms) - 3:
                        linewidth = 2.0
                    else:
                        linewidth = 1.0
                    ax.step(leftSubBinEdges,
                            histogramList,
                            'black',
                            alpha=1.0,
                            linewidth=linewidth,
                            where='post',
                            zorder=1,
                            # color=baseColorNames[5],
                            label=labels[i])
                else:
                    ax.step(leftSubBinEdges,
                            histogramList,
                            'g',
                            alpha=0.5,
                            where='post',
                            zorder=-1,
                            color=color,
                            label='Replicas')
            else:
                # color = baseColorNames[i]
                ax.step(leftSubBinEdges,
                        histogramList,
                        where='post',
                        # color=color,
                        label=labels[i])
        else:
            for j in range(len(leftBinEdges)):
                leftSubBinEdges[j] += binWidths[j] / len(histograms) * i
            ax.bar(leftSubBinEdges,
                   histograms[i],
                   subBinWidths,
                   color=colorCycle[i],
                   linewidth=0,
                   log=logarithmic,
                   label=labels[i])

    if plotMode == 'replica':
        x = []
        sigmaMinusStdDevs = []
        sigmaPlusStdDevs = []
        for leftBinEdge, binWidth, upper, lower in zip(leftBinEdges,
                                                       binWidths,
                                                       histograms[-2],
                                                       histograms[-1]):
            x.append(leftBinEdge)
            x.append(leftBinEdge + binWidth)
            sigmaMinusStdDevs.append(lower)
            sigmaMinusStdDevs.append(lower)
            sigmaPlusStdDevs.append(upper)
            sigmaPlusStdDevs.append(upper)
        ax.fill_between(x,
                        sigmaPlusStdDevs,
                        sigmaMinusStdDevs,
                        # facecolor=baseColorNames[5],
                        edgecolor='none',
                        alpha=0.3,
                        zorder=0)

    if logarithmic and style == 'step':
        plt.yscale('log')

    xLogarithmic = readConfiguration(configuration,
                                     key='xLogarithmic',
                                     default=False)
    if xLogarithmic:
        plt.xscale('log')

    if yRange:
        plt.ylim(yRange)

    if xRange:
        plt.xlim(xRange)

    legendHidden = readConfiguration(configuration,
                                     key='legendHidden',
                                     default=False)

    if not legendHidden:
        legendLocation = readConfiguration(configuration,
                                           key='legendLocation',
                                           default='best')
        handles, labels = ax.get_legend_handles_labels()
        if plotMode == 'replica':
            handles = handles[-4:-2]
            labels = ["Replica", "$\mu \pm \sigma$"]
        ax.legend(handles, labels, loc=legendLocation)
        ax.get_legend().set_title('Closure Tests')
        # ax.get_legend().set_bbox_to_anchor((0.1,0.075))

    title = readConfiguration(configuration,
                              key='title',
                              default=None)
    if title is not None:
        plt.suptitle(title)

    xLabel = readConfiguration(configuration,
                               key='xLabel',
                               default=None)
    if xLabel is not None:
        ax.set_xlabel(xLabel)
    yLabel = readConfiguration(configuration,
                               key='yLabel',
                               default=None)
    if yLabel is not None:
        ax.set_ylabel(yLabel)

    ax.get_yaxis().get_major_formatter().set_useOffset(False)

    plt.savefig(os.path.splitext(configurationFileName)[0] + ".pdf")

# TODO: Include parsing of arguments
# These should have a higher priority compared to the configuration file

parser = argparse.ArgumentParser()
parser.add_argument("files", default="Plot.yaml", nargs='*')
parser.add_argument("-n", "--normalize", action="store_true")
parser.add_argument("-d", "--delete-bins", type=int, default=0)
parser.add_argument("-l", "--usetex", action="store_true")
args = parser.parse_args()

if args.usetex:
    def figure_size_from_width(width):
        """Returns a single plot figure size in inches given a width in points"""
        inches_per_point = 1.0/72.27
        golden_mean = (math.sqrt(5)-1.0)/2.0
        inches_width = width * inches_per_point
        fig_height = inches_width*golden_mean
        return [inches_width,fig_height]

    # Got from LaTeX output using \showthe\linewidth
    line_width = 510.0  

    # Fits two figs on one page with 3 lines of caption each
    fig_line_width_fraction = 0.85  

    fig_size = figure_size_from_width(line_width * fig_line_width_fraction)

    # Set parameters
    params = {'font.size': 9,
              # 'axes.labelsize': 10,
              # 'text.fontsize': 10,
              # 'legend.fontsize': 10,
              # 'xtick.labelsize': 8,
              # 'ytick.labelsize': 8,
              'font.family': 'serif',
              'text.usetex': True,
              'text.latex.preamble': [r"\usepackage{amstext}",
                                      r"\usepackage{siunitx}"],
              'figure.figsize': fig_size}
    plt.rcParams.update(params)

# baseColorNames = [grapefruit.Color.RgbToHtml(100./255.,177./255.,209./255.),  # '#a6cee3',
#                   '#1f78b4',
#                   grapefruit.Color.RgbToHtml(142./255.,223./255.,87./255.),   # '#b2df8a',
#                   '#33a02c',
#                   '#fb9a99',
#                   '#e31a1c',
#                   '#fdbf6f',
#                   '#ff7f00',
#                   grapefruit.Color.RgbToHtml(220./255.,105./255.,220./255.),
#                   grapefruit.Color.RgbToHtml(160./255.,105./255.,214./255.),  # '#cab2d6',
#                   '#6a3d9a',
#                   # '#ffff99',
#                   grapefruit.Color.RgbToHtml(165./255.,119./255.,6./255.),
#                   '#b15928',
#                   grapefruit.Color.RgbToHtml(128./255.,64./255.,0./255.)
#                   ]

sns.set_context('paper')  # Contexts are paper, notebook, talk and poster
for fileName in args.files:
    plot(fileName)
