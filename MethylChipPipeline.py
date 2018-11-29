#!/usr/bin/env python3

import re,time,os,sys
import traceback
import argparse
import subprocess
import configparser
import logging
from logging.handlers import RotatingFileHandler
import glob
from datetime import datetime
import gzip
import yaml


class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def MethylChip(args, scriptpath):
    config = args.config
    yamldict = dict()
    with open(config, 'r', encoding = "utf-8") as f:
        conf = f.read()
        yamldict = yaml.load(conf)
    print(yamldict)
    Import = str(yamldict["import"]).lower()
    matrixfile = yamldict["matrixfile"]
    Analysis = yamldict["analysis"]
    detPfile = yamldict["detPfile"]
    chipType = yamldict["chiptype"]
    Pdfile = yamldict["pdfile"]
    Directory = yamldict["directory"]
    Loadmethod = yamldict["loadmethod"]
    detPcut = yamldict["detPcut"]
    probeCutoff = yamldict["probecutoff"]
    sampleCutoff = yamldict["samplecutoff"]
    Population = str(yamldict["population"]).upper()
    filterXY = str(yamldict["filterXY"]).upper()
    normMethod = yamldict["normMethod"]
    Cores = yamldict["cores"]
    diffexprProbeP = yamldict["diffexprProbeP"]
    dmrMethod = yamldict["dmrMethod"]
    minProbe = yamldict["minProbe"]
    diffexprRegionP = yamldict["diffexprRegionP"]
    dmrCores = yamldict["dmrcores"]
    maxGap = yamldict["maxGap"]
    adjPprobe = yamldict["adjPprobe"]
    bpSpan = yamldict["bpSpan"]
    compGroup = " ".join(yamldict["compGroup"])
    gseaMethod = yamldict["gseaMethod"]
    gseaPvalue = yamldict["gseaPvalue"]
    controlGroup = yamldict["controlGroup"]
    freqThreshold = yamldict["freqThreshold"]
    sampleType = yamldict["sampleType"]
    print("Rscript {scriptpath}/Methyl_Chip_Analysis.R --import {Import} --pdfile {Pdfile} --matrixfile "
                     "{matrixfile} --arraytype {chipType} --directory {Directory} --loadmethod {Loadmethod} "
                     "--probecutoff {probeCutoff} --samplecutoff {sampleCutoff} --detPcut {detPcut} --population {Population} "
                     "--filterXY {filterXY} --detPfile {detPfile} --analysis {Analysis} --normmethod {normMethod} "
                     "--cores {Cores} --diffexprprobeP {diffexprProbeP} --dmrmethod {dmrMethod} --minProbe {minProbe} "
                     "--diffexprregionP {diffexprRegionP} --dmrcores {dmrCores} --maxGap {maxGap} --adjPprobe {adjPprobe} "
                     "--bpSpan {bpSpan} --compareGroup '{compGroup}' --gseaMethod {gseaMethod} --gseaPvalue {gseaPvalue} "
                     "--controlGroup {controlGroup} --freqThreshold {freqThreshold} --sampleType {sampleType}".format(scriptpath = scriptpath, Import = Import,
                      Pdfile = Pdfile, matrixfile = matrixfile, chipType = chipType, Directory = Directory, Loadmethod = Loadmethod,
                      probeCutoff = probeCutoff, sampleCutoff = sampleCutoff, detPcut = detPcut, Population = Population,
                      filterXY = filterXY, detPfile = detPfile, Analysis = Analysis, normMethod = normMethod,
                      Cores = Cores, diffexprProbeP = diffexprProbeP, dmrMethod = dmrMethod, minProbe = minProbe,
                      diffexprRegionP = diffexprRegionP, dmrCores = dmrCores, maxGap = maxGap, adjPprobe = adjPprobe,
                      bpSpan = bpSpan, compGroup = compGroup, gseaMethod = gseaMethod, gseaPvalue = gseaPvalue,
                      controlGroup = controlGroup, freqThreshold = freqThreshold, sampleType = sampleType))
    subprocess.Popen("Rscript {scriptpath}/Methyl_Chip_Analysis.R --import {Import} --pdfile {Pdfile} --matrixfile "
                     "{matrixfile} --arraytype {chipType} --directory {Directory} --loadmethod {Loadmethod} "
                     "--probecutoff {probeCutoff} --samplecutoff {sampleCutoff} --detPcut {detPcut} --population {Population} "
                     "--filterXY {filterXY} --detPfile {detPfile} --analysis {Analysis} --normmethod {normMethod} "
                     "--cores {Cores} --diffexprprobeP {diffexprProbeP} --dmrmethod {dmrMethod} --minProbe {minProbe} "
                     "--diffexprregionP {diffexprRegionP} --dmrcores {dmrCores} --maxGap {maxGap} --adjPprobe {adjPprobe} "
                     "--bpSpan {bpSpan} --compareGroup '{compGroup}' --gseaMethod {gseaMethod} --gseaPvalue {gseaPvalue} "
                     "--controlGroup {controlGroup} --freqThreshold {freqThreshold} --sampleType {sampleType}".format(scriptpath = scriptpath, Import = Import,
                      Pdfile = Pdfile, matrixfile = matrixfile, chipType = chipType, Directory = Directory, Loadmethod = Loadmethod,
                      probeCutoff = probeCutoff, sampleCutoff = sampleCutoff, detPcut = detPcut, Population = Population,
                      filterXY = filterXY, detPfile = detPfile, Analysis = Analysis, normMethod = normMethod,
                      Cores = Cores, diffexprProbeP = diffexprProbeP, dmrMethod = dmrMethod, minProbe = minProbe,
                      diffexprRegionP = diffexprRegionP, dmrCores = dmrCores, maxGap = maxGap, adjPprobe = adjPprobe,
                      bpSpan = bpSpan, compGroup = compGroup, gseaMethod = gseaMethod, gseaPvalue = gseaPvalue,
                      controlGroup = controlGroup, freqThreshold = freqThreshold, sampleType = sampleType), shell=True)

def main(args, scriptpath):
    args = args
    scriptpath = scriptpath
    MethylChip(args, scriptpath)

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
This script is used for methyl-chip data analysis.
Usage:

python3 {scriptpath}/MethylChipPipeline.py -config config.yaml


'''.format(scriptpath = scriptpath))

    parse.add_argument('-config', '--yamlconfig', required = True, dest = "config", help = "Config file contain neccessary parameters", type = str, nargs = '?')

    args = parse.parse_args()

    main(args, scriptpath)
