#!/usr/bin/env python3

### Script to generate ROOT file with detector geometry from detector xml file,
### from https://fccsw.web.cern.ch/fccsw/tutorials/static/python/dd4hep2root
### usage described in FCC software tutorial:
### https://hep-fcc.github.io/fcc-tutorials/master/full-detector-simulations/Visualization/Visualization.html#detector-geometry

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Convert detector")
    parser.add_argument(
        "-c", "--compact", help="Compact file location(s)", required=True, type=str, nargs="+"
    )
    parser.add_argument(
        "-o", "--out", help="Converted file path", default="detector.root", type=str
    )
    args = parser.parse_args()

    convert(args.compact, args.out)


def convert(compact_files, out_path):
    print("INFO: Converting following compact file(s):")
    for cfile in compact_files:
        print("      " + cfile)

    import ROOT

    ROOT.gSystem.Load("libDDCore")
    description = ROOT.dd4hep.Detector.getInstance()
    for cfile in compact_files:
        description.fromXML(cfile)

    ROOT.gGeoManager.SetVisLevel(9)
    ROOT.gGeoManager.SetVisOption(0)
    ROOT.gGeoManager.Export(out_path)


if __name__ == "__main__":
    main()
