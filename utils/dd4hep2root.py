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

    predefinedColors = dict((color.GetNumber(), color) for color in ROOT.gROOT.GetListOfColors() if color)
    cachedColors = dict((num, num) for num in predefinedColors)

    def mapColor(colorNumber):
        if colorNumber not in cachedColors:
            color = ROOT.gROOT.GetColor(colorNumber)

            # Guard clause in case GetColor() doesn't return a valid color
            if not color:
                return colorNumber

            r, g, b = color.GetRed(), color.GetGreen(), color.GetBlue()

            # Score every predefined ROOT color by square distance
            mostSimilar = min(
                predefinedColors,
                key=lambda num: (
                    (predefinedColors[num].GetRed() - r) ** 2
                    + (predefinedColors[num].GetGreen() - g) ** 2
                    + (predefinedColors[num].GetBlue() - b) ** 2
                ),
            )

            cachedColors[colorNumber] = mostSimilar

        return cachedColors[colorNumber]

    # Map volumes to the closest predefined ROOT colors
    for volume in ROOT.gGeoManager.GetListOfVolumes():
        volume.SetLineColor(mapColor(volume.GetLineColor()))
        volume.SetFillColor(mapColor(volume.GetFillColor()))

    ROOT.gGeoManager.SetVisLevel(9)
    ROOT.gGeoManager.SetVisOption(0)
    ROOT.gGeoManager.Export(out_path)


if __name__ == "__main__":
    main()
