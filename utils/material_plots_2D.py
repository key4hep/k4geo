"""
This script must be called with python: 'python material_scan_2D.py --{argument} {value}'.
The output files are saved in data/{outputDir}/name.suffix.
If no outputDir is specified, it will be data/plots/name.suffix.
"""

from __future__ import print_function

import argparse
import math
import sys
from os import fspath
from os.path import expandvars
from pathlib import Path

import ROOT

sys.path.append(expandvars("$FCCSW") + "/Examples/scripts")
from plotstyle import FCCStyle


def main():
    parser = argparse.ArgumentParser(description="Material Plotter")
    parser.add_argument(
        "--inputFile", "-f", type=str, help="relative path to the input file"
    )
    parser.add_argument(
        "--angleMin",
        default=6,
        type=float,
        help="Minimum eta/theta/cosTheta",
    )
    parser.add_argument(
        "--angleMax",
        default=6,
        type=float,
        help="Maximum eta/theta/cosTheta",
    )
    parser.add_argument(
        "--angleDef",
        default="eta",
        choices=["eta", "theta", "cosTheta", "thetaRad"],
        type=str,
        help="Angle definition to use: eta, theta, thetaRad or cosTheta; default: eta",
    )
    parser.add_argument(
        "--angleBinning",
        "-b",
        default=0.05,
        type=float,
        help="Eta/theta/cosTheta bin width",
    )
    parser.add_argument(
        "--nPhiBins",
        default=100,
        type=int,
        help="Number of bins in phi",
    )
    parser.add_argument("--x0Max", "-x", default=0.0, type=float, help="Max of x0")
    parser.add_argument(
        "--outputDir",
        "-o",
        type=str,
        default="plots",
        help="Directory to store output files in",
    )
    args = parser.parse_args()

    output_dir = Path("data") / args.outputDir
    output_dir.mkdir(
        parents=True, exist_ok=True
    )  # Create the directory if it doesn't exist

    ROOT.gStyle.SetNumberContours(100)

    f = ROOT.TFile.Open(fspath(Path(args.inputFile).with_suffix(".root")), "read")
    tree = f.Get("materials")

    ROOT.gROOT.SetBatch(1)

    h_x0 = ROOT.TH2F(
        "h_x0",
        "h_x0",
        int((args.angleMax - args.angleMin) / args.angleBinning),
        args.angleMin,
        args.angleMax,
        args.nPhiBins,
        -math.pi,
        math.pi,
    )
    h_lambda = ROOT.TH2F(
        "h_lambda",
        "h_lambda",
        int((args.angleMax - args.angleMin) / args.angleBinning),
        args.angleMin,
        args.angleMax,
        args.nPhiBins,
        -math.pi,
        math.pi,
    )
    h_depth = ROOT.TH2F(
        "h_depth",
        "h_depth",
        int((args.angleMax - args.angleMin) / args.angleBinning),
        args.angleMin,
        args.angleMax,
        args.nPhiBins,
        -math.pi,
        math.pi,
    )

    for angleBinning, entry in enumerate(tree):
        nMat = entry.nMaterials

        entry_x0, entry_lambda, entry_depth = 0.0, 0.0, 0.0
        for i in range(nMat):
            if entry.material.at(i) == "Air":
                continue

            entry_x0 += entry.nX0.at(i) * 100.0
            entry_lambda += entry.nLambda.at(i)
            entry_depth += entry.matDepth.at(i)

        h_x0.Fill(tree.angle, tree.phi, entry_x0)
        h_lambda.Fill(tree.angle, tree.phi, entry_lambda)
        h_depth.Fill(tree.angle, tree.phi, entry_depth)

    # go through the plots
    plots = ["x0", "lambda", "depth"]
    histograms = [h_x0, h_lambda, h_depth]
    axis_titles = [
        "Material budget x/X_{0} [%]",
        "Number of #lambda",
        "Material depth [cm]",
    ]
    for i, plot in enumerate(plots):
        cv = ROOT.TCanvas("", "", 800, 600)
        cv.SetRightMargin(0.18)
        histograms[i].Draw("COLZ")

        if args.angleDef == "eta":
            title = "#eta"
        elif args.angleDef == "theta":
            title = "#theta [#circ]"
        elif args.angleDef == "thetaRad":
            title = "#theta [rad]"
        elif args.angleDef == "cosTheta":
            title = "cos(#theta)"
        histograms[i].GetXaxis().SetTitle(title)
        histograms[i].GetYaxis().SetTitle("#phi")

        histograms[i].GetZaxis().SetTitle(axis_titles[i])

        if args.x0Max != 0.0 and plot == "x0":
            histograms[i].SetMaximum(args.x0Max)

        histograms[i].GetXaxis().SetRangeUser(args.angleMin, args.angleMax)

        ROOT.gStyle.SetPadRightMargin(0.5)
        output_path = output_dir / plot
        cv.Print(fspath(output_path.with_suffix(".pdf")))
        cv.Print(fspath(output_path.with_suffix(".png")))
        cv.SaveAs(fspath(output_path.with_suffix(".root")))


if __name__ == "__main__":
    FCCStyle.initialize()
    main()
