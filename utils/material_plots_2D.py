"""
This script must be called with python: 'python material_scan_2D.py --{argument} {value}'.
The output files are saved in data/{outputDir}/name.suffix.
If no outputDir is specified, it will be data/plots/name.suffix.
"""

from __future__ import print_function

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.append(os.path.expandvars("$FCCSW") + "/Examples/scripts")
from plotstyle import FCCStyle
import ROOT


def create_histogram(
    name_and_title: str,
    angle_min: float,
    angle_max: float,
    angle_binning: float,
    n_phi_bins: int,
) -> ROOT.TH2F:
    num_bins = int((angle_max - angle_min) / angle_binning)
    return ROOT.TH2F(
        name_and_title,
        name_and_title,
        num_bins,
        angle_min,
        angle_max,
        n_phi_bins,
        -math.pi,
        math.pi,
    )


def main():
    parser = argparse.ArgumentParser(description="Material Plotter")
    parser.add_argument(
        "--inputFile", "--fname", "-f", type=str, help="relative path to the input file"
    )
    parser.add_argument(
        "--angleMin", dest="angleMin", default=6, type=float, help="minimum eta/theta/cosTheta"
    )
    parser.add_argument(
        "--angleMax", dest="angleMax", default=6, type=float, help="maximum eta/theta/cosTheta"
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
    parser.add_argument("--nPhiBins", default=100, type=int, help="Number of bins in phi")
    parser.add_argument("--x0max", "-x", default=0.0, type=float, help="Max of x0")
    parser.add_argument(
        "--outputDir",
        "-o",
        type=str,
        default="plots",
        help="Directory to store output files in",
    )
    parser.add_argument(
        "--ignoreMats",
        "-i",
        dest="ignoreMats",
        nargs="+",
        default=[],
        help="List of materials that should be ignored",
    )

    args = parser.parse_args()

    output_dir = Path("data") / args.outputDir
    output_dir.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

    ROOT.gStyle.SetNumberContours(100)

    f = ROOT.TFile.Open(os.fspath(Path(args.inputFile).with_suffix(".root")), "read")
    tree = f.Get("materials")

    ROOT.gROOT.SetBatch(1)

    h_x0 = create_histogram("h_x0", args.angleMin, args.angleMax, args.angleBinning, args.nPhiBins)
    h_lambda = create_histogram(
        "h_lambda", args.angleMin, args.angleMax, args.angleBinning, args.nPhiBins
    )
    h_depth = create_histogram(
        "h_depth", args.angleMin, args.angleMax, args.angleBinning, args.nPhiBins
    )

    for angleBinning, entry in enumerate(tree):
        nMat = entry.nMaterials

        entry_x0, entry_lambda, entry_depth = 0.0, 0.0, 0.0
        for i in range(nMat):
            # Ignore certain materials if specified
            if entry.material.at(i) in args.ignoreMats:
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
    axis_titles = ["Material budget x/X_{0} [%]", "Number of #lambda", "Material depth [cm]"]
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

        if args.x0max != 0.0 and plot == "x0":
            histograms[i].SetMaximum(args.x0max)

        histograms[i].GetXaxis().SetRangeUser(args.angleMin, args.angleMax)

        ROOT.gStyle.SetPadRightMargin(0.5)
        output_path = output_dir / plot
        cv.Print(os.fspath(output_path.with_suffix(".pdf")))
        cv.Print(os.fspath(output_path.with_suffix(".png")))
        cv.SaveAs(os.fspath(output_path.with_suffix(".root")))


if __name__ == "__main__":
    FCCStyle.initialize()
    main()
