from __future__ import print_function
import argparse

import sys, os

sys.path.append(os.path.expandvars("$FCCSW") + "/Examples/scripts")
from plotstyle import FCCStyle

import ROOT


def main():
    parser = argparse.ArgumentParser(description="Material Plotter")
    parser.add_argument("--fname", "-f", dest="fname", type=str, help="name of file to read")
    parser.add_argument("--outputDir", "-o", dest="outputDir", default=".", type=str, help="Output directory for plots")
    parser.add_argument(
        "--angleMin", dest="angleMin", default=-6, type=float, help="minimum eta/theta/cosTheta"
    )
    parser.add_argument(
        "--angleMax", dest="angleMax", default=6, type=float, help="maximum eta/theta/cosTheta"
    )
    parser.add_argument(
        "--angleDef",
        dest="angleDef",
        default="eta",
        type=str,
        help="angle definition to use: eta, theta, cosTheta or thetaRad, default: eta",
    )
    parser.add_argument(
        "--angleBinning",
        "-b",
        dest="angleBinning",
        default=0.05,
        type=float,
        help="eta/theta/cosTheta/thetaRad bin width",
    )
    parser.add_argument("--x0max", "-x", dest="x0max", default=0.0, type=float, help="Max of x0")
    parser.add_argument(
        "--removeMatsSubstrings",
        dest="removeMatsSubstrings",
        nargs="+",
        default=[],
        help="Substrings to be removed from materials strings (e.g. '66D' for reduced density materials)",
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

    f = ROOT.TFile.Open(args.fname, "read")
    tree = f.Get("materials")
    histDict = {}

    ROOT.gROOT.SetBatch(1)

    # go through the eta/theta/cosTheta/thetaRad bins and fill the histograms in the histDict, skipping air
    # keys in the histDict are the material names
    for angleBinning, entry in enumerate(tree):
        nMat = entry.nMaterials
        for i in range(nMat):
            material = entry.material.at(i)

            # Removing substrings from materials
            for substring in args.removeMatsSubstrings:
                material = material.replace(substring, "")

            # Ignore certain materials if specified
            if material in args.ignoreMats:
                continue

            if material not in histDict.keys():
                histDict[material] = {
                    "x0": ROOT.TH1F(
                        "",
                        "",
                        (int)((args.angleMax - args.angleMin) / args.angleBinning),
                        args.angleMin,
                        args.angleMax,
                    ),
                    "lambda": ROOT.TH1F(
                        "",
                        "",
                        (int)((args.angleMax - args.angleMin) / args.angleBinning),
                        args.angleMin,
                        args.angleMax,
                    ),
                    "depth": ROOT.TH1F(
                        "",
                        "",
                        (int)((args.angleMax - args.angleMin) / args.angleBinning),
                        args.angleMin,
                        args.angleMax,
                    ),
                }
            hs = histDict[material]
            hs["x0"].SetBinContent(
                angleBinning + 1,
                hs["x0"].GetBinContent(angleBinning + 1) + entry.nX0.at(i) * 100.0,
            )
            hs["lambda"].SetBinContent(
                angleBinning + 1,
                hs["lambda"].GetBinContent(angleBinning + 1) + entry.nLambda.at(i),
            )
            hs["depth"].SetBinContent(
                angleBinning + 1,
                hs["depth"].GetBinContent(angleBinning + 1) + entry.matDepth.at(i),
            )

    axis_titles = ["Material budget x/X_{0} [%] ", "Number of #lambda", "Material depth [cm]"]

    # This loop does the drawing, sets the style and saves the pdf files
    for plot, title in zip(["x0", "lambda", "depth"], axis_titles):
        if args.angleDef == "eta":
            xtitle = "#eta"
            legend = ROOT.TLegend(0.2, 0.6, 0.5, 0.94)
        elif args.angleDef == "theta":
            xtitle = "#theta"
            legend = ROOT.TLegend(0.5, 0.6, 0.8, 0.94)
        elif args.angleDef == "thetaRad":
            xtitle = "#theta [rad]"
            legend = ROOT.TLegend(0.5, 0.6, 0.8, 0.94)
        elif args.angleDef == "cosTheta":
            xtitle = "cos(#theta)"
            legend = ROOT.TLegend(0.2, 0.6, 0.5, 0.94)

        legend.SetLineColor(0)
        ths = ROOT.THStack()

        # sort by list of values: To do so set reorder = True and define your desired order of the materials in the "order" list
        reorder = False
        histDict_ordered = {}
        if reorder:
            order = [
                "Silicon",
                "CarbonFiber",
                "CarbonFleece",
                "Rohacell",
                "Aluminum",
                "GlueEcobond45",
                "Kapton",
                "Water",
            ]
            ordered_list = sorted(histDict.items(), key=lambda pair: order.index(pair[0]))

            for key, value in ordered_list:
                histDict_ordered[key] = value
            print(histDict_ordered)
        else:
            histDict_ordered = histDict

        # Make the plots
        for i, material in enumerate(histDict_ordered.keys()):
            linecolor = 1
            fillcolor = FCCStyle.fillcolors[i if i < 7 else 0]

            # If you want to map a material to a specific color, do that here
            match material:
                case "CarbonFiber":
                    fillcolor = FCCStyle.fillcolors[0]
                case "CarbonFoam":
                    fillcolor = FCCStyle.fillcolors[0]
                case "CarbonFleece":
                    fillcolor = ROOT.kBlack
                case "Rohacell":
                    fillcolor = FCCStyle.fillcolors[4]
                case "Silicon":
                    fillcolor = FCCStyle.fillcolors[2]
                case "Aluminum":
                    fillcolor = FCCStyle.fillcolors[1]
                case "Kapton":
                    fillcolor = FCCStyle.fillcolors[3]
                case "GlueEcobond45":
                    fillcolor = FCCStyle.fillcolors[6]
                case "Water":
                    fillcolor = FCCStyle.fillcolors[5]
                case "PCB":
                    fillcolor = ROOT.kGreen

            histDict_ordered[material][plot].SetLineColor(linecolor)
            histDict_ordered[material][plot].SetFillColor(fillcolor)
            histDict_ordered[material][plot].SetLineWidth(1)
            histDict_ordered[material][plot].SetFillStyle(1001)

            ths.Add(histDict_ordered[material][plot])

        for i, material in enumerate(reversed(histDict_ordered.keys())):
            legend.AddEntry(histDict_ordered[material][plot], material, "f")

        cv = ROOT.TCanvas()
        ths.Draw()
        ths.GetXaxis().SetTitle(xtitle)
        ths.GetYaxis().SetTitle(title)

        legend.SetTextSize(0.04)
        legend.Draw()
        if args.x0max != 0.0 and plot == "x0":
            ths.SetMaximum(args.x0max)

        ths.GetXaxis().SetRangeUser(args.angleMin, args.angleMax)
        output_path = os.path.join(args.outputDir, plot)
        cv.Print(output_path + ".pdf")
        cv.Print(output_path + ".png")
        cv.SaveAs(output_path + ".root")


if __name__ == "__main__":
    FCCStyle.initialize()
    main()
