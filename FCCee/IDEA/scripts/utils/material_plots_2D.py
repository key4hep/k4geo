from __future__ import print_function
import argparse

from plotstyle import FCCStyle
import math
import ROOT

def main():
    parser = argparse.ArgumentParser(description='Material Plotter')
    parser.add_argument('--fname', "-f", dest='fname', type=str, help="name of file to read")
    parser.add_argument('--angleMin', dest='angleMin', default=6, type=float, help="minimum eta/theta/cosTheta")
    parser.add_argument('--angleMax', dest='angleMax', default=6, type=float, help="maximum eta/theta/cosTheta")
    parser.add_argument('--angleDef', dest='angleDef', default="eta", type=str, help="angle definition to use: eta, theta or cosTheta, default: eta")
    parser.add_argument('--angleBinning', "-b", dest='angleBinning', default=0.05, type=float, help="eta/theta/cosTheta bin width")
    parser.add_argument('--nPhiBins', dest='nPhiBins', default=100, type=int, help="number of bins in phi")
    parser.add_argument('--x0max', "-x", dest='x0max', default=0.0, type=float, help="Max of x0")                                                                                                                                                
    args = parser.parse_args()

    ROOT.gStyle.SetNumberContours(100)

    f = ROOT.TFile.Open(args.fname, "read")
    tree = f.Get("materials")
    histDict = {}

    ROOT.gROOT.SetBatch(1)

    h_x0 = ROOT.TH2F("h_x0","h_x0", int((args.angleMax-args.angleMin)/args.angleBinning),args.angleMin,args.angleMax,args.nPhiBins,-math.pi,math.pi)
    h_lambda = ROOT.TH2F("h_lambda","h_lambda", int((args.angleMax-args.angleMin)/args.angleBinning),args.angleMin,args.angleMax,args.nPhiBins,-math.pi,math.pi)
    h_depth = ROOT.TH2F("h_depth","h_depth", int((args.angleMax-args.angleMin)/args.angleBinning),args.angleMin,args.angleMax,args.nPhiBins,-math.pi,math.pi)

    for angleBinning, entry in enumerate(tree):
        nMat = entry.nMaterials

        entry_x0, entry_lambda, entry_depth = 0.0, 0.0, 0.0
        for i in range(nMat):
            if entry.material.at(i) == "Air": continue

            entry_x0        += entry.nX0.at(i)*100.0
            entry_lambda    += entry.nLambda.at(i)
            entry_depth     += entry.matDepth.at(i)

        h_x0.Fill(tree.angle,tree.phi,entry_x0)
        h_lambda.Fill(tree.angle,tree.phi,entry_lambda)
        h_depth.Fill(tree.angle,tree.phi,entry_depth)

    # go through the 
    plots = ["x0", "lambda", "depth"]
    histograms = [h_x0, h_lambda, h_depth]
    axis_titles = ["Material budget x/X_{0} [%]", "Number of #lambda", "Material depth [cm]"]
    for i in range(len(plots)):
        cv = ROOT.TCanvas("","",800,600)
        cv.SetRightMargin(0.18)
        histograms[i].Draw("COLZ")

        if args.angleDef=="eta":
            title="#eta"
        elif args.angleDef=="theta":
            title="#theta [#circ]"
        elif args.angleDef=="cosTheta":
            title="cos(#theta)"
        histograms[i].GetXaxis().SetTitle(title)
        histograms[i].GetYaxis().SetTitle("#phi")

        histograms[i].GetZaxis().SetTitle(axis_titles[i])

        if args.x0max != 0.0 and plots[i]=="x0":
            histograms[i].SetMaximum(args.x0max)

        histograms[i].GetXaxis().SetRangeUser(args.angleMin, args.angleMax)

        ROOT.gStyle.SetPadRightMargin(0.5)
        cv.Print(plots[i] + ".pdf")
        cv.Print(plots[i] + ".png")
        cv.SaveAs(plots[i] + ".root")

if __name__ == "__main__":
    FCCStyle.initialize()
    main()
