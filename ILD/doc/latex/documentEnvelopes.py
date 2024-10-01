"""
Create Latex documents with drawings and tables
documenting envelope parameters

@author  F.Gaede, CERN/DESY
@version 1.0

"""

import sys
import copy
import cStringIO


# --- define the envelope parameters for every subdetector -----

envDict = {}
envDict["VXD"] = [
    "VXD_inner_radius",
    "VXD_outer_radius",
    "VXD_half_length",
    "VXD_cone_min_z",
    "VXD_cone_max_z",
    "VXD_inner_radius_1",
]

envDict["SIT"] = [
    "SIT_inner_radius",
    "SIT_outer_radius",
    "SIT_half_length",
    "SIT_outer_radius_1",
    "SIT_half_length_1",
]

envDict["TPC"] = ["TPC_inner_radius", "TPC_outer_radius", "TPC_half_length"]

envDict["FTD"] = [
    "FTD_inner_radius",
    "FTD_outer_radius",
    "FTD_half_length",
    "FTD_outer_radius_1",
    "FTD_outer_radius_2",
    "FTD_min_z_0",
    "FTD_min_z_1",
    "FTD_min_z_2",
    "FTD_cone_min_z",
    "FTD_cone_radius",
]

envDict["SET"] = ["SET_inner_radius", "SET_outer_radius", "SET_half_length"]


envDict["Ecal"] = [
    "Ecal_Hcal_symmetry",
    "Ecal_inner_radius",
    "Ecal_outer_radius",
    "Ecal_half_length",
    "Ecal_symmetry",
]

envDict["EcalEndcap"] = [
    "EcalEndcap_inner_radius",
    "EcalEndcap_outer_radius",
    "EcalEndcap_min_z",
    "EcalEndcap_max_z",
]

envDict["EcalEndcapRing"] = [
    "EcalEndcapRing_inner_radius",
    "EcalEndcapRing_outer_radius",
    "EcalEndcapRing_min_z",
    "EcalEndcapRing_max_z",
]

envDict["Hcal"] = [
    "Hcal_inner_radius",
    "Hcal_outer_radius",
    "Hcal_half_length",
    "Hcal_inner_symmetry",
]

envDict["HcalEndcap"] = [
    "EcalEndcap_symmetry",
    "HcalEndcap_inner_radius",
    "HcalEndcap_outer_radius",
    "HcalEndcap_min_z",
    "HcalEndcap_max_z",
]

envDict["HcalEndcapRing"] = [
    "HcalEndcapRing_inner_radius",
    "HcalEndcapRing_outer_radius",
    "HcalEndcapRing_min_z",
    "HcalEndcapRing_max_z",
    "HcalEndcapRing_symmetry",
]

envDict["Coil"] = ["Coil_inner_radius", "Coil_outer_radius", "Coil_half_length"]

envDict["Yoke"] = ["Yoke_inner_radius", "Yoke_outer_radius", "Yoke_half_length", "Yoke_symmetry"]

envDict["YokeEndcap"] = [
    "YokeEndcap_inner_radius",
    "YokeEndcap_outer_radius",
    "YokeEndcap_min_z",
    "YokeEndcap_max_z",
    "YokeEndcap_symmetry",
]

envDict["YokeEndcapPlug"] = [
    "YokeEndcapPlug_inner_radius",
    "YokeEndcapPlug_outer_radius",
    "YokeEndcapPlug_min_z",
    "YokeEndcapPlug_max_z",
    "YokeEndcapPlug_symmetry",
]

envDict["BeamCal"] = [
    "BeamCal_inner_radius",
    "BeamCal_outer_radius",
    "BeamCal_min_z",
    "BeamCal_max_z",
    "BeamCal_thickness",
    "BeamCal_tubeIncoming_radius",
]

envDict["LumiCal"] = [
    "LumiCal_inner_radius",
    "LumiCal_outer_radius",
    "LumiCal_min_z",
    "LumiCal_max_z",
    "LumiCal_thickness",
]

envDict["LHCal"] = [
    "LHCal_inner_radius",
    "LHCal_outer_radius",
    "LHCal_min_z",
    "LHCal_max_z",
    "LHCal_thickness",
]


# ----- define the envelope shape points in rz ------------

envRZDict = {}

envRZDict["VXD"] = [
    ("0", "VXD_inner_radius"),
    ("0", "VXD_outer_radius"),
    ("VXD_half_length", "VXD_outer_radius"),
    ("VXD_half_length", "VXD_inner_radius_1"),
    ("VXD_cone_max_z", "VXD_inner_radius_1"),
    ("VXD_cone_min_z", "VXD_inner_radius"),
    ("0", "VXD_inner_radius"),
]

envRZDict["SIT"] = [
    ("0", "SIT_inner_radius"),
    ("0", "SIT_outer_radius"),
    ("SIT_half_length", "SIT_outer_radius"),
    ("SIT_half_length", "SIT_outer_radius_1"),
    ("SIT_half_length_1", "SIT_outer_radius_1"),
    ("SIT_half_length_1", "SIT_inner_radius"),
    ("0", "SIT_inner_radius"),
]

envRZDict["FTD"] = [
    ("FTD_min_z_0", "FTD_inner_radius"),
    ("FTD_min_z_0", "FTD_outer_radius_1"),
    ("FTD_min_z_1", "FTD_outer_radius_1"),
    ("FTD_min_z_1", "FTD_outer_radius_2"),
    ("FTD_min_z_2", "FTD_outer_radius_2"),
    ("FTD_min_z_2", "FTD_outer_radius"),
    ("FTD_half_length", "FTD_outer_radius"),
    ("FTD_half_length", "FTD_cone_radius"),
    ("FTD_cone_min_z", "FTD_inner_radius"),
    ("FTD_min_z_0", "FTD_inner_radius"),
]

envRZDict["SET"] = [
    ("0", "SET_inner_radius"),
    ("0", "SET_outer_radius"),
    ("SET_half_length", "SET_outer_radius"),
    ("SET_half_length", "SET_inner_radius"),
    ("0", "SET_inner_radius"),
]

envRZDict["TPC"] = [
    ("0", "TPC_inner_radius"),
    ("0", "TPC_outer_radius"),
    ("TPC_half_length", "TPC_outer_radius"),
    ("TPC_half_length", "TPC_inner_radius"),
    ("0", "TPC_inner_radius"),
]

envRZDict["Ecal"] = [
    ("0", "Ecal_inner_radius"),
    ("0", "Ecal_outer_radius"),
    ("Ecal_half_length", "Ecal_outer_radius"),
    ("Ecal_half_length", "Ecal_inner_radius"),
    ("0", "Ecal_inner_radius"),
]

envRZDict["EcalEndcap"] = [
    ("EcalEndcap_min_z", "EcalEndcap_inner_radius"),
    ("EcalEndcap_min_z", "EcalEndcap_outer_radius"),
    ("EcalEndcap_max_z", "EcalEndcap_outer_radius"),
    ("EcalEndcap_max_z", "EcalEndcap_inner_radius"),
    ("EcalEndcap_min_z", "EcalEndcap_inner_radius"),
]

envRZDict["EcalEndcapRing"] = [
    ("EcalEndcapRing_min_z", "EcalEndcapRing_inner_radius"),
    ("EcalEndcapRing_min_z", "EcalEndcapRing_outer_radius"),
    ("EcalEndcapRing_max_z", "EcalEndcapRing_outer_radius"),
    ("EcalEndcapRing_max_z", "EcalEndcapRing_inner_radius"),
    ("EcalEndcapRing_min_z", "EcalEndcapRing_inner_radius"),
]

envRZDict["Hcal"] = [
    ("0", "Hcal_inner_radius"),
    ("0", "Hcal_outer_radius"),
    ("Hcal_half_length", "Hcal_outer_radius"),
    ("Hcal_half_length", "Hcal_inner_radius"),
    ("0", "Hcal_inner_radius"),
]


envRZDict["HcalEndcap"] = [
    ("HcalEndcap_min_z", "HcalEndcap_inner_radius"),
    ("HcalEndcap_min_z", "HcalEndcap_outer_radius"),
    ("HcalEndcap_max_z", "HcalEndcap_outer_radius"),
    ("HcalEndcap_max_z", "HcalEndcap_inner_radius"),
    ("HcalEndcap_min_z", "HcalEndcap_inner_radius"),
]

envRZDict["HcalEndcapRing"] = [
    ("HcalEndcapRing_min_z", "HcalEndcapRing_inner_radius"),
    ("HcalEndcapRing_min_z", "HcalEndcapRing_outer_radius"),
    ("HcalEndcapRing_max_z", "HcalEndcapRing_outer_radius"),
    ("HcalEndcapRing_max_z", "HcalEndcapRing_inner_radius"),
    ("HcalEndcapRing_min_z", "HcalEndcapRing_inner_radius"),
]


envRZDict["Yoke"] = [
    ("0", "Yoke_inner_radius"),
    ("0", "Yoke_outer_radius"),
    ("Yoke_half_length", "Yoke_outer_radius"),
    ("Yoke_half_length", "Yoke_inner_radius"),
    ("0", "Yoke_inner_radius"),
]


envRZDict["YokeEndcap"] = [
    ("YokeEndcap_min_z", "YokeEndcap_inner_radius"),
    ("YokeEndcap_min_z", "YokeEndcap_outer_radius"),
    ("YokeEndcap_max_z", "YokeEndcap_outer_radius"),
    ("YokeEndcap_max_z", "YokeEndcap_inner_radius"),
    ("YokeEndcap_min_z", "YokeEndcap_inner_radius"),
]

envRZDict["YokeEndcapPlug"] = [
    ("YokeEndcapPlug_min_z", "YokeEndcapPlug_inner_radius"),
    ("YokeEndcapPlug_min_z", "YokeEndcapPlug_outer_radius"),
    ("YokeEndcapPlug_max_z", "YokeEndcapPlug_outer_radius"),
    ("YokeEndcapPlug_max_z", "YokeEndcapPlug_inner_radius"),
    ("YokeEndcapPlug_min_z", "YokeEndcapPlug_inner_radius"),
]

envRZDict["Coil"] = [
    ("0", "Coil_inner_radius"),
    ("0", "Coil_outer_radius"),
    ("Coil_half_length", "Coil_outer_radius"),
    ("Coil_half_length", "Coil_inner_radius"),
    ("0", "Coil_inner_radius"),
]


envRZDict["BeamCal"] = [
    ("BeamCal_min_z", "BeamCal_inner_radius"),
    ("BeamCal_min_z", "BeamCal_outer_radius"),
    ("BeamCal_max_z", "BeamCal_outer_radius"),
    ("BeamCal_max_z", "BeamCal_inner_radius"),
    ("BeamCal_min_z", "BeamCal_inner_radius"),
]

envRZDict["LumiCal"] = [
    ("LumiCal_min_z", "LumiCal_inner_radius"),
    ("LumiCal_min_z", "LumiCal_outer_radius"),
    ("LumiCal_max_z", "LumiCal_outer_radius"),
    ("LumiCal_max_z", "LumiCal_inner_radius"),
    ("LumiCal_min_z", "LumiCal_inner_radius"),
]

envRZDict["LHCal"] = [
    ("LHCal_min_z", "LHCal_inner_radius"),
    ("LHCal_min_z", "LHCal_outer_radius"),
    ("LHCal_max_z", "LHCal_outer_radius"),
    ("LHCal_max_z", "LHCal_inner_radius"),
    ("LHCal_min_z", "LHCal_inner_radius"),
]


# -----------------------------------------------

try:
    dictFile = sys.argv[1]

except IndexError:
    print(" usage:  python documentEnvelopes.py pyDict.py ")
    print("    pyDict.py : python file with a data dictionary (created with extractParameters)")
    print()
    sys.exit(1)


# ------ read dictionary 'values' from file
execfile(dictFile)
values["0"] = 0


# -----------------------------------------------
def run():
    writeTexFile("VXD", "_rz_envelope", getRZEnvCmds, 20)
    writeTexFile("SIT", "_rz_envelope", getRZEnvCmds, 20)
    writeTexFile("FTD", "_rz_envelope", getRZEnvCmds, 20)
    writeTexFile("TPC", "_rz_envelope", getRZEnvCmds, 20)
    writeTexFile("SET", "_rz_envelope", getRZEnvCmds, 20)

    writeTexFile("ILD", "_rz_quadrant", getILDRZQuadrantCmds, 20)

    writeILDEnvTable("ILD_enevelope_table.tex")


# -----------------------------------------------
def writeILDEnvTable(file):
    fn = "./ILD_envelopeTable.tex"
    of = open(fn, "w")
    cmds = []

    cmds.extend(getDocHeaderCmds("article", ["multirow"]))

    dets = [
        "VXD",
        "FTD",
        "SIT",
        "TPC",
        "SET",
        "Ecal",
        "EcalEndcap",
        "EcalEndcapRing",
        "Hcal",
        "HcalEndcap",
        "HcalEndcapRing",
        "Coil",
        "Yoke",
        "YokeEndcap",
        "YokeEndcapPlug",
        "BeamCal",
        "LHCal",
        "LumiCal",
    ]

    cmds.extend(getEnvelopeTableCmds(dets, "\\large{Envelope parameters for ILD\_o1\_v05}"))

    cmds.extend(getDocFooterCmds())

    for cmd in cmds:
        print(cmd, file=of)

    of.close()


# -----------------------------------------------
def getEnvelopeTableCmds(dets, title):
    cmds = []
    cmds.append("\\begin{tabular}{|l | c | c | c | l r |}")
    cmds.append("\\hline")

    if len(title):
        cmds.append("\\multicolumn{6}{|c|}{} \\\\")
        cmds.append("\\multicolumn{6}{|c|}{" + title + "} \\\\")
        cmds.append("\\multicolumn{6}{|c|}{} \\\\")
        cmds.append("\\hline")

    cmds.append(
        " detector & inner radius & outer radius & half length  & \multicolumn{2}{c|}{additional parameters} \\\\"
    )
    cmds.append("          &              &              & min z, max z &          &        \\\\")
    cmds.append("\\hline")

    for d in dets:
        cmds.extend(getTableLinesCmds(d))

    cmds.append("\\end{tabular}")
    return cmds


# -----------------------------------------------
def getTableLinesCmds(det):
    cmds = []
    params = copy.deepcopy(envDict[det])

    ri = det + "_inner_radius"
    ro = det + "_outer_radius"
    hl = det + "_half_length"
    zs = det + "_min_z"
    ze = det + "_max_z"

    line = det + " & "

    if ri in params:
        line += ("%.1f" % values[ri]) + " & "
        params.remove(ri)
    else:
        line += " -  & "

    if ro in params:
        line += ("%.1f" % values[ro]) + " & "
        params.remove(ro)
    else:
        line += " -  & "

    if hl in params:
        line += ("%.1f" % values[hl]) + " & "
        params.remove(hl)
    else:
        line += ("%.1f" % values[zs]) + ", " + ("%.1f" % values[ze]) + " & "
        params.remove(zs)
        params.remove(ze)

    # --- first extra parameter - if any
    if len(params) > 0:
        p = params[0]
        line += "\\small{\\verb#" + p + "#} & " + ("%.1f" % values[p])
        params.remove(p)
    else:
        line += " & "

    line += "  \\\\ "
    cmds.append(line)

    # --- other extra parameters - if any - need extra line

    while len(params) > 0:
        line = " & & & & "
        p = params[0]
        line += "\\small{\\verb#" + p + "#} & " + ("%.1f" % values[p])
        params.remove(p)
        line += "  \\\\ "
        cmds.append(line)

    cmds.append("\\hline")

    return cmds


# -----------------------------------------------
def getILDRZQuadrantCmds(det, width):
    cmds = []

    cmds.extend(getColorCmds())

    cmds.append("\\begin{tikzpicture}")

    scale = 0.01
    # ---------------------------------------------------

    cmds.append(lineOStr("[fill=VXDcol]", getEnvPoints("VXD", scale)))
    cmds.append(lineOStr("[fill=SITcol]", getEnvPoints("SIT", scale)))
    cmds.append(lineOStr("[fill=FTDcol]", getEnvPoints("FTD", scale)))
    cmds.append(lineOStr("[fill=TPCcol]", getEnvPoints("TPC", scale)))
    cmds.append(lineOStr("[fill=ECALcol]", getEnvPoints("Ecal", scale)))
    cmds.append(lineOStr("[fill=ECALcol]", getEnvPoints("EcalEndcap", scale)))
    cmds.append(lineOStr("[fill=ECALcol]", getEnvPoints("EcalEndcapRing", scale)))
    cmds.append(lineOStr("[fill=HCALcol]", getEnvPoints("Hcal", scale)))
    cmds.append(lineOStr("[fill=HCALcol]", getEnvPoints("HcalEndcap", scale)))
    cmds.append(lineOStr("[fill=HCALcol]", getEnvPoints("HcalEndcapRing", scale)))
    cmds.append(lineOStr("[fill=YOKEcol]", getEnvPoints("Yoke", scale)))
    cmds.append(lineOStr("[fill=YOKEcol]", getEnvPoints("YokeEndcap", scale)))
    cmds.append(lineOStr("[fill=YOKEcol]", getEnvPoints("YokeEndcapPlug", scale)))
    cmds.append(lineOStr("[fill=COILcol]", getEnvPoints("Coil", scale)))
    cmds.append(lineOStr("[fill=SITcol]", getEnvPoints("BeamCal", scale)))
    cmds.append(lineOStr("[fill=SITcol]", getEnvPoints("LumiCal", scale)))
    cmds.append(lineOStr("[fill=SITcol]", getEnvPoints("LHCal", scale)))

    cmds.append("\\end{tikzpicture}")
    return cmds


# -----------------------------------------------
def getEnvPoints(det, scale):
    points = []
    env = envRZDict[det]
    for ep in env:
        p = (scale * values[ep[0]], scale * values[ep[1]])
        print(det, " point: ", p)
        points.append(p)
    return points


# -----------------------------------------------
def writeTexFile(det, fileExt, method, width):
    fn = "./figs/" + det + fileExt + ".tex"

    of = open(fn, "w")

    cmds = []

    cmds.extend(getDocHeaderCmds("standalone", ["tikz", "graphicx"]))

    cmds.extend(method(det, width))

    cmds.extend(getDocFooterCmds())

    for cmd in cmds:
        print(cmd, file=of)

    of.close()


# -----------------------------------------------
def lineStr(tl, opt=""):
    o = cStringIO.StringIO()

    print("\draw ", opt, file=o)

    i = 0
    for t in tl:
        if i > 0:
            print(" -- ", file=o)
        print("(", t[0], ",", t[1], ") ", file=o)
        i += 1
    print(";", file=o)
    str = o.getvalue()
    o.close()
    return str


# -----------------------------------------------


def lineOStr(opt, tl):
    return lineStr(tl, opt)


# -----------------------------------------------


def getDocHeaderCmds(docClass, packages):
    cmds = []

    cmds.append("\\documentclass[a4]{" + docClass + "}")
    if docClass == "standalone":
        cmds.append("\\standaloneconfig{border=20pt}")

    for p in packages:
        cmds.append("\\usepackage{" + p + "}")

    cmds.append("\\usepackage[a4paper,top=3cm, bottom=2.5cm, left=1cm, right=2cm ]{geometry}")

    cmds.append("\\begin{document}")

    return cmds


# -----------------------------------------------
def getColorCmds():
    cmds = []
    cmds.append("\\definecolor{VXDcol}{RGB}{255,255,255}")
    cmds.append("\\definecolor{SITcol}{RGB}{221,221,221}")
    cmds.append("\\definecolor{SETcol}{RGB}{221,221,221}")
    cmds.append("\\definecolor{TPCcol}{RGB}{245,243,0}")
    cmds.append("\\definecolor{ECALcol}{RGB}{123,243,0}")
    cmds.append("\\definecolor{HCALcol}{RGB}{196,194,49}")
    cmds.append("\\definecolor{YOKEcol}{RGB}{24,194,196}")
    cmds.append("\\definecolor{COILcol}{RGB}{73,73,221}")
    cmds.append("\\definecolor{FTDcol}{RGB}{101,28,147}")
    cmds.append("\\definecolor{FCALcol}{RGB}{171,170,171}")
    return cmds


def getDocFooterCmds():
    cmds = []
    cmds.append("\\end{document}")
    return cmds


# -----------------------------------------------
def fixstr(aStr):
    l = len(aStr)
    s = ""
    for i in range(0, l + 2):
        s = s + " "
    s = s + aStr
    return s


# -----------------------------------------------


def getRZEnvCmds(det, width):
    cmds = []
    envPoints = envRZDict[det]

    vals = values

    cmds.append("\\begin{tikzpicture}")

    xmaxOrg = -1e99
    # ---- compute the scale such that  ymax == width
    for ep in envPoints:
        if vals[ep[0]] > xmaxOrg:
            xmaxOrg = vals[ep[0]]

    scale = width / xmaxOrg

    # ---------------------------------------------------

    points = []
    xvals = []
    yvals = []
    xmax, ymax = -1.0e99, -1e99

    for ep in envPoints:
        p = (scale * vals[ep[0]], scale * vals[ep[1]])
        points.append(p)

        x, y = p
        if x > xmax:
            xmax = x
        if y > ymax:
            ymax = y

        if ep[0] not in xvals:
            xvals.append(ep[0])
            p0 = (scale * vals[ep[0]], 0)
            cmds.append(lineOStr("[dashed]", [p, p0]))
            cmds.append(
                "\\node [rotate=-90] at ("
                + str(p0[0])
                + ","
                + str(p0[1])
                + ") {{\\verb#"
                + fixstr(ep[0])
                + "#}};"
            )

        if ep[1] not in yvals:
            yvals.append(ep[1])
            p1 = (0, scale * vals[ep[1]])
            cmds.append(lineOStr("[dashed]", [p, p1]))
            cmds.append(
                "\\node [left] at ("
                + str(p1[0])
                + ","
                + str(p1[1])
                + ") {{\\verb#"
                + str(ep[1])
                + "#}};"
            )

    cmds.append(lineOStr("[ultra thick]", points))

    cmds.append(lineOStr("[<->,thick]", ((0, 1.25 * ymax), (0, 0), (1.25 * xmax, 0))))

    cmds.append("\\end{tikzpicture}")

    return cmds


# -----------------------------------------------


# -----------------------------------------------
if __name__ == "__main__":
    run()
# -----------------------------------------------
