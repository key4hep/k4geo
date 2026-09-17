#!/usr/bin/env python
# coding: utf-8

from dd4hep import Detector
import graphviz
import re
from collections import Counter

from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    "--compactFile", help="DD4hep compact description xml", required=True
)
parser.add_argument(
    "--maxDepth",
    help="Maximum traversal depth of the detector tree",
    default=10,
    type=int,
)
parser.add_argument(
    "--maxEdges", help="Maximum number of edges per connection", default=5, type=int
)
parser.add_argument(
    "--drawList",
    action="extend",
    nargs="+",
    metavar=("Vertex", "InnerTrackers"),
    help="Names of the immediate children of the start element that should be drawn.",
    default=[],
)
parser.add_argument(
    "--startNode", help="Start of the traversal", default="world"
)

args = parser.parse_args()

theDetector = Detector.getInstance()
theDetector.fromXML(args.compactFile)

# take part between last / and before the .xml
detector_name = args.compactFile.split("/")[-1].split(".")[0]

dot = graphviz.Digraph("detector_name")

start = theDetector.detector(args.startNode)
# start = theDetector.world()

dot.node(start.name(), start.name())


def process_name(raw_name):
    name = re.sub(r"\d+", "X", raw_name)
    return name


edge_counter = Counter()


def add_children(detElement, depth=0, children=[]):
    depth += 1
    if not children:
        children = detElement.children()
    for raw_name, child in children:
        name = process_name(raw_name)
        parent_name = process_name(detElement.name())
        sens = child.volume().isSensitive()
        sens_lbl = "sensitive" if sens else ""
        dot.node(name, name, {"xlabel": sens_lbl})
        # dot.node(name, name)
        if edge_counter[(parent_name, name)] < args.maxEdges:
            dot.edge(parent_name, name)
            edge_counter[(parent_name, name)] += 1
        if depth < args.maxDepth:
            add_children(child, depth)


start_children = [(name, start.child(name)) for name in args.drawList]
add_children(start, children=start_children)

dot.render(f"{detector_name}_{args.startNode}_{'-'.join(args.drawList)}.gv")
