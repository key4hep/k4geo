#!/usr/bin/env python3

"""Check that every volume placed in the world lies inside the world box."""

import sys
import argparse
from array import array

import ROOT
import dd4hep


def extent(node, matrix, ext):
    """Update ext with the largest |x|, |y|, |z| of the node's shape vertices in world coordinates"""
    matrix = ROOT.TGeoHMatrix(matrix)
    matrix.Multiply(node.GetMatrix())
    volume = node.GetVolume()
    if volume.IsAssembly():  # an assembly has no shape of its own: look at its daughters
        for i in range(volume.GetNdaughters()):
            extent(volume.GetNode(i), matrix, ext)
        return
    shape = volume.GetShape()
    points = array("d", [0.0] * (3 * shape.GetNmeshVertices()))
    shape.SetPoints(points)
    local, master = array("d", [0.0] * 3), array("d", [0.0] * 3)
    for i in range(0, len(points), 3):
        local[0], local[1], local[2] = points[i], points[i + 1], points[i + 2]
        matrix.LocalToMaster(local, master)
        for j in range(3):
            ext[j] = max(ext[j], abs(master[j]))


def main():
    parser = argparse.ArgumentParser(
        description="Load a DD4hep compact geometry and check that nothing extrudes from the world."
    )
    parser.add_argument("-c", "--compact", help="Compact file location", required=True, type=str)
    args = parser.parse_args()

    # getInstance, not make_unique: some segmentations read their metadata from the global instance
    description = dd4hep.Detector.getInstance()
    description.fromXML(args.compact)
    description.manager().SetNsegments(360)  # vertices of curved shapes every degree

    world = description.manager().GetTopVolume()
    world_half = [world.GetShape().GetDX(), world.GetShape().GetDY(), world.GetShape().GetDZ()]

    failed = False
    for i in range(world.GetNdaughters()):
        node = world.GetNode(i)
        ext = [0.0, 0.0, 0.0]
        extent(node, ROOT.TGeoHMatrix(), ext)
        excess = max(ext[j] - world_half[j] for j in range(3))
        if excess > 0:
            print(f"ERROR: '{node.GetName()}' extends {excess * 10:.1f} mm beyond the world")
            failed = True

    if failed:
        sys.exit(1)
    print(f"OK: all {world.GetNdaughters()} volumes lie inside the world")


if __name__ == "__main__":
    main()
