#!/usr/bin/env python3

"""Check that a loaded DD4hep geometry has the expected detector name in its header."""

import sys
import argparse

import dd4hep


def main():
    parser = argparse.ArgumentParser(
        description="Load a DD4hep compact geometry and verify its header name."
    )
    parser.add_argument("-c", "--compact", help="Compact file location", required=True, type=str)
    parser.add_argument("-n", "--name", help="Expected detector name", required=True, type=str)
    args = parser.parse_args()

    description = dd4hep.Detector.make_unique("detector")
    description.fromXML(args.compact)

    actual_name = description.header().name()
    if actual_name != args.name:
        print(f"ERROR: Detector name mismatch: expected '{args.name}', got '{actual_name}'")
        sys.exit(1)

    print(f"OK: Detector name '{actual_name}' matches expected name '{args.name}'")


if __name__ == "__main__":
    main()
