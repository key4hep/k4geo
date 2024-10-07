from argparse import ArgumentParser
import re


def main():
    parser = ArgumentParser(description="Printout overlaps from overlap dump file")
    parser.add_argument("filename")
    args = parser.parse_args()

    pattern = re.compile(r".+\ G4Exception-START\ .+(\n.+)+\ G4Exception-END\ .+")
    with open(args.filename) as file:
        matches = pattern.finditer(file.read())
        for match in matches:
            print(match.group())


if __name__ == "__main__":
    main()
