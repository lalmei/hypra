"""Write a Shields.io endpoint for the total coverage reported by pytest-cov."""

import argparse
import json
from pathlib import Path
from xml.etree import ElementTree


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coverage_xml", type=Path)
    parser.add_argument("badge_json", type=Path)
    args = parser.parse_args()

    report = ElementTree.parse(args.coverage_xml).getroot()
    covered = int(report.attrib["lines-covered"]) + int(report.attrib["branches-covered"])
    total = int(report.attrib["lines-valid"]) + int(report.attrib["branches-valid"])
    percent = 100 * covered / total
    color = "red"
    if percent >= 90:
        color = "brightgreen"
    elif percent >= 70:
        color = "yellow"

    badge = {
        "schemaVersion": 1,
        "label": "coverage",
        "message": f"{percent:.2f}%",
        "color": color,
    }
    args.badge_json.write_text(json.dumps(badge) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
