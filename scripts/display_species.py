#!/usr/bin/env python3
"""Filter Kraken reports to display species above a percentage threshold."""

import argparse
from pathlib import Path
from typing import List, Tuple


def parse_report(report_path: Path, threshold: float) -> List[Tuple[float, str, str]]:
    entries: List[Tuple[float, str, str]] = []
    with report_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            try:
                percent = float(fields[0])
            except ValueError:
                continue
            rank = fields[3].strip()
            if rank.upper() != "S":
                continue
            if percent < threshold:
                continue
            taxid = fields[4].strip()
            name = fields[5].strip()
            entries.append((percent, taxid, name))
    entries.sort(key=lambda x: x[0], reverse=True)
    return entries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, help="Kraken report file (kreport format)")
    parser.add_argument("--out", required=True, help="Output summary text file")
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.0,
        help="Minimum percentage (0-100) of reads assigned at species level",
    )
    args = parser.parse_args()

    report_path = Path(args.report)
    if not report_path.exists():
        raise SystemExit(f"Report file not found: {report_path}")

    entries = parse_report(report_path, args.threshold)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = f"Species with > {args.threshold}% abundance" if args.threshold else "Species abundance"

    if entries:
        lines = [header]
        for percent, taxid, name in entries:
            lines.append(f"- {percent:.2f}%\t{taxid}\t{name}")
    else:
        lines = [header, "(no species met the threshold)"]

    out_path.write_text("\n".join(lines) + "\n")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
