#!/usr/bin/env python3
import argparse
import shlex
from collections import OrderedDict
from pathlib import Path

def transform(cmd_str: str) -> str:
    tokens = shlex.split(cmd_str)

    # Ensure DASHBOARD_ENABLED=True is in the --config block
    try:
        idx = tokens.index("--config")
        j = idx + 1
        k = j
        while k < len(tokens) and not tokens[k].startswith("-"):
            k += 1
        replaced = False
        for m in range(j, k):
            if tokens[m].startswith("DASHBOARD_ENABLED="):
                tokens[m] = "DASHBOARD_ENABLED=True"
                replaced = True
                break
        if not replaced:
            tokens.insert(k, "DASHBOARD_ENABLED=True")
    except ValueError:
        tokens.extend(["--config", "DASHBOARD_ENABLED=True"])

    # Ensure --forcerun dashboard_serve is present (merge if already there)
    if "--forcerun" in tokens:
        i = tokens.index("--forcerun")
        if i + 1 < len(tokens) and not tokens[i+1].startswith("-"):
            val = tokens[i+1]
            if "dashboard_serve" not in val.split(","):
                tokens[i+1] = val + ",dashboard_serve"
        else:
            tokens.insert(i + 1, "dashboard_serve")
    else:
        tokens.extend(["--forcerun", "dashboard_serve"])

    return " ".join(tokens)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmd-log", required=True, help="Path to WEPP cmd log")
    ap.add_argument("--out", required=True, help="Path to write dashboard how-to text")
    args = ap.parse_args()

    cmd_log = Path(args.cmd_log)
    howto_out = Path(args.out)

    latest = OrderedDict()
    if cmd_log.exists():
        with cmd_log.open() as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                parts = line.split(" : ", 1)
                if len(parts) != 2:
                    continue
                name, cmd_str = parts
                if name in latest:
                    del latest[name]  # keep only the last for each name
                latest[name] = cmd_str

    header = (
        "If you want to run WEPP with dashboard, go to the WEPP main dir:\n"
        "  cd WEPP\n"
        "and run the following commands (one per pathogen):\n"
    )

    lines_out = [header]
    for name, cmd_str in latest.items():
        lines_out.append(f"{name} : {transform(cmd_str)}")

    text = "\n".join(lines_out) + "\n"
    print(text)

    howto_out.parent.mkdir(parents=True, exist_ok=True)
    with howto_out.open("w") as out:
        out.write(text)

if __name__ == "__main__":
    main()
