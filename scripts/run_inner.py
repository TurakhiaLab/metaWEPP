#!/usr/bin/env python3
import argparse, subprocess, os, re, shlex, sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--dir", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--primer_bed", required=True)
parser.add_argument("--min_af", required=True)
parser.add_argument("--min_q", required=True)
parser.add_argument("--max_reads", required=True)
parser.add_argument("--tree", required=True)
parser.add_argument("--ref", required=True)
parser.add_argument("--clade_idx", nargs="?", required=False, default=None)
parser.add_argument("--snakefile", required=True)
parser.add_argument("--workdir", required=True)
parser.add_argument("--cfgfile", required=True)
parser.add_argument("--customconfig", required=False)
parser.add_argument("--cores", default="32")
parser.add_argument("--sequencing_type", required=True)
parser.add_argument("--clade_list", nargs="?", required=False, default=None)
parser.add_argument("--min_prop", required=False, default=None)
parser.add_argument("--min_len", required=False, default=None)
parser.add_argument("--cmd_log", required=True)
parser.add_argument("--pathogens_name", required=True)
parser.add_argument("--dashboard_enabled", action="store_true")
parser.add_argument("--taxonium_file", required=False, default="")
args = parser.parse_args()

def upsert_cmd_log(path, name, cmd_str):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    if path.exists():
        with open(path) as f:
            for line in f:
                k = line.rstrip("\n").split(" : ", 1)
                if len(k) == 2 and k[0] == name:
                    continue
                lines.append(line.rstrip("\n"))
    lines.append(f"{name} : {cmd_str}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.replace(tmp, path)

def parse_config_file(path):
    cfg = {}
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r'^([^:#\s]+)\s*:\s*["\']?(.*?)["\']?\s*$', line)
            if not m:
                raise ValueError(f"Invalid line in config file: {line}")
            k, v = m.groups()
            cfg[k] = v
    return cfg

def csv_tokens(s: str, keep_empty: bool = False):
    if s is None:
        return []
    parts = [p.strip() for p in str(s).split(",")]
    if keep_empty:
        return parts
    return [p for p in parts if p != ""]

def to_csv_str(tokens):
    return ",".join(tokens)

def get_cfg(main_cfg, overlay, key, fallback=None):
    if overlay and key in overlay:
        return overlay[key]
    if fallback is not None:
        return fallback
    return main_cfg.get(key, fallback)

def replicate_or_align_csv(raw, n, field_name):
    vals = csv_tokens(raw or "", keep_empty=False)
    if len(vals) == 1:
        return vals * n
    if len(vals) == n:
        return vals
    raise SystemExit(
        f"[config] {field_name}: expected 1 value or {n} values; got {len(vals)}. Raw='{raw}'"
    )

def resolve_clade_lists_for_pathogens(pathogens, clade_idx_list, clade_list_stream_raw):
    stream = csv_tokens(clade_list_stream_raw or "", keep_empty=True)
    pos = 0
    out = {}
    for i, name in enumerate(pathogens):
        try:
            idx = int(clade_idx_list[i])
        except ValueError:
            raise SystemExit(f"[config] CLADE_IDX for '{name}' is not an integer: '{clade_idx_list[i]}'")

        if idx < -1:
            raise SystemExit(f"[config] CLADE_IDX supports only -1 or >=0; got {idx} for '{name}'")

        if idx == -1:
            if pos < len(stream) and stream[pos] == "":
                pos += 1
            out[name] = ""
            continue

        need = idx + 1
        if pos + need > len(stream):
            have = len(stream) - pos
            raise SystemExit(
                f"[config] CLADE_LIST ran out of items for pathogen '{name}' (need {need}, have {have}). "
                f"CLADE_LIST='{clade_list_stream_raw}' CLADE_IDX='{','.join(map(str,clade_idx_list))}'"
            )
        group = stream[pos:pos + need]
        pos += need
        group = [g for g in group if g != ""]
        out[name] = to_csv_str(group) if group else ""

    leftover = stream[pos:]
    if any(tok != "" for tok in leftover):
        print(f"[warn] CLADE_LIST has {len(leftover)} leftover token(s): {leftover}", file=sys.stderr)

    return out

def choose_for_current(name, mapping, pathogens):
    if name in mapping:
        return mapping[name]
    if "default" in mapping:
        return mapping["default"]
    return mapping[pathogens[0]]

def index_for_current(name, pathogens):
    if name in pathogens:
        return pathogens.index(name)
    if "default" in pathogens:
        return pathogens.index("default")
    raise SystemExit(
        f"[config] PATHOGENS must list '{name}' or include a 'default' entry."
    )

def _normalize_clade_token(token):
    if not token:
        return ""
    parts = [p for p in re.split(r":", token) if p]
    return ",".join(parts)

def per_pathogen_clade_mapping(pathogens, raw):
    tokens_with_empty = [t.strip() for t in str(raw or "").split(",")]

    if len(pathogens) == 1:
        return {pathogens[0]: _normalize_clade_token(tokens_with_empty[0] if tokens_with_empty else "")}

    if len(tokens_with_empty) == len(pathogens):
        return {
            name: _normalize_clade_token(tok)
            for name, tok in zip(pathogens, tokens_with_empty)
        }

    nonempty_tokens = [tok for tok in tokens_with_empty if tok]
    if not nonempty_tokens:
        return {name: "" for name in pathogens}

    if len(nonempty_tokens) == 1 and len(tokens_with_empty) >= 1:
        normalized = _normalize_clade_token(nonempty_tokens[0])
        return {name: normalized for name in pathogens}

    return None

main_cfg = parse_config_file(args.cfgfile)
overlay = parse_config_file(args.customconfig) if args.customconfig else {}

PRIMER_BED = get_cfg(main_cfg, overlay, "PRIMER_BED", args.primer_bed)
MIN_AF = get_cfg(main_cfg, overlay, "MIN_AF", args.min_af)
MIN_Q = get_cfg(main_cfg, overlay, "MIN_Q", args.min_q)
MAX_READS = get_cfg(main_cfg, overlay, "MAX_READS", args.max_reads)
SEQUENCING_TYPE = get_cfg(main_cfg, overlay, "SEQUENCING_TYPE", args.sequencing_type)

MIN_PROP = get_cfg(main_cfg, overlay, "MIN_PROP", args.min_prop)
MIN_LEN = get_cfg(main_cfg, overlay, "MIN_LEN", args.min_len)

# Require explicit PATHOGENS roster so per-pathogen settings resolve deterministically.
PATHOGENS_raw = get_cfg(main_cfg, overlay, "PATHOGENS", None)
if PATHOGENS_raw is None:
    raise SystemExit("[config] PATHOGENS must list one or more pathogen names (comma-separated).")
PATHOGENS = [p for p in csv_tokens(PATHOGENS_raw, keep_empty=False) if p]
if not PATHOGENS:
    raise SystemExit("[config] PATHOGENS must list one or more pathogen names (comma-separated).")

current_name = args.pathogens_name
if current_name not in PATHOGENS and "default" not in PATHOGENS:
    raise SystemExit(
        f"[config] PATHOGENS does not include '{current_name}' and no 'default' fallback is defined."
    )

pathogen_index = index_for_current(current_name, PATHOGENS)

CLADE_IDX_raw_cfg = get_cfg(main_cfg, overlay, "CLADE_IDX", args.clade_idx)
if CLADE_IDX_raw_cfg is None:
    raise SystemExit("[config] CLADE_IDX is required (CSV). Example: '1,0,0'")
CLADE_IDX_list = replicate_or_align_csv(CLADE_IDX_raw_cfg, len(PATHOGENS), "CLADE_IDX")
CLADE_IDX_CUR = CLADE_IDX_list[pathogen_index]

if args.clade_list:
    CLADE_LIST_raw_cfg = args.clade_list
else:
    CLADE_LIST_raw_cfg = get_cfg(main_cfg, overlay, "CLADE_LIST", "")

per_pathogen_clade_list = per_pathogen_clade_mapping(PATHOGENS, CLADE_LIST_raw_cfg)
if per_pathogen_clade_list is None:
    per_pathogen_clade_list = resolve_clade_lists_for_pathogens(
        PATHOGENS,
        CLADE_IDX_list,
        CLADE_LIST_raw_cfg,
    )
CLADE_LIST_CUR = choose_for_current(current_name, per_pathogen_clade_list, PATHOGENS)

MIN_DEPTH_raw = get_cfg(main_cfg, overlay, "MIN_DEPTH_FOR_WEPP", "")
if MIN_DEPTH_raw == "":
    fallback_depth = main_cfg.get("MIN_DEPTH", "10")
    MIN_DEPTH_raw = str(fallback_depth)
MIN_DEPTH_list = replicate_or_align_csv(MIN_DEPTH_raw, len(PATHOGENS), "MIN_DEPTH_FOR_WEPP")
MIN_DEPTH_CUR = MIN_DEPTH_list[pathogen_index]

cmd = [
    "snakemake",
    "--snakefile",
    args.snakefile,
    "--directory",
    args.workdir,
    "--cores",
    args.cores,
    "--use-conda",
    "--config",
    f"DIR={args.dir}",
    f"FILE_PREFIX={args.prefix}",
    f"TREE={args.tree}",
    f"REF={args.ref}",
    f"PRIMER_BED={PRIMER_BED}",
    f"MIN_AF={MIN_AF}",
    f"MIN_Q={MIN_Q}",
    f"MAX_READS={MAX_READS}",
    f"SEQUENCING_TYPE={SEQUENCING_TYPE}",
    f"PATHOGENS={current_name}",
    f"CLADE_LIST={CLADE_LIST_CUR}",
    f"CLADE_IDX={CLADE_IDX_CUR}",
    f"MIN_DEPTH={MIN_DEPTH_CUR}",
]

if MIN_PROP is not None:
    cmd.append(f"MIN_PROP={MIN_PROP}")
if MIN_LEN is not None:
    cmd.append(f"MIN_LEN={MIN_LEN}")

if args.taxonium_file:
    cmd.append(f"TAXONIUM_FILE={args.taxonium_file}")

cmd_log_path = Path(args.cmd_log)

if not args.dashboard_enabled:
    upsert_cmd_log(cmd_log_path, current_name, " ".join(cmd))
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)
else:
    if not cmd_log_path.exists():
        raise SystemExit(f"cmd_log not found: {cmd_log_path}")
    matched = None
    with open(cmd_log_path) as f:
        for line in reversed(f.readlines()):
            k = line.rstrip("\n").split(" : ", 1)
            if len(k) == 2 and k[0] == current_name:
                matched = k[1]
                break
    if matched is None:
        raise SystemExit(f"No command found in cmd_log for PATHOGENS_FOR_DASHBOARD '{current_name}'")
    tokens = shlex.split(matched)
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
        if args.taxonium_file:
            found = False
            for m in range(j, k):
                if tokens[m].startswith("TAXONIUM_FILE="):
                    tokens[m] = f"TAXONIUM_FILE={args.taxonium_file}"
                    found = True
                    break
            if not found:
                tokens.insert(k, f"TAXONIUM_FILE={args.taxonium_file}")
    except ValueError:
        tokens.extend(["--config", "DASHBOARD_ENABLED=True"])
    tokens.extend(["--forcerun", "dashboard_serve"])
    print("Running command (dashboard enabled):", " ".join(tokens))
    subprocess.run(tokens, check=True)
