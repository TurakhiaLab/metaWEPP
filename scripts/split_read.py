#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession is in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is not in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, json, os, sys, subprocess, tempfile, shutil, re
from pathlib import Path

# ────────────────────────── helpers ──────────────────────────
def open_in(fn):
    if fn.endswith(".gz"):
        return subprocess.Popen(["pigz", "-dc", fn], stdout=subprocess.PIPE, text=True).stdout
    return open(fn, "r")

def read_id(h):
    rid = h[1:].split()[0]
    if rid.endswith("/1") or rid.endswith("/2"):
        rid = rid[:-2]
    return rid

def load_taxid_map(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                acc, taxid = line.split(None, 1)
                taxid = taxid.strip()
                acc = acc.strip()
                d.setdefault(taxid, []).append(acc)
    return d


def load_classifications(path, tax2acc):
    m = {}
    opn = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with opn(path, mode) as f:
        for line in f:
            if line.startswith("C\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    _, rid, taxid = parts[:3]
                    accs = tax2acc.get(taxid)
                    if accs:
                        m[rid] = accs[:]  
    return m 

def load_ref_accessions(s):
    if s is None:
        return set()
    p = Path(s)
    if p.exists():
        return {l.strip() for l in p.read_text().splitlines() if l.strip()}
    return {x.strip() for x in s.split(",") if x.strip()}

def out_root_for(acc, base_out, acc2dir, is_ref, acc2taxid, taxid2name):
    if is_ref:
        return Path(base_out) / acc2dir.get(acc, acc)
    taxid = acc2taxid.get(acc)
    if taxid:
        name = taxid2name.get(taxid)
        if name:
            return Path(base_out) / "Other_Pathogens" / sanitize_folder_name(name)
    acc_name = acc.rsplit("|", 1)[-1]
    return Path(base_out) / "Other_Pathogens" / acc_name

def writer_for(acc, mate, ext, out_root, refs, acc2dir, acc2taxid, taxid2name,
               writers, processes, pigz_threads, batch_dir):
    if acc in writers:
        return writers[acc]
    root = out_root_for(acc, out_root, acc2dir, acc in refs, acc2taxid, taxid2name)
    root.mkdir(parents=True, exist_ok=True)
    fn = root / f"{acc.split('.',1)[0]}_{mate}{ext}"
    if batch_dir is not None:
        tmp_fn = Path(batch_dir) / (fn.as_posix().replace("/", "__"))
        w = open(tmp_fn, "w")
        writers[acc] = (w, fn, tmp_fn)
        return writers[acc]
    p = subprocess.Popen(
        ["pigz", f"-p{pigz_threads}", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(fn, "wb"),
        text=True, 
        encoding="utf-8",
        bufsize=1 << 20
    )
    writers[acc] = (p.stdin, fn, None)
    processes[acc] = p
    return writers[acc]

def close_writers(writers, processes):
    for (h, _, tmp) in writers.values():
        if h and not h.closed:
            h.close()
    for acc, p in processes.items():
        rc = p.wait()
        if rc != 0:
            raise RuntimeError(f"pigz failed for {acc}")

def batch_compress(writers, pigz_threads):
    groups = {}
    for (handle, final_fn, tmp_fn) in writers.values():
        if tmp_fn is None:
            continue
        groups.setdefault(final_fn.parent, []).append((final_fn, tmp_fn))
    for parent, lst in groups.items():
        parent.mkdir(parents=True, exist_ok=True)
        for final_fn, tmp_fn in lst:
            subprocess.check_call(["pigz", f"-p{pigz_threads}", "-c", tmp_fn], stdout=open(str(final_fn), "wb"))
    # cleanup tmp files
    for (_, _, tmp_fn) in writers.values():
        if tmp_fn and Path(tmp_fn).exists():
            Path(tmp_fn).unlink()

def sanitize_folder_name(name: str) -> str:
    # collapse whitespace → "_", strip, and remove path separators
    s = re.sub(r"\s+", "_", name.strip())
    return s.replace("/", "_").replace("\\", "_")

def load_taxon_names(kraken_report_path):
    """
    Return {taxid: scientific_name} parsed from a Kraken2 report.
    Works for tab-separated standard Kraken reports; falls back to regex.
    """
    taxid2name = {}
    with open(kraken_report_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                taxid = parts[4].strip()
                name  = parts[5].strip()
                taxid2name[taxid] = name
            else:
                # fallback for space-separated display
                m = re.match(r"\s*\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(.*)$", line)
                if m:
                    taxid = m.group(1)
                    name  = m.group(2).strip()
                    taxid2name[taxid] = name
    return taxid2name

def parse_report_tree(report_path):
    """
    Build taxonomy tree from Kraken2 report using indent depth of the name field.
    Returns: (parent: {taxid->parent_taxid or None}, children: {taxid->[child_taxid...]})
    """
    parent = {}
    children = {}
    stack = []  # [(taxid, indent)]
    with open(report_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            taxid = parts[4].strip()
            name_field = parts[5]  # includes leading spaces for depth
            indent = len(name_field) - len(name_field.lstrip())
            # maintain stack by indent
            while stack and stack[-1][1] >= indent:
                stack.pop()
            if stack:
                p = stack[-1][0]
                parent[taxid] = p
                children.setdefault(p, []).append(taxid)
            else:
                parent[taxid] = None
            stack.append((taxid, indent))
    return parent, children

def descendants_of(taxid, children):
    """All descendants (not including taxid itself)."""
    out = set()
    todo = list(children.get(taxid, []))
    while todo:
        cur = todo.pop()
        if cur in out:
            continue
        out.add(cur)
        todo.extend(children.get(cur, []))
    return out

def build_recipients(acc2dir, acc2taxid, parent, children):
    """
    For each panel accession (key in acc2dir), compute taxids it should receive:
    own taxid + descendants + immediate parent.
    Returns: {taxid -> set(accessions)}
    """
    rec = {}
    for acc in acc2dir.keys():
        t = acc2taxid.get(acc)
        if not t:
            continue
        capture = {t}
        capture |= descendants_of(t, children)
        p = parent.get(t)
        if p:
            capture.add(p)
        for tv in capture:
            rec.setdefault(tv, set()).add(acc)
    return rec

def load_classified_taxids(path):
    """Return {read_id -> classified_taxid} from Kraken output."""
    m = {}
    opn = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with opn(path, mode) as f:
        for line in f:
            if line.startswith("C\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    _, rid, taxid = parts[:3]
                    m[rid] = taxid
    return m


# ───────────────────────── splitting ─────────────────────────
def split_fastq(
    fq, mate, read2taxid, refs, acc2dir, acc2taxid, taxid2name,
    out_root, pigz_threads, batch_dir, recipients_by_taxid, tax2acc
):
    writers = {}
    processes = {}
    ext = ".fq.gz"
    with open_in(fq) as ih:
        it = iter(ih)
        while True:
            h = next(it, None)
            if h is None:
                break
            s = next(it, "")
            p = next(it, "")
            q = next(it, "")
            if not q:
                break
            rid = read_id(h)
            t = read2taxid.get(rid)
            if not t:
                continue

            # Preferred: send to all recipient panel accessions whose capture set contains this taxid
            acc_targets = list(recipients_by_taxid.get(t, []))

            # Fallback: if no recipient claims it, send to any accession(s) that map directly to this taxid
            if not acc_targets:
                acc_targets = tax2acc.get(t, []) or []

            # Nothing to write? skip.
            if not acc_targets:
                continue

            for acc in acc_targets:
                w, _, _ = writer_for(
                    acc, mate, ext, out_root, refs, acc2dir, acc2taxid, taxid2name,
                    writers, processes, pigz_threads, batch_dir
                )
                w.write(h); w.write(s); w.write(p); w.write(q)

    close_writers(writers, processes)
    if batch_dir is not None:
        batch_compress(writers, pigz_threads)

# ─────────────────────────── CLI ────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", "--kraken-out", required=True)
    ap.add_argument("-m", "--mapping", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2")
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--ref-accessions")
    ap.add_argument("--acc2dir", required=True)
    ap.add_argument("--pigz-threads", type=int, default=4)
    ap.add_argument("--batch-compress", action="store_true")
    ap.add_argument("--kraken-report", required=True)   # ← NEW
    return ap.parse_args()

def main():
    a = parse_args()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)

    refs = load_ref_accessions(a.ref_accessions)
    acc2dir = json.load(open(a.acc2dir))

    tax2acc = load_taxid_map(a.mapping)  # {taxid -> [acc, ...]}
    acc2taxid = {acc: taxid for taxid, accs in tax2acc.items() for acc in accs}

    parent, children = parse_report_tree(a.kraken_report)
    recipients_by_taxid = build_recipients(acc2dir=acc2dir, acc2taxid=acc2taxid,
                                           parent=parent, children=children)

    taxid2name = load_taxon_names(a.kraken_report)
    read2taxid = load_classified_taxids(a.kraken_out)
    if not read2taxid:
        sys.exit("No classified reads matched mapping.")

    batch_dir = tempfile.mkdtemp(prefix="split_tmp_") if a.batch_compress else None
    try:
        if a.r2:
            pid = os.fork()
            if pid == 0:
                split_fastq(
                    a.r2, "R2", read2taxid, refs, acc2dir, acc2taxid, taxid2name,
                    a.out_dir, a.pigz_threads, batch_dir, recipients_by_taxid, tax2acc
                )
                os._exit(0)
            split_fastq(
                a.r1, "R1", read2taxid, refs, acc2dir, acc2taxid, taxid2name,
                a.out_dir, a.pigz_threads, batch_dir, recipients_by_taxid, tax2acc
            )
            os.waitpid(pid, 0)
        else:
            split_fastq(
                a.r1, "R1", read2taxid, refs, acc2dir, acc2taxid, taxid2name,
                a.out_dir, a.pigz_threads, batch_dir, recipients_by_taxid, tax2acc
            )
    finally:
        if batch_dir and Path(batch_dir).exists():
            shutil.rmtree(batch_dir, ignore_errors=True)

if __name__ == "__main__":
    main()