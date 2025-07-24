#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession is in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is not in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, json, os, sys, subprocess, tempfile, shutil
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
                d[taxid.strip()] = acc.strip()
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
                    acc = tax2acc.get(taxid)
                    if acc:
                        m[rid] = acc
    return m

def load_ref_accessions(s):
    if s is None:
        return set()
    p = Path(s)
    if p.exists():
        return {l.strip() for l in p.read_text().splitlines() if l.strip()}
    return {x.strip() for x in s.split(",") if x.strip()}

def out_root_for(acc, base_out, acc2dir, is_ref):
    if is_ref:
        return Path(base_out) / acc2dir.get(acc, acc)
    # take the token after the last '|'
    acc_name = acc.rsplit("|", 1)[-1]
    return Path(base_out) / "Other_Pathogens" / acc_name

def writer_for(acc, mate, ext, out_root, refs, acc2dir, writers, processes, pigz_threads, batch_dir):
    if acc in writers:
        return writers[acc]
    root = out_root_for(acc, out_root, acc2dir, acc in refs)
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

# ───────────────────────── splitting ─────────────────────────
def split_fastq(fq, mate, read2acc, refs, acc2dir, out_root, pigz_threads, batch_dir):
    writers = {}
    processes = {}
    get = read2acc.get
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
            acc = get(read_id(h))
            if not acc:
                continue
            w, _, _ = writer_for(acc, mate, ext, out_root, refs, acc2dir, writers, processes, pigz_threads, batch_dir)
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
    return ap.parse_args()

def main():
    a = parse_args()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    refs = load_ref_accessions(a.ref_accessions)
    tax2acc = load_taxid_map(a.mapping)
    read2acc = load_classifications(a.kraken_out, tax2acc)
    if not read2acc:
        sys.exit("No classified reads matched mapping.")
    acc2dir = json.load(open(a.acc2dir))
    batch_dir = None
    if a.batch_compress:
        batch_dir = tempfile.mkdtemp(prefix="split_tmp_")
    try:
        if a.r2:
            pid = os.fork()
            if pid == 0:
                split_fastq(a.r2, "R2", read2acc, refs, acc2dir, a.out_dir, a.pigz_threads, batch_dir)
                os._exit(0)
            split_fastq(a.r1, "R1", read2acc, refs, acc2dir, a.out_dir, a.pigz_threads, batch_dir)
            os.waitpid(pid, 0)
        else:
            split_fastq(a.r1, "R1", read2acc, refs, acc2dir, a.out_dir, a.pigz_threads, batch_dir)
    finally:
        if batch_dir and Path(batch_dir).exists():
            shutil.rmtree(batch_dir, ignore_errors=True)

if __name__ == "__main__":
    main()

