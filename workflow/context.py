from __future__ import annotations

import json
import os
import re
import sys
import gzip
import bz2
import lzma


from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional


class ContextError(RuntimeError):
    """Raised when the workflow context cannot be constructed."""


@dataclass
class MetaWEPPContext:
    """Aggregate and expose derived configuration for the metaWEPP workflow."""

    config: Dict[str, object]
    config_dir: Path = field(default_factory=lambda: Path("config"))

    def __post_init__(self) -> None:
        self.dir_name: str = self._require_config("DIR")

        self.wepp_root = Path("WEPP")
        self.wepp_workflow = self.wepp_root / "workflow" / "Snakefile"
        self.wepp_results_dir = Path("WEPP") / "results"
        self.wepp_data_dir = Path("WEPP") / "data"
        self.wepp_config = Path("WEPP") / "config" / "config.yaml"
        self.wepp_cmd_log = Path("WEPP") / "cmd_log"
        self.config_path = Path("config") / "config.yaml"

        self.dashboard_enabled: bool = bool(self.config.get("DASHBOARD_ENABLED", False))
        pathogen_str = str(self.config.get("PATHOGENS_FOR_DASHBOARD", ""))
        self.pathogens_for_dashboard: List[str] = [p.strip() for p in pathogen_str.split(",") if p.strip()]

        self.kraken_db: Path = Path(self._require_config("KRAKEN_DB"))
        self.original_map = self.kraken_db / "seqid2taxid.map"
        self.taxid_map_path = self.config_dir / "accession2taxid.map"
        self.config_dir.mkdir(parents=True, exist_ok=True)
        self._convert_taxid_map()

        self.acc2taxid = self._load_acc2taxid()

        self.is_single_end: bool = str(self.config.get("SEQUENCING_TYPE", "d")).lower() in {"s", "n"}
        self.fasta_reads_root = Path("data") / self.dir_name
        if not self.fasta_reads_root.exists():
            raise ContextError(f"Input folder {self.dir_name} does not exist under data/")

        self.fq1, self.fq2 = self._select_fastqs()

        self.pathogen_root = Path("data") / "pathogens_for_wepp"

        self.sample_tag = self.dir_name
        self.out_root = Path("results") / Path(self.fq1).parent.name

        self.acc2dir_json_path = self.config_dir / "acc2dirname.json"
        self.acc2classifieddir_json_path = self.config_dir / "acc2classified_dir.json"
        self.acc2covered_json_path = self.config_dir / "acc2covered.json"
        self.add_ref_sentinel = self.out_root / ".add_ref_mat.done"
        self.split_sentinel = self.out_root / ".split_done"
        self.kraken_out = self.out_root / "kraken_output.txt"
        self.kraken_report = self.out_root / "kraken_report.txt"
        self.visualization = self.out_root / "classification_proportions.png"

        self.acc2fasta: Dict[str, str] = {}
        self.acc2pb: Dict[str, str] = {}
        self.acc2jsonl: Dict[str, str] = {}
        self.acc2dir: Dict[str, str] = {}
        self.headers: List[str] = []
        self.header2taxid: Dict[str, str] = {}
        self.ref_accessions: List[str] = []
        self.exclude_taxids_str = ""
        self._warned_empty = False

        self.refresh_inventory()

        self.min_depth = float(self.config.get("MIN_DEPTH_FOR_WEPP", 0))

        self.fq1_gz = self._gz_path(self.fq1)
        self.fq2_gz = None if self.is_single_end else self._gz_path(self.fq2)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _require_config(self, key: str) -> str:
        if key not in self.config:
            raise ContextError(
                f"Missing required config entry '{key}'. "
                "Supply it via config.yaml or --config."
            )
        value = self.config[key]
        if value is None or value == "":
            raise ContextError(f"Config entry '{key}' is empty")
        return str(value)

    def _convert_taxid_map(self) -> None:
        if not self.original_map.exists():
            raise ContextError(
                f"Kraken2 seqid2taxid.map not found at {self.original_map}. "
                "Run add_ref_mat.py or install the database first."
            )
        with self.original_map.open() as infile, self.taxid_map_path.open("w") as outfile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                full_id, taxid = line.split()[:2]
                accession = full_id.split("|")[-1]
                outfile.write(f"{accession}\t{taxid}\n")

    def _load_acc2taxid(self) -> Dict[str, str]:
        acc2taxid: Dict[str, str] = {}
        with self.original_map.open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                full_id, taxid = parts[:2]
                accession = full_id.split("|")[-1]
                acc2taxid[accession] = taxid
        return acc2taxid

    def _select_fastqs(self) -> tuple[Path, Optional[Path]]:
        patterns = ["*_R1.fastq*", "*_R1.fq*", "*.fastq*"]
        r1_candidates: List[Path] = []
        for pat in patterns:
            r1_candidates.extend(sorted(self.fasta_reads_root.glob(pat)))
        if not r1_candidates:
            raise ContextError(
                f"No R1 FASTQ found in {self.fasta_reads_root} (looked for *_R1.fastq*, *_R1.fq*, *.fastq*)"
            )
        fq1 = r1_candidates[0]
        if self.is_single_end:
            print(
                f"Input FASTQs chosen:\n  FQ1 = {fq1}\n  FQ2 = ",
                file=sys.stderr,
            )
            return fq1, None

        r2_patterns = ["*_R2.fastq*", "*_R2.fq*"]
        r2_candidates: List[Path] = []
        for pat in r2_patterns:
            r2_candidates.extend(sorted(self.fasta_reads_root.glob(pat)))
        if not r2_candidates:
            raise ContextError(
                f"No R2 FASTQ found in {self.fasta_reads_root} for paired-end mode (looked for *_R2.fastq*, *_R2.fq*)"
            )
        fq2 = r2_candidates[0]
        print(
            f"Input FASTQs chosen:\n  FQ1 = {fq1}\n  FQ2 = {fq2}",
            file=sys.stderr,
        )
        return fq1, fq2

    def refresh_inventory(self) -> None:
        """Rescan pathogens_for_wepp for references and auxiliary files."""

        fastas = sorted(self.pathogen_root.glob("*/*.[fF][aAn]*"))
        if not fastas and not getattr(self, "_warned_empty", False):
            print(
                "\nWARNING: No pathogens given for variant analysis with WEPP!\n",
                file=sys.stderr,
            )
            self._warned_empty = True

        acc2fasta: Dict[str, str] = {}
        acc2dir: Dict[str, str] = {}
        acc2pb: Dict[str, str] = {}
        acc2jsonl: Dict[str, str] = {}
        headers: List[str] = []
        header2taxid: Dict[str, str] = {}
        ref_accessions: List[str] = []

        for fasta in fastas:
            header_line = self._read_first_line(fasta)
            header = header_line[1:] if header_line.startswith(">") else header_line
            accession = header.split()[0]
            dir_name = fasta.parent.name

            acc2fasta[accession] = str(fasta.resolve())
            acc2dir[accession] = dir_name
            headers.append(header)
            ref_accessions.append(accession)

            taxid = self.acc2taxid.get(accession)
            if taxid:
                header2taxid[header] = taxid
            else:
                print(
                    f"WARNING: no taxid for accession {accession} (header '{header}') in {self.original_map}",
                    file=sys.stderr,
                )

            pb_path = self._infer_pb_path(fasta.parent)
            if pb_path:
                acc2pb[accession] = pb_path

            jsonl_path = self._infer_jsonl_path(fasta.parent)
            if jsonl_path:
                acc2jsonl[accession] = jsonl_path

        self.acc2fasta = acc2fasta
        self.acc2dir = acc2dir
        self.acc2pb = acc2pb
        self.acc2jsonl = acc2jsonl
        self.headers = headers
        self.header2taxid = header2taxid
        self.ref_accessions = ref_accessions

        taxids = {
            str(int(t))
            for t in self.header2taxid.values()
            if str(t).strip().isdigit()
        }
        self.exclude_taxids_str = re.sub(
            r"[^0-9,]", "", ",".join(sorted(taxids, key=int))
        )

        self._write_json(self.acc2dir_json_path, self.acc2dir)

    def _infer_pb_path(self, directory: Path) -> Optional[str]:
        match = self._first_match(directory, ("*.pb.gz", "*.pb"))
        if match:
            return str(match.resolve())
        fallback = directory / "viz.pb.gz"
        return str(fallback.resolve()) if fallback.exists() else str(fallback)

    def _infer_jsonl_path(self, directory: Path) -> Optional[str]:
        match = self._first_match(directory, ("*.jsonl.gz", "*.jsonl"))
        if match:
            return str(match.resolve())
        return None

    @staticmethod
    def _read_first_line(path: Path) -> str:
        p = str(path)
        if p.endswith(".gz"):
            opener = lambda: gzip.open(p, "rt", encoding="utf-8")
        elif p.endswith(".bz2"):
            opener = lambda: bz2.open(p, "rt", encoding="utf-8")
        elif p.endswith(".xz"):
            opener = lambda: lzma.open(p, "rt", encoding="utf-8")
        else:
            opener = lambda: open(p, "rt", encoding="utf-8")

        try:
            with opener() as handle:
                line = handle.readline().strip()
        except UnicodeDecodeError as e:
            raise ContextError(f"Could not decode {path} (is it compressed?): {e}") from e

        if not line:
            raise ContextError(f"{path} appears to be empty")
        return line

    @staticmethod
    def _first_match(directory: Path, patterns: Iterable[str]) -> Optional[Path]:
        for pattern in patterns:
            for candidate in directory.glob(pattern):
                return candidate
        return None

    @staticmethod
    def _write_json(path: Path, payload: Dict[str, str]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as handle:
            json.dump(payload, handle, indent=2)

    @staticmethod
    def _gz_path(path: Path | None) -> Optional[str]:
        if path is None:
            return None
        text = str(path)
        return text if text.endswith(".gz") else f"{text}.gz"

    # ------------------------------------------------------------------
    # Public helpers consumed by Snakemake rules
    # ------------------------------------------------------------------
    def tag(self, accession: str) -> str:
        return f"{self.dir_for_acc(accession)}_{self.sample_tag}"

    def dir_for_acc(self, accession: str) -> str:
        self.refresh_inventory()
        if accession not in self.acc2dir:
            raise ContextError(
                f"Unknown accession '{accession}'. Rerun add_ref_mat.py to add it."
            )
        return self.acc2dir[accession]

    def wepp_tree(self, accession: str) -> str:
        dir_tag = self.tag(accession)
        tree_name = os.path.basename(self.pb_for(accession))
        return str(self.wepp_data_dir / dir_tag / tree_name)

    def wepp_ref(self, accession: str) -> str:
        dir_tag = self.tag(accession)
        fasta_name = os.path.basename(self.fasta_for(accession))
        return str(self.wepp_data_dir / dir_tag / fasta_name)

    def wepp_jsonl(self, accession: str) -> str:
        src = self.jsonl_for(accession)
        if not src:
            return ""
        dir_tag = self.tag(accession)
        return str(self.wepp_data_dir / dir_tag / os.path.basename(src))

    def fasta_for(self, accession: str) -> str:
        self.refresh_inventory()
        try:
            return self.acc2fasta[accession]
        except KeyError as exc:
            raise ContextError(
                f"No FASTA found for accession '{accession}'."
            ) from exc

    def pb_for(self, accession: str) -> str:
        self.refresh_inventory()
        path = self.acc2pb.get(accession)
        if not path:
            raise ContextError(
                f"No MAT (.pb) file found next to accession '{accession}'. "
                "Run scripts/add_ref_mat.py to install it."
            )
        if not Path(path).exists():
            raise ContextError(
                f"Expected MAT (.pb) file for accession '{accession}' at {path}, "
                "but it does not exist."
            )
        return path

    def jsonl_for(self, accession: str) -> str:
        self.refresh_inventory()
        return self.acc2jsonl.get(accession, "")

    def get_ref_accessions(self) -> List[str]:
        self.refresh_inventory()
        seen = set()
        unique: List[str] = []
        for accession in self.ref_accessions:
            if accession not in seen:
                seen.add(accession)
                unique.append(accession)
        return unique

    def build_ref_arg(self) -> str:
        refs = self.get_ref_accessions()
        if not refs:
            return ""
        return " --ref-accessions " + ",".join(refs)

    def classified_json_empty(self) -> bool:
        try:
            with self.acc2classifieddir_json_path.open() as handle:
                return handle.read().strip() == "{}"
        except FileNotFoundError:
            return True

    # ------------------------------------------------------------------
    # Snakemake callback utilities (checkpoint dependent)
    # ------------------------------------------------------------------
    def split_fastqs_for_coverage(self, checkpoints, wc) -> List[str]:
        ckpt = checkpoints.build_acc2classified_dir.get(**wc)
        classified_json = ckpt.output[0]
        with open(classified_json) as handle:
            acc2dir = json.load(handle)
        r1s: List[str] = []
        for acc, out_dir in acc2dir.items():
            acc_stem = acc.split(".", 1)[0]
            r1s.append(str(self.out_root / out_dir / f"{acc_stem}_R1.fq.gz"))
        return r1s

    def final_targets(self, checkpoints, _wc) -> List[str]:
        cls_ckpt = checkpoints.build_acc2classified_dir.get()
        cov_ckpt = checkpoints.coverage_calculate.get()
        try:
            with open(cls_ckpt.output[0]) as handle:
                acc2classified = json.load(handle)
        except FileNotFoundError:
            acc2classified = {}

        try:
            with open(cov_ckpt.output[0]) as handle:
                acc2covered = json.load(handle)
        except FileNotFoundError:
            acc2covered = {}

        files: List[str] = []
        for accession in self.get_ref_accessions():
            out_dir_cov = acc2covered.get(accession)
            out_dir_cls = acc2classified.get(accession)
            if out_dir_cov and out_dir_cls and out_dir_cov == out_dir_cls:
                dir_tag = f"{out_dir_cov}_{self.sample_tag}"
                files.append(str(self.wepp_results_dir / dir_tag / f"{accession}_run.txt"))

        files.append(str(self.fq1))
        if not self.is_single_end and self.fq2:
            files.append(str(self.fq2))
        return files

    def final_targets_dashboard(self, checkpoints, _wc) -> List[str]:
        checkpoints.build_acc2classified_dir.get()
        wanted: List[str] = []

        if self.acc2classifieddir_json_path.exists():
            with self.acc2classifieddir_json_path.open() as handle:
                acc2dir = json.load(handle)
        else:
            acc2dir = {}

        for accession in self.get_ref_accessions():
            out_dir = acc2dir.get(accession)
            if not out_dir:
                continue
            dir_tag = f"{out_dir}_{self.sample_tag}"
            wanted.append(str(self.wepp_results_dir / dir_tag / f"{accession}_dashboard_run.txt"))

        wanted.append(str(self.fq1))
        if not self.is_single_end and self.fq2:
            wanted.append(str(self.fq2))
        return wanted

    def final_results_files(self, checkpoints, wc) -> List[str]:
        ckpt = checkpoints.build_acc2classified_dir.get(**wc)
        with open(ckpt.output[0]) as handle:
            acc2dir = json.load(handle)

        wanted: List[str] = []
        for accession in self.get_ref_accessions():
            out_dir = acc2dir.get(accession)
            if not out_dir:
                continue
            acc_stem = accession.split(".", 1)[0]
            wanted.append(str(self.out_root / out_dir / f"{acc_stem}_R1.fq.gz"))
            if not self.is_single_end:
                wanted.append(str(self.out_root / out_dir / f"{acc_stem}_R2.fq.gz"))
        return wanted

    def run_txts_for_all(self, checkpoints, _wc) -> List[str]:
        checkpoints.build_acc2classified_dir.get()
        if self.acc2covered_json_path.exists():
            with self.acc2covered_json_path.open() as handle:
                acc2dir = json.load(handle)
        else:
            acc2dir = {}

        files: List[str] = []
        for accession in self.get_ref_accessions():
            directory = acc2dir.get(accession)
            if directory:
                dir_tag = f"{directory}_{self.sample_tag}"
                files.append(str(self.wepp_results_dir / dir_tag / f"{accession}_run.txt"))
        return files

    def dashboard_howto_path(self, _wc) -> str:
        return str(self.wepp_cmd_log / f"{self.sample_tag}_dashboard_howto.txt")

    # ------------------------------------------------------------------
    # Logging helpers
    # ------------------------------------------------------------------
    def log_summary(self) -> None:
        self.refresh_inventory()
        print(f"Dashboard = {self.dashboard_enabled}", file=sys.stderr)
        print("PATHOGENS_FOR_DASHBOARD:", ", ".join(self.pathogens_for_dashboard))
        print(f"IS_SINGLE_END = {self.is_single_end}", file=sys.stderr)
        if self.ref_accessions:
            print("Headers of reference fa file:", ", ".join(self.ref_accessions))
            print("Found accession → folder mappings:")
            for acc, dirname in self.acc2dir.items():
                print(f"  {acc} -> {dirname}")


def build_context(config: Dict[str, object]) -> MetaWEPPContext:
    ctx = MetaWEPPContext(config)
    ctx.log_summary()
    return ctx
