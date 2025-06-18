# data_ingestion.py
from __future__ import annotations

import json
import logging
import os
from typing import Generator, Optional, Union

import pandas as pd

logger = logging.getLogger("data_ingestion")

# ──────────────────────────────────────────────────────────────
# Soft-imports heavy deps (только warning, проект не падает)
# ──────────────────────────────────────────────────────────────
try:
    from cyvcf2 import VCF           # noqa: N812  (pep-8 alias)
except ImportError:                                      # pragma: no cover
    VCF = None
    logger.warning("'cyvcf2' not installed → VCF loading disabled")

try:
    import pysam
except ImportError:                                      # pragma: no cover
    pysam = None
    logger.warning("'pysam' not installed → BAM loading disabled")

try:
    import gffutils
except ImportError:                                      # pragma: no cover
    gffutils = None
    logger.warning("'gffutils' not installed → GFF loading disabled")

try:
    import h5py
except ImportError:                                      # pragma: no cover
    h5py = None
    logger.warning("'h5py' not installed → HDF5 loading disabled")

try:
    from .DataIngestionFunc.auto_gff_loader import auto_load_gff
except ImportError:                                      # pragma: no cover
    auto_load_gff = None
    logger.warning("'auto_load_gff' unavailable → GFF fallback disabled")

# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────
def _ensure_file(path: str) -> None:
    """Проверяет наличие и непустоту файла; при ошибке — вызывает `SystemExit(1)`."""
    if not os.path.exists(path):
        logger.error("File '%s' not found", path)
        raise SystemExit(1)
    if os.path.getsize(path) == 0:
        logger.error("File '%s' is empty", path)
        raise SystemExit(1)


def _return_if_not_empty(
    df: Optional[pd.DataFrame],
) -> Optional[pd.DataFrame]:
    """Возврат `None`, если DataFrame пуст / None."""
    if df is None or df.empty:
        logger.warning("Loaded data is empty")
        return None
    return df


# ──────────────────────────────────────────────────────────────
# Loaders
# ──────────────────────────────────────────────────────────────
def load_csv(
    path: str,
    *,
    chunksize: int | None = None,
    **read_csv_kwargs,
) -> Union[pd.DataFrame, Generator[pd.DataFrame, None, None], None]:
    """
    Read a **CSV** file.

    Parameters
    ----------
    path : str
    chunksize : int, optional
        Stream rows by chunks (memory-friendly).
    **read_csv_kwargs:
        Passed straight into ``pandas.read_csv``.

    Examples
    --------
    >>> df = load_csv("samples.csv")
    >>> for chunk in load_csv("huge.csv", chunksize=10_000):
    ...     print(len(chunk))
    """
    try:
        _ensure_file(path)
        if chunksize:
            return pd.read_csv(path, chunksize=chunksize, **read_csv_kwargs)
        df = pd.read_csv(path, **read_csv_kwargs)
        return _return_if_not_empty(df)
    except pd.errors.EmptyDataError:
        logger.error("CSV '%s' contains no data (EmptyDataError)", path)
    except pd.errors.ParserError as e:
        logger.error("CSV parsing error '%s': %s", path, e)
    except Exception as e:
        logger.exception("Unexpected CSV error '%s': %s", path, e)
    return None


def load_excel(path: str) -> Optional[pd.DataFrame]:
    """Read **XLS/XLSX** via *pandas*."""
    try:
        _ensure_file(path)
        return _return_if_not_empty(pd.read_excel(path))
    except Exception as e:
        logger.exception("Excel error '%s': %s", path, e)
    return None


def load_json(path: str) -> Optional[pd.DataFrame]:
    """Read **JSON** into flattened DataFrame."""
    try:
        _ensure_file(path)
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return _return_if_not_empty(pd.json_normalize(data))
    except Exception as e:
        logger.exception("JSON error '%s': %s", path, e)
    return None


def load_fasta(path: str) -> Optional[pd.DataFrame]:
    """Read **FASTA** → DataFrame(sequence=str)."""
    try:
        _ensure_file(path)
        seqs: list[str] = []
        with open(path) as fh:
            buff = ""
            for line in fh:
                if line.startswith(">"):
                    if buff:
                        seqs.append(buff)
                        buff = ""
                else:
                    buff += line.strip()
            if buff:
                seqs.append(buff)
        return _return_if_not_empty(pd.DataFrame({"sequence": seqs}))
    except Exception as e:
        logger.exception("FASTA error '%s': %s", path, e)
    return None


def load_vcf(path: str) -> Optional[pd.DataFrame]:
    """Read **VCF** (requires *cyvcf2*)."""
    if VCF is None:  # pragma: no cover
        logger.error("Requested VCF but 'cyvcf2' is missing")
        return None
    try:
        _ensure_file(path)
        rows = [
            {
                "CHROM": r.CHROM,
                "POS": r.POS,
                "ID": r.ID,
                "REF": r.REF,
                "ALT": [str(a) for a in r.ALT],
                "QUAL": r.QUAL,
                "FILTER": r.FILTER,
            }
            for r in VCF(path)
        ]
        return _return_if_not_empty(pd.DataFrame(rows))
    except Exception as e:
        logger.exception("VCF error '%s': %s", path, e)
    return None


def load_bam(path: str) -> Optional[pd.DataFrame]:
    """Read **BAM** (requires *pysam*)."""
    if pysam is None:  # pragma: no cover
        logger.error("Requested BAM but 'pysam' is missing")
        return None
    try:
        _ensure_file(path)
        bam = pysam.AlignmentFile(path, "rb")
        rows = [
            {
                "QNAME": r.query_name,
                "FLAG": r.flag,
                "RNAME": bam.get_reference_name(r.reference_id),
                "POS": r.reference_start,
                "MAPQ": r.mapping_quality,
                "CIGAR": r.cigarstring,
                "SEQ": r.query_sequence,
                "QUAL": r.query_qualities,
            }
            for r in bam.fetch()
        ]
        bam.close()
        return _return_if_not_empty(pd.DataFrame(rows))
    except Exception as e:
        logger.exception("BAM error '%s': %s", path, e)
    return None


def load_hdf5(path: str) -> Optional[pd.DataFrame]:
    """Read **HDF5** (all datasets → columns)."""
    if h5py is None:  # pragma: no cover
        logger.error("Requested HDF5 but 'h5py' is missing")
        return None
    try:
        _ensure_file(path)
        with h5py.File(path) as hdf:
            if not hdf.keys():
                logger.warning("HDF5 '%s' has no datasets", path)
                return None
            data = {k: hdf[k][()] for k in hdf.keys()}
        return _return_if_not_empty(pd.DataFrame(data))
    except Exception as e:
        logger.exception("HDF5 error '%s': %s", path, e)
    return None


def load_gff(
    path: str,
) -> Union[pd.DataFrame, Generator[pd.DataFrame, None, None], None]:
    """Read **GFF/GTF** (streaming large files)."""
    if auto_load_gff is None or gffutils is None:  # pragma: no cover
        logger.error("Requested GFF but deps missing")
        return None
    try:
        _ensure_file(path)
        return auto_load_gff(path)
    except Exception as e:
        logger.exception("GFF error '%s': %s", path, e)
    return None


# ──────────────────────────────────────────────────────────────
# Public façade
# ──────────────────────────────────────────────────────────────
def load_data(
    path: str,
    ftype: str,
    **kwargs,
) -> Union[pd.DataFrame, Generator[pd.DataFrame, None, None], None]:
    """
    Universal data loader.

    Parameters
    ----------
    path : str
        File path.
    ftype : str
        One of {'csv','excel','xls','xlsx','json','fasta',
        'vcf','bam','gff','hdf5'}.
    **kwargs :
        Passed verbatim to the underlying `load_*` (e.g. ``chunksize=``).

    Examples
    --------
    >>> df = load_data("data.csv", "csv")
    >>> for chunk in load_data("big.csv", "csv", chunksize=1_000):
    ...     print(len(chunk))
    """
    ftype = ftype.lower()
    loaders = {
        "csv": load_csv,
        "excel": load_excel,
        "xls": load_excel,
        "xlsx": load_excel,
        "json": load_json,
        "fasta": load_fasta,
        "vcf": load_vcf,
        "bam": load_bam,
        "gff": load_gff,
        "hdf5": load_hdf5,
    }
    if ftype not in loaders:
        logger.error("Unsupported file type '%s'", ftype)
        return None
    return loaders[ftype](path, **kwargs)
