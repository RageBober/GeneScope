# data_ingestion.py
from __future__ import annotations

import json
import logging
import os
from collections.abc import Callable
from collections.abc import Iterator
from typing import Any
from typing import cast

import pandas as pd

logger = logging.getLogger("data_ingestion")

# ═══════════════════════════════════════════════════════════════
#  Soft-imports «heavy» зависимостей (только warning, проект не падает)
# ═══════════════════════════════════════════════════════════════

# --- cyvcf2 ------------------------------------------------------
try:
    from cyvcf2 import VCF

    VCF: Optional[Any] = VCF  # либо types.ModuleType  | None
except ImportError:
    VCF = None  #   └── объединённый тип = OK
    logger.warning("'cyvcf2' not installed → VCF loading disabled")

# --- pysam -------------------------------------------------------
try:
    import pysam as _pysam

    pysam: Any = _pysam
except ImportError:
    pysam = None  # тип уже Any → нет ошибки
    logger.warning("'pysam' not installed → BAM loading disabled")

# --- gffutils ----------------------------------------------------
try:
    import gffutils as _gffutils

    gffutils: Any = _gffutils
except ImportError:
    gffutils = None
    logger.warning("'gffutils' not installed → GFF loading disabled")

# --- h5py --------------------------------------------------------
try:
    import h5py as _h5py

    h5py: Any = _h5py
except ImportError:
    h5py = None
    logger.warning("'h5py' not installed → HDF5 loading disabled")

# --- auto_gff_loader --------------------------------------------
try:
    from .DataIngestionFunc.auto_gff_loader import auto_load_gff as _auto_load_gff

    auto_load_gff: (
        Callable[..., pd.DataFrame | Iterator[pd.DataFrame] | None] | None
    ) = _auto_load_gff
except ImportError:
    auto_load_gff = None
    logger.warning("'auto_load_gff' unavailable → GFF fallback disabled")

# ═══════════════════════════════════════════════════════════════
#  Тип-алиасы
# ═══════════════════════════════════════════════════════════════
Loader = Callable[
    ...,
    pd.DataFrame | Iterator[pd.DataFrame] | None,
]


# ═══════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════
def _ensure_file(path: str) -> None:
    """
    Гарантирует, что путь указывает **именно на существующий файл**, а не каталог.
    *Пустые* файлы допускаем – это валидный кейс для `pd.read_csv`, который
    сам выбросит `EmptyDataError`, и мы аккуратно вернём `None`.
    """
    if not os.path.exists(path) or not os.path.isfile(path):
        logger.error("Path '%s' is not a regular file", path)
        raise ValueError(f"Path '{path}' is not a regular file")


def _return_if_not_empty(df: pd.DataFrame | None) -> pd.DataFrame | None:
    """Возвращает DataFrame, либо `None`, если он пуст / None."""
    if df is None or df.empty:
        logger.warning("Loaded data is empty")
        return None
    return df


# ═══════════════════════════════════════════════════════════════
#  CSV / Excel / JSON / FASTA
# ═══════════════════════════════════════════════════════════════
def load_csv(
    path: str,
    *,
    chunksize: Optional[int] = None,
    **read_csv_kwargs: Any,
) -> pd.DataFrame | Iterator[pd.DataFrame] | None:
    try:
        _ensure_file(path)
        if chunksize:
            reader = pd.read_csv(path, chunksize=chunksize, **read_csv_kwargs)
            return cast("Iterator[pd.DataFrame]", reader)
        df = pd.read_csv(path, **read_csv_kwargs)
        return _return_if_not_empty(df)
    except pd.errors.EmptyDataError:
        logger.error("CSV '%s' contains no data (EmptyDataError)", path)
    except pd.errors.ParserError as exc:
        logger.error("CSV parsing error '%s': %s", path, exc)
    except Exception as exc:
        logger.exception("Unexpected CSV error '%s': %s", path, exc)
    return None


def load_excel(path: str) -> pd.DataFrame | None:
    try:
        _ensure_file(path)
        return _return_if_not_empty(pd.read_excel(path))
    except Exception as exc:
        logger.exception("Excel error '%s': %s", path, exc)
    return None


def load_json(path: str) -> pd.DataFrame | None:
    try:
        _ensure_file(path)
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return _return_if_not_empty(pd.json_normalize(data))
    except Exception as exc:
        logger.exception("JSON error '%s': %s", path, exc)
    return None


def load_fasta(path: str) -> pd.DataFrame | None:
    try:
        _ensure_file(path)
        seqs: List[str] = []
        with open(path) as fh:
            buf = ""
            for line in fh:
                if line.startswith(">"):
                    if buf:
                        seqs.append(buf)
                        buf = ""
                else:
                    buf += line.strip()
            if buf:
                seqs.append(buf)
        return _return_if_not_empty(pd.DataFrame({"sequence": seqs}))
    except Exception as exc:
        logger.exception("FASTA error '%s': %s", path, exc)
    return None


# ═══════════════════════════════════════════════════════════════
#  Bio-форматы (VCF, BAM, HDF5, GFF)
# ═══════════════════════════════════════════════════════════════
def load_vcf(path: str) -> pd.DataFrame | None:
    if VCF is None:  # pragma: no cover
        logger.error("Requested VCF but 'cyvcf2' is missing")
        return None
    try:
        _ensure_file(path)
        rows = [
            {
                "CHROM": rec.CHROM,
                "POS": rec.POS,
                "ID": rec.ID,
                "REF": rec.REF,
                "ALT": [str(a) for a in rec.ALT],
                "QUAL": rec.QUAL,
                "FILTER": rec.FILTER,
            }
            for rec in VCF(path)
        ]
        return _return_if_not_empty(pd.DataFrame(rows))
    except Exception as exc:
        logger.exception("VCF error '%s': %s", path, exc)
    return None


def load_bam(path: str) -> pd.DataFrame | None:
    if pysam is None:  # pragma: no cover
        logger.error("Requested BAM but 'pysam' is missing")
        return None
    try:
        _ensure_file(path)
        
        # Use context manager for proper resource cleanup
        with pysam.AlignmentFile(path, "rb") as bam:
            rows = [
                {
                    "QNAME": rec.query_name,
                    "FLAG": rec.flag,
                    "RNAME": bam.get_reference_name(rec.reference_id),
                    "POS": rec.reference_start,
                    "MAPQ": rec.mapping_quality,
                    "CIGAR": rec.cigarstring,
                    "SEQ": rec.query_sequence,
                    "QUAL": rec.query_qualities,
                }
                for rec in bam.fetch()
            ]
        
        return _return_if_not_empty(pd.DataFrame(rows))
    except pysam.utils.SamtoolsError as exc:
        logger.error("SAM/BAM format error '%s': %s", path, exc)
    except MemoryError:
        logger.error("Not enough memory to load BAM file '%s'", path)
    except Exception as exc:
        logger.exception("BAM error '%s': %s", path, exc)
    return None


def load_hdf5(path: str) -> pd.DataFrame | None:
    if h5py is None:  # pragma: no cover
        logger.error("Requested HDF5 but 'h5py' is missing")
        return None
    try:
        _ensure_file(path)
        
        # Use context manager for proper resource cleanup
        with h5py.File(path, 'r') as hdf:
            if not hdf.keys():
                logger.warning("HDF5 '%s' has no datasets", path)
                return None
            
            # Safely read data with memory consideration
            data = {}
            for key in hdf.keys():
                dataset = hdf[key]
                # Check dataset size to prevent memory issues
                if dataset.size > 1000000:  # Limit to 1M elements
                    logger.warning("HDF5 dataset '%s' too large (%d elements), skipping", key, dataset.size)
                    continue
                data[key] = dataset[()]
        
        return _return_if_not_empty(pd.DataFrame(data))
    except OSError as exc:
        logger.error("HDF5 file access error '%s': %s", path, exc)
    except MemoryError:
        logger.error("Not enough memory to load HDF5 file '%s'", path)
    except Exception as exc:
        logger.exception("HDF5 error '%s': %s", path, exc)
    return None


def load_gff(
    path: str,
) -> pd.DataFrame | Iterator[pd.DataFrame] | None:
    if auto_load_gff is None or gffutils is None:  # pragma: no cover
        logger.error("Requested GFF but deps missing")
        return None
    try:
        _ensure_file(path)
        # auto_load_gff имеет тип Any → приводим к Loader
        return auto_load_gff(path)
    except Exception as exc:
        logger.exception("GFF error '%s': %s", path, exc)
    return None


# ═══════════════════════════════════════════════════════════════
#  Facade
# ═══════════════════════════════════════════════════════════════
def load_data(
    path: str,
    ftype: str,
    **kwargs: Any,
) -> pd.DataFrame | Iterator[pd.DataFrame] | None:
    """
    Universal data loader.

    Examples
    --------
    >>> df = load_data("data.csv", "csv")
    >>> for chunk in load_data("big.csv", "csv", chunksize=1_000):
    ...     print(len(chunk))
    """
    ftype = ftype.lower()
    loaders: Dict[str, Loader] = {
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
