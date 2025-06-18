
from __future__ import annotations

"""Utilities for working with nucleotide sequences.

Public API
----------
* ``encode_sequences`` – k‑mer bag‑of‑words encoding with CountVectorizer
* ``gc_fraction``       – GC fraction (0‑1) for one sequence
* ``gc_content``        – vectorised GC fraction for a DataFrame column
"""

from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer

__all__ = [
    "encode_sequences",
    "gc_fraction",
    "gc_content",
]


def encode_sequences(sequences: List[str] | pd.Series, k: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """Encode DNA/RNA sequences into a *k*-mer count matrix.

    Parameters
    ----------
    sequences : list[str] | pandas.Series
        Iterable with nucleotide strings. Case‑insensitive; non‑ACGT letters are
        kept as‑is (``CountVectorizer`` treats every char).
    k : int, default ``3``
        Length of the *k*-mer. If *k* is larger than *len(seq)* the corresponding
        row in the output matrix will be all zeros (no k‑mers can be formed).

    Returns
    -------
    X : numpy.ndarray, shape ``(n_sequences, n_kmers)``
        Dense count matrix (can be empty ``(0, 0)`` if ``sequences`` is empty
        or *k*‑mer extraction failed).
    feature_names : numpy.ndarray, shape ``(n_kmers,)``
        Array with the *k*-mer strings in the order of columns in *X*.
    """
    # ---- early exit: empty input -----------------------------------------
    if sequences is None or len(sequences) == 0:
        return np.empty((0, 0), dtype=int), np.array([], dtype=str)

    # Coerce to list[str] & uppercase for deterministic vocabulary
    seqs = [str(s).upper() for s in sequences]

    vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k))
    try:
        X = vectorizer.fit_transform(seqs).toarray()
        feature_names = vectorizer.get_feature_names_out()
    except ValueError:
        # Happens when **none** of the sequences can yield a single k‑mer
        X = np.zeros((len(seqs), 0), dtype=int)
        feature_names = np.array([], dtype=str)

    return X, feature_names


def gc_fraction(seq: str) -> float:  # noqa: D401 – simple helper
    """Return GC fraction (*0 ≤ x ≤ 1*) for *seq* (string).

    Empty sequence → ``0.0``.
    """
    if not seq:
        return 0.0
    seq_up = seq.upper()
    gc = sum(base in "GC" for base in seq_up)
    return gc / len(seq_up)


def gc_content(df: pd.DataFrame, col: str = "sequence") -> pd.Series:
    """Vectorised *GC* fraction for a DataFrame column.

    Parameters
    ----------
    df : pandas.DataFrame
    col : str, default ``"sequence"``
        Column with nucleotide strings.

    Returns
    -------
    pandas.Series
        GC fraction for every row (aligned with ``df.index``).
    """
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame")

    return df[col].astype(str).apply(gc_fraction)
