from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from sklearn.feature_extraction.text import CountVectorizer

"""Utilities for working with nucleotide sequences.

Public API
----------
* ``encode_sequences`` – k-mer bag-of-words encoding with CountVectorizer
* ``gc_fraction``       – GC fraction (0-1) for one sequence
* ``gc_content``        – vectorised GC fraction for a DataFrame column
"""

__all__ = [
    "encode_sequences",
    "gc_content",
    "gc_fraction",
]


def encode_sequences(
    sequences: Sequence[str] | pd.Series[str],
    k: int = 3,
) -> tuple[NDArray[np.int_], NDArray[np.str_]]:
    """Encode DNA/RNA sequences into a *k*-mer count matrix.

    Parameters
    ----------
    sequences : Sequence[str] | pandas.Series[str]
        Iterable with nucleotide strings.
    k : int, default ``3``
        Length of the *k*-mer.

    Returns
    -------
    X : NDArray[np.int_], shape ``(n_sequences, n_kmers)``
    feature_names : NDArray[np.str_], shape ``(n_kmers,)``
    """
    # ---- early exit ------------------------------------------------------
    if not sequences:  # empty list / Series
        return np.empty((0, 0), dtype=np.int_), np.empty(0, dtype=np.str_)

    seqs: list[str] = [str(s).upper() for s in sequences]

    vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k))
    try:
        X: NDArray[np.int_] = vectorizer.fit_transform(seqs).toarray()
        feature_names: NDArray[np.str_] = vectorizer.get_feature_names_out()
    except ValueError:
        # Ни один k-мер не сформировался (например, k > len(seq) для всех)
        X = np.zeros((len(seqs), 0), dtype=np.int_)
        feature_names = np.empty(0, dtype=np.str_)

    return X, feature_names


def gc_fraction(seq: str) -> float:
    """Return GC fraction (*0 ≤ x ≤ 1*) for ``seq``."""
    if not seq:
        return 0.0
    seq_up = seq.upper()
    gc = sum(base in "GC" for base in seq_up)
    return gc / len(seq_up)


def gc_content(
    df: pd.DataFrame,
    col: str = "sequence",
) -> pd.Series[float]:
    """Vectorised GC fraction for column ``df[col]``."""
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame")

    return df[col].astype(str).apply(gc_fraction).astype(float)
