from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any, Literal

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# 1. базовый анализ
# -----------------------------------------------------------------------------
def analyze_data(df: pd.DataFrame) -> None:
    """Вывод describe() + корреляционная матрица."""
    try:
        print(df.describe())
        sns.heatmap(df.corr(numeric_only=True), annot=True, cmap="coolwarm")
        plt.show()
    except Exception as exc:
        print(f"Error analyzing data: {exc}")


def generate_statistics(df: pd.DataFrame) -> dict[str, Any]:
    """Статистика по числовым / категориальным столбцам."""
    stats: dict[str, Any] = {}
    try:
        stats["numeric"] = df.describe().to_dict()
        cat_cols = df.select_dtypes(include=["object", "category"]).columns
        stats["categorical"] = {c: df[c].value_counts().to_dict() for c in cat_cols}
    except Exception as exc:
        print(f"Error generating statistics: {exc}")
    return stats


def correlation_analysis(df: pd.DataFrame) -> None:
    """Тепловая карта корреляций."""
    try:
        corr = df.corr(numeric_only=True)
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap="coolwarm", fmt=".2f")
        plt.title("Correlation Matrix")
        plt.show()
    except Exception as exc:
        print(f"Error in correlation analysis: {exc}")


# -----------------------------------------------------------------------------
# 2. PCA
# -----------------------------------------------------------------------------
def extract_pca(df: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """Возвращает DataFrame с главными компонентами с улучшенной валидацией."""
    from sklearn.preprocessing import StandardScaler

    num_cols = df.select_dtypes(include=["float64", "int"]).columns
    if len(num_cols) == 0:
        raise ValueError("Нет числовых столбцов для PCA")

    X = df[num_cols].copy()
    X = X.dropna()

    if X.empty:
        raise ValueError("Нет данных после удаления NaN")

    if X.shape[0] < 2:
        raise ValueError("Недостаточно строк для PCA (нужно минимум 2)")

    if X.shape[1] < n_components:
        n_components = X.shape[1]
        logger.warning(f"Снижено количество компонент до {n_components}")

    # Нормализация данных
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    if n_components < 1:
        raise ValueError("Не осталось переменных для PCA")

    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X_scaled)

    result = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_components)])
    result.attrs["explained_variance_ratio"] = pca.explained_variance_ratio_

    return result


# -----------------------------------------------------------------------------
# 3. LabelEncoder wrapper (для Pipeline)
# -----------------------------------------------------------------------------
class LabelEncoderWrapper(LabelEncoder):  # type: ignore[misc]
    def fit(self, X: Sequence[Any], y: Any | None = None) -> LabelEncoderWrapper:
        super().fit(X)
        return self

    def transform(self, X: Sequence[Any]) -> Any:
        return super().transform(X).reshape(-1, 1)

    def fit_transform(self, X: Sequence[Any], y: Any | None = None) -> Any:
        return super().fit_transform(X).reshape(-1, 1)


# -----------------------------------------------------------------------------
# 4. отбор признаков
# -----------------------------------------------------------------------------
def select_features(
    df: pd.DataFrame,
    target_column: str,
    *,
    k: int = 10,
    encoding_method: Literal["onehot", "label"] = "onehot",
) -> list[str]:
    """Возвращает список из *k* лучших признаков."""
    X = df.drop(columns=[target_column])
    y = df[target_column]

    X = X.loc[:, X.var(numeric_only=True) > 0]  # zero-variance filter
    if len(set(y)) < 2:
        raise ValueError("Целевая переменная должна содержать ≥ 2 класса.")

    cat_cols = X.select_dtypes(include=["object", "category"]).columns.tolist()
    num_cols = X.select_dtypes(include=["float64", "int"]).columns.tolist()

    cat_transformer = (
        OneHotEncoder(drop="first", handle_unknown="ignore")
        if encoding_method == "onehot"
        else Pipeline(steps=[("label_encoder", LabelEncoderWrapper())])
    )

    preprocessor = ColumnTransformer(
        [("num", "passthrough", num_cols), ("cat", cat_transformer, cat_cols)]
    )

    pipe = Pipeline(
        [
            ("prep", preprocessor),
            ("var_thr", VarianceThreshold(0.0)),
            ("kbest", SelectKBest(f_classif, k=k)),
        ]
    )

    if y.dtype.kind in {"O", "U"} or y.dtype.name == "category":
        y = LabelEncoder().fit_transform(y)

    pipe.fit(X, y)

    all_feature_names = get_feature_names(
        preprocessor, X, num_cols, cat_cols, encoding_method
    )
    mask = pipe.named_steps["kbest"].get_support()
    return [feat for feat, keep in zip(all_feature_names, mask) if keep]


# -----------------------------------------------------------------------------
# 5. восстановление имён после ColumnTransformer
# -----------------------------------------------------------------------------
def get_feature_names(
    preprocessor: ColumnTransformer,
    X: pd.DataFrame,
    numeric_cols: Sequence[str],
    categorical_cols: Sequence[str],
    encoding_method: Literal["onehot", "label"],
) -> list[str]:
    """Имена фичей на выходе ColumnTransformer."""
    names: list[str] = list(numeric_cols)

    try:
        if encoding_method == "onehot":
            cat_tr = next(
                tr for name, tr, _ in preprocessor.transformers if name == "cat"
            )
            if not hasattr(cat_tr, "categories_"):
                cat_tr.fit(X[categorical_cols])
            ohe_names = cat_tr.get_feature_names_out(categorical_cols)
            names.extend(ohe_names.tolist())
        else:  # label-encoding
            names.extend(categorical_cols)
    except Exception as exc:
        print(f"Error in get_feature_names: {exc}")

    return names
