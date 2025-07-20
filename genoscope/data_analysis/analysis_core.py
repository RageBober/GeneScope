from __future__ import annotations

# ─── std / typing ────────────────────────────────────────────────────────────
from typing import Any, Dict, List, Literal, Sequence

# ─── third-party ─────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

# -----------------------------------------------------------------------------
# 1. базовый анализ
# -----------------------------------------------------------------------------
def analyze_data(df: pd.DataFrame) -> None:
    """Вывод describe() + корреляционная матрица."""
    try:
        print(df.describe())
        sns.heatmap(df.corr(numeric_only=True), annot=True, cmap="coolwarm")
        plt.show()
    except Exception as exc:  # noqa: BLE001
        print(f"Error analyzing data: {exc}")


def generate_statistics(df: pd.DataFrame) -> Dict[str, Any]:
    """Статистика по числовым / категориальным столбцам."""
    stats: Dict[str, Any] = {}
    try:
        stats["numeric"] = df.describe().to_dict()
        cat_cols = df.select_dtypes(include=["object", "category"]).columns
        stats["categorical"] = {c: df[c].value_counts().to_dict() for c in cat_cols}
    except Exception as exc:  # noqa: BLE001
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
    except Exception as exc:  # noqa: BLE001
        print(f"Error in correlation analysis: {exc}")

# -----------------------------------------------------------------------------
# 2. PCA
# -----------------------------------------------------------------------------
def extract_pca(df: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """Возвращает DataFrame с главными компонентами."""
    num_cols = df.select_dtypes(include=["float64", "int"]).columns
    X = df[num_cols]
    if X.shape[0] < n_components or X.shape[1] < n_components:
        raise ValueError(
            f"Недостаточно данных для извлечения {n_components} главных компонентов."
        )
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X)
    return pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_components)])

# -----------------------------------------------------------------------------
# 3. LabelEncoder wrapper (для Pipeline)
# -----------------------------------------------------------------------------
class LabelEncoderWrapper(LabelEncoder):  # type: ignore[misc]
    def fit(self, X: Sequence[Any], y: Any | None = None) -> "LabelEncoderWrapper":  # noqa: N802
        super().fit(X)
        return self

    def transform(self, X: Sequence[Any]) -> Any:  # noqa: N802
        return super().transform(X).reshape(-1, 1)

    def fit_transform(self, X: Sequence[Any], y: Any | None = None) -> Any:  # noqa: N802
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
) -> List[str]:
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

    all_feature_names = get_feature_names(preprocessor, X, num_cols, cat_cols, encoding_method)
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
) -> List[str]:
    """Имена фичей на выходе ColumnTransformer."""
    names: List[str] = list(numeric_cols)

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
    except Exception as exc:  # noqa: BLE001
        print(f"Error in get_feature_names: {exc}")

    return names
