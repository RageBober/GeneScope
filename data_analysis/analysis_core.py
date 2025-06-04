# analysis_core.py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif, VarianceThreshold
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline


def analyze_data(df: pd.DataFrame) -> None:
    """
    Выполняет базовый анализ данных:
    - Выводит описательную статистику по DataFrame,
    - Строит корреляционную матрицу.

    Параметры:
        df (pd.DataFrame): Данные для анализа.
    """
    try:
        print(df.describe())
        sns.heatmap(df.corr(), annot=True, cmap="coolwarm")
        plt.show()
    except Exception as e:
        print(f"Error analyzing data: {e}")


def generate_statistics(df: pd.DataFrame) -> dict[str, object]:
    """
    Генерирует базовую статистику по числовым и категориальным признакам.

    Параметры:
        df (pd.DataFrame): Исходные данные.

    Возвращает:
        dict: Статистика по числовым и категориальным столбцам.
    """
    stats = {}
    try:
        stats["numeric"] = df.describe().to_dict()
        categorical_cols = df.select_dtypes(include=["object", "category"]).columns
        stats["categorical"] = {
            col: df[col].value_counts().to_dict() for col in categorical_cols
        }
        return stats
    except Exception as e:
        print(f"Error generating statistics: {e}")
        return stats


def correlation_analysis(df: pd.DataFrame, target_column: str) -> None:
    """
    Выполняет корреляционный анализ и строит тепловую карту корреляций по DataFrame.

    Параметры:
        df (pd.DataFrame): Данные для анализа.
        target_column (str): Не используется (оставлено для совместимости).
    """
    try:
        corr_matrix = df.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
        plt.title("Correlation Matrix")
        plt.show()
    except Exception as e:
        print(f"Error in correlation analysis: {e}")


def extract_pca(df: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """
    Извлекает главные компоненты (PCA) из числовых данных.

    Параметры:
        df (pd.DataFrame): Входные данные.
        n_components (int): Количество главных компонент.

    Возвращает:
        pd.DataFrame: Главные компоненты (PC1, PC2 и т.д.)

    Исключения:
        ValueError: Если данных недостаточно для PCA.
    """
    numeric_cols = df.select_dtypes(include=["float64", "int"]).columns
    X = df[numeric_cols]
    if X.shape[0] < n_components or X.shape[1] < n_components:
        raise ValueError(
            f"Недостаточно данных для извлечения {n_components} главных компонентов."
        )
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(X)
    pc_columns = [f"PC{i+1}" for i in range(n_components)]
    return pd.DataFrame(data=principal_components, columns=pc_columns)


class LabelEncoderWrapper(LabelEncoder):
    """
    Обёртка для LabelEncoder, совместимая с Pipeline.

    Позволяет корректно работать с категориальными признаками в scikit-learn пайплайнах.
    """

    def fit(self, X, y=None):
        """
        Обучает LabelEncoder на данных X.
        """
        return super().fit(X)

    def transform(self, X):
        """
        Преобразует категориальные признаки в числовые метки.

        Возвращает:
            np.ndarray: Массив преобразованных значений, форма (n_samples, 1)
        """
        return super().transform(X).reshape(-1, 1)

    def fit_transform(self, X, y=None):
        """
        Обучает и преобразует данные X.
        """
        return super().fit_transform(X).reshape(-1, 1)


def select_features(
    df: pd.DataFrame, target_column: str, k: int = 10, encoding_method: str = "onehot"
) -> list[str]:
    """
    Выбирает k лучших признаков для предсказания целевой переменной.

    Использует SelectKBest (f_classif) и предварительную обработку категориальных признаков.

    Параметры:
        df (pd.DataFrame): Исходные данные.
        target_column (str): Имя целевого столбца.
        k (int): Количество лучших признаков.
        encoding_method (str): Метод кодирования категориальных признаков ('onehot' или 'label').

    Возвращает:
        list: Список имён выбранных признаков.

    Исключения:
        ValueError: Если не хватает классов в целевой переменной или не поддерживается метод кодирования.
    """
    X = df.drop(columns=[target_column])
    y = df[target_column]
    # Удаляем признаки с нулевой дисперсией
    X = X.loc[:, X.var() > 0]
    if len(set(y)) < 2:
        raise ValueError("Целевая переменная должна содержать как минимум два класса.")
    categorical_cols = X.select_dtypes(include=["object", "category"]).columns.tolist()
    numeric_cols = X.select_dtypes(include=["float64", "int"]).columns.tolist()
    if encoding_method == "onehot":
        categorical_transformer = OneHotEncoder(drop="first", handle_unknown="ignore")
    elif encoding_method == "label":
        categorical_transformer = Pipeline(
            steps=[("label_encoder", LabelEncoderWrapper())]
        )
    else:
        raise ValueError("Unsupported encoding method. Choose 'onehot' or 'label'.")
    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_cols),
            ("cat", categorical_transformer, categorical_cols),
        ]
    )
    pipeline = Pipeline(
        steps=[
            ("preprocessor", preprocessor),
            ("variance_threshold", VarianceThreshold(threshold=0.0)),
            ("feature_selection", SelectKBest(score_func=f_classif, k=k)),
        ]
    )
    if y.dtype == "object" or y.dtype.name == "category":
        le = LabelEncoder()
        y = le.fit_transform(y)
    try:
        pipeline.fit(X, y)
    except ValueError as e:
        raise ValueError(f"Ошибка при обучении пайплайна: {e}")
    feature_names = get_feature_names(
        preprocessor, X, numeric_cols, categorical_cols, encoding_method
    )
    selector = pipeline.named_steps["feature_selection"]
    mask = selector.get_support()
    selected_features = [feature for feature, m in zip(feature_names, mask) if m]
    return selected_features


def get_feature_names(
    preprocessor,
    X: pd.DataFrame,
    numeric_cols: list[str],
    categorical_cols: list[str],
    encoding_method: str,
) -> list[str]:
    """
    Возвращает имена признаков после преобразования ColumnTransformer.

    Параметры:
        preprocessor (ColumnTransformer): Трансформер для обработки признаков.
        X (pd.DataFrame): Исходные данные.
        numeric_cols (list): Список числовых столбцов.
        categorical_cols (list): Список категориальных столбцов.
        encoding_method (str): Метод кодирования ('onehot' или 'label').

    Возвращает:
        list: Список имён признаков после преобразования.
    """
    feature_names = []
    try:
        feature_names.extend(numeric_cols)
        if encoding_method == "onehot":
            for name, transformer, cols in preprocessor.transformers:
                if name == "cat":
                    if not hasattr(transformer, "categories_"):
                        transformer.fit(X[cols])
                    ohe_feature_names = transformer.get_feature_names_out(cols)
                    feature_names.extend(ohe_feature_names)
                    break
        elif encoding_method == "label":
            feature_names.extend(categorical_cols)
        return feature_names
    except Exception as e:
        print(f"Error in get_feature_names: {e}")
        return feature_names
