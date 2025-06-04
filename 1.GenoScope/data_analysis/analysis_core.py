# analysis_core.py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif, VarianceThreshold
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

def analyze_data(df):
    """
    Выполняет базовый анализ данных: вывод описательной статистики и построение корреляционной матрицы.
    """
    try:
        print(df.describe())
        sns.heatmap(df.corr(), annot=True, cmap='coolwarm')
        plt.show()
    except Exception as e:
        print(f"Error analyzing data: {e}")

def generate_statistics(df):
    """
    Генерирует базовую статистику для числовых и категориальных признаков.
    """
    stats = {}
    try:
        stats['numeric'] = df.describe().to_dict()
        categorical_cols = df.select_dtypes(include=['object', 'category']).columns
        stats['categorical'] = {col: df[col].value_counts().to_dict() for col in categorical_cols}
        return stats
    except Exception as e:
        print(f"Error generating statistics: {e}")
        return stats

def correlation_analysis(df, target_column):
    """
    Выполняет корреляционный анализ и строит корреляционную матрицу.
    """
    try:
        corr_matrix = df.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
        plt.title("Correlation Matrix")
        plt.show()
    except Exception as e:
        print(f"Error in correlation analysis: {e}")

def extract_pca(df, n_components=2):
    """
    Извлекает главные компоненты из числовых данных.
    """
    numeric_cols = df.select_dtypes(include=['float64', 'int']).columns
    X = df[numeric_cols]
    if X.shape[0] < n_components or X.shape[1] < n_components:
        raise ValueError(f"Недостаточно данных для извлечения {n_components} главных компонентов.")
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(X)
    pc_columns = [f'PC{i+1}' for i in range(n_components)]
    return pd.DataFrame(data=principal_components, columns=pc_columns)

class LabelEncoderWrapper(LabelEncoder):
    """
    Обёртка для LabelEncoder, предназначенная для использования внутри Pipeline.
    """
    def fit(self, X, y=None):
        return super().fit(X)

    def transform(self, X):
        return super().transform(X).reshape(-1, 1)

    def fit_transform(self, X, y=None):
        return super().fit_transform(X).reshape(-1, 1)

def select_features(df, target_column, k=10, encoding_method='onehot'):
    """
    Выбирает k лучших признаков для целевой переменной с использованием SelectKBest.
    """
    X = df.drop(columns=[target_column])
    y = df[target_column]
    # Удаляем признаки с нулевой дисперсией
    X = X.loc[:, X.var() > 0]
    if len(set(y)) < 2:
        raise ValueError("Целевая переменная должна содержать как минимум два класса.")
    categorical_cols = X.select_dtypes(include=['object', 'category']).columns.tolist()
    numeric_cols = X.select_dtypes(include=['float64', 'int']).columns.tolist()
    if encoding_method == 'onehot':
        categorical_transformer = OneHotEncoder(drop='first', handle_unknown='ignore')
    elif encoding_method == 'label':
        categorical_transformer = Pipeline(steps=[('label_encoder', LabelEncoderWrapper())])
    else:
        raise ValueError("Unsupported encoding method. Choose 'onehot' or 'label'.")
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', 'passthrough', numeric_cols),
            ('cat', categorical_transformer, categorical_cols)
        ]
    )
    pipeline = Pipeline(steps=[
        ('preprocessor', preprocessor),
        ('variance_threshold', VarianceThreshold(threshold=0.0)),
        ('feature_selection', SelectKBest(score_func=f_classif, k=k))
    ])
    if y.dtype == 'object' or y.dtype.name == 'category':
        le = LabelEncoder()
        y = le.fit_transform(y)
    try:
        pipeline.fit(X, y)
    except ValueError as e:
        raise ValueError(f"Ошибка при обучении пайплайна: {e}")
    feature_names = get_feature_names(preprocessor, X, numeric_cols, categorical_cols, encoding_method)
    selector = pipeline.named_steps['feature_selection']
    mask = selector.get_support()
    selected_features = [feature for feature, m in zip(feature_names, mask) if m]
    return selected_features

def get_feature_names(preprocessor, X, numeric_cols, categorical_cols, encoding_method):
    """
    Возвращает имена признаков после преобразования ColumnTransformer.
    """
    feature_names = []
    try:
        feature_names.extend(numeric_cols)
        if encoding_method == 'onehot':
            for name, transformer, cols in preprocessor.transformers:
                if name == 'cat':
                    if not hasattr(transformer, 'categories_'):
                        transformer.fit(X[cols])
                    ohe_feature_names = transformer.get_feature_names_out(cols)
                    feature_names.extend(ohe_feature_names)
                    break
        elif encoding_method == 'label':
            feature_names.extend(categorical_cols)
        return feature_names
    except Exception as e:
        print(f"Error in get_feature_names: {e}")
        return feature_names

