from sklearn.feature_extraction.text import CountVectorizer
import numpy as np


def encode_sequences(sequences: list[str], k: int = 3) -> tuple[np.ndarray, list[str]]:
    """
    Кодирует последовательности ДНК/РНК в числовые признаки с использованием k-mers.

    Каждый k-mer (подстрока длины k) преобразуется в отдельный признак (столбец).
    Используется CountVectorizer из sklearn для построения матрицы признаков по частоте встречаемости k-mer'ов.

    Параметры:
        sequences (list или np.ndarray): Список строковых последовательностей (ДНК или РНК).
        k (int, optional): Размерность k-mer'а (длина подстроки, по умолчанию 3).

    Возвращает:
        encoded (np.ndarray): Матрица признаков (n_samples x n_kmers), где каждое значение — количество k-mer в последовательности.
        feature_names (np.ndarray): Массив имён признаков (k-mer'ов).
    """
    vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k))
    encoded = vectorizer.fit_transform(sequences).toarray()
    feature_names = vectorizer.get_feature_names_out()
    return encoded, feature_names
