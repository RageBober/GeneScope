"""
Улучшенные тесты для модуля очистки данных.
"""

import pandas as pd
import pytest

from genoscope.data_analysis.data_cleaning import _check_numeric_columns
from genoscope.data_analysis.data_cleaning import _fill_missing_with_ml
from genoscope.data_analysis.data_cleaning import detect_outliers
from genoscope.data_analysis.data_cleaning import handle_missing_values
from genoscope.data_analysis.data_cleaning import remove_duplicates


class TestRemoveDuplicates:
    """Тесты функции remove_duplicates."""

    def test_remove_duplicates_empty_df(self):
        """Тест с пустым DataFrame."""
        df = pd.DataFrame()
        result = remove_duplicates(df)
        assert result.empty
        assert isinstance(result, pd.DataFrame)

    def test_remove_duplicates_no_duplicates(self):
        """Тест без дубликатов."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = remove_duplicates(df)
        pd.testing.assert_frame_equal(result, df)
        assert len(result) == 3

    def test_remove_duplicates_with_duplicates(self):
        """Тест с дубликатами."""
        df = pd.DataFrame({"A": [1, 2, 2, 3], "B": [4, 5, 5, 6]})
        result = remove_duplicates(df)
        expected = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=[0, 1, 3])
        pd.testing.assert_frame_equal(result, expected)
        assert len(result) == 3

    def test_remove_duplicates_all_duplicates(self):
        """Тест где все строки - дубликаты."""
        df = pd.DataFrame({"A": [1, 1, 1], "B": [2, 2, 2]})
        result = remove_duplicates(df)
        expected = pd.DataFrame({"A": [1], "B": [2]}, index=[0])
        pd.testing.assert_frame_equal(result, expected)
        assert len(result) == 1


class TestHandleMissingValues:
    """Тесты функции handle_missing_values."""

    def test_handle_missing_empty_df(self):
        """Тест с пустым DataFrame."""
        df = pd.DataFrame()
        result = handle_missing_values(df, method="mean")
        assert result.empty

    def test_handle_missing_no_missing(self):
        """Тест без пропущенных значений."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = handle_missing_values(df, method="mean")
        pd.testing.assert_frame_equal(result, df)

    def test_handle_missing_ffill(self):
        """Тест метода ffill."""
        df = pd.DataFrame({"A": [1, None, 3], "B": [None, 5, 6]})
        result = handle_missing_values(df, method="ffill")
        expected = pd.DataFrame({"A": [1.0, 1.0, 3.0], "B": [None, 5.0, 6.0]})
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_bfill(self):
        """Тест метода bfill."""
        df = pd.DataFrame({"A": [1, None, 3], "B": [None, 5, 6]})
        result = handle_missing_values(df, method="bfill")
        expected = pd.DataFrame({"A": [1.0, 3.0, 3.0], "B": [5.0, 5.0, 6.0]})
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_mean(self):
        """Тест метода mean."""
        df = pd.DataFrame({"A": [1, None, 3], "B": [2, 4, None]})
        result = handle_missing_values(df, method="mean")
        expected = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0],  # mean = (1+3)/2 = 2
                "B": [2.0, 4.0, 3.0],  # mean = (2+4)/2 = 3
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_median(self):
        """Тест метода median."""
        df = pd.DataFrame({"A": [1, None, 3, 5], "B": [2, 4, None, 8]})
        result = handle_missing_values(df, method="median")
        expected = pd.DataFrame(
            {
                "A": [1.0, 3.0, 3.0, 5.0],  # median = 3
                "B": [2.0, 4.0, 4.0, 8.0],  # median = 4
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_mode(self):
        """Тест метода mode."""
        df = pd.DataFrame({"A": [1, None, 1, 2], "B": ["a", "b", None, "b"]})
        result = handle_missing_values(df, method="mode")
        expected = pd.DataFrame(
            {
                "A": [1.0, 1.0, 1.0, 2.0],  # mode = 1
                "B": ["a", "b", "b", "b"],  # mode = "b"
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_interpolate(self):
        """Тест метода interpolate."""
        df = pd.DataFrame(
            {
                "A": [1, None, 3],
                "B": [2, None, None],  # не числовая колонка будет игнорироваться
            }
        )
        result = handle_missing_values(df, method="interpolate")
        expected = pd.DataFrame(
            {"A": [1.0, 2.0, 3.0], "B": [2, None, None]}  # interpolated  # unchanged
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_specific_columns(self):
        """Тест с указанием конкретных колонок."""
        df = pd.DataFrame({"A": [1, None, 3], "B": [None, 5, 6], "C": [7, None, 9]})
        result = handle_missing_values(df, method="mean", columns=["A", "C"])
        expected = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0],  # filled
                "B": [None, 5.0, 6.0],  # unchanged
                "C": [7.0, 8.0, 9.0],  # filled
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_nonexistent_columns(self):
        """Тест с несуществующими колонками."""
        df = pd.DataFrame({"A": [1, None, 3], "B": [4, 5, 6]})
        result = handle_missing_values(df, method="mean", columns=["A", "C", "D"])
        # Должно обработать только существующую колонку A
        expected = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [4.0, 5.0, 6.0]})
        pd.testing.assert_frame_equal(result, expected)

    def test_handle_missing_ml_method(self):
        """Тест ML метода."""
        df = pd.DataFrame({"A": [1, None, 3, 4], "B": [2, 3, None, 5]})
        result = handle_missing_values(df, method="ml")
        # Проверяем что пропуски заполнены
        assert not result.isna().any().any()
        assert len(result) == len(df)

    def test_handle_missing_unsupported_method(self):
        """Тест с неподдерживаемым методом."""
        df = pd.DataFrame({"A": [1, None, 3]})
        result = handle_missing_values(df, method="unsupported")
        # При ошибке должен вернуться исходный DataFrame
        pd.testing.assert_frame_equal(result, df)


class TestDetectOutliers:
    """Тесты функции detect_outliers."""

    def test_detect_outliers_empty_df(self):
        """Тест с пустым DataFrame."""
        df = pd.DataFrame()
        result = detect_outliers(df)
        assert result.empty

    def test_detect_outliers_iqr_method(self):
        """Тест IQR метода."""
        # Данные с очевидным выбросом
        df = pd.DataFrame(
            {"A": [1, 2, 3, 4, 100], "B": [10, 20, 30, 40, 50]}  # 100 - выброс
        )
        result = detect_outliers(df, method="iqr")
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (5, 2)
        # Выброс должен быть помечен как True
        assert result.loc[4, "A"]

    def test_detect_outliers_zscore_method(self):
        """Тест Z-score метода."""
        df = pd.DataFrame(
            {"A": [1, 2, 3, 4, 100], "B": [10, 20, 30, 40, 50]}  # 100 - выброс
        )
        result = detect_outliers(df, method="z-score", threshold=2)
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (5, 2)

    def test_detect_outliers_specific_column(self):
        """Тест с конкретной колонкой."""
        df = pd.DataFrame({"A": [1, 2, 3, 4, 100], "B": [10, 20, 30, 40, 50]})
        result = detect_outliers(df, column="A", method="iqr")
        assert isinstance(result, pd.Series)
        assert len(result) == 5
        assert result.iloc[4]  # выброс

    def test_detect_outliers_mahalanobis_method(self):
        """Тест Mahalanobis метода."""
        df = pd.DataFrame({"A": [1, 2, 3, 4, 100], "B": [10, 20, 30, 40, 50]})
        result = detect_outliers(df, method="mahalanobis", threshold=5)
        assert isinstance(result, pd.Series)
        assert len(result) == 5

    def test_detect_outliers_isolation_forest(self):
        """Тест Isolation Forest метода."""
        df = pd.DataFrame({"A": [1, 2, 3, 4, 100], "B": [10, 20, 30, 40, 50]})
        result = detect_outliers(df, method="isolation_forest", threshold=0.1)
        assert isinstance(result, pd.Series)
        assert len(result) == 5

    def test_detect_outliers_unsupported_method(self):
        """Тест с неподдерживаемым методом."""
        df = pd.DataFrame({"A": [1, 2, 3, 4, 5]})
        with pytest.raises(ValueError):
            detect_outliers(df, method="unsupported")


class TestCheckNumericColumns:
    """Тесты функции _check_numeric_columns."""

    def test_check_numeric_columns_all_numeric(self):
        """Тест с полностью числовыми колонками."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [1.1, 2.2, 3.3]})
        result = _check_numeric_columns(df, ["A", "B"], False, "test")
        assert result == ["A", "B"]

    def test_check_numeric_columns_mixed_types(self):
        """Тест со смешанными типами."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": ["a", "b", "c"], "C": [1.1, 2.2, 3.3]})
        result = _check_numeric_columns(df, ["A", "B", "C"], False, "test")
        assert result == ["A", "C"]

    def test_check_numeric_columns_allow_mixed(self):
        """Тест с разрешением смешанных типов."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": ["a", "b", "c"]})
        result = _check_numeric_columns(df, ["A", "B"], True, "test")
        assert result == ["A"]  # B все равно не числовая

    def test_check_numeric_columns_nonexistent(self):
        """Тест с несуществующими колонками."""
        df = pd.DataFrame({"A": [1, 2, 3]})
        result = _check_numeric_columns(df, ["A", "B"], False, "test")
        assert result == ["A"]


class TestFillMissingWithML:
    """Тесты функции _fill_missing_with_ml."""

    def test_fill_missing_with_ml_numeric(self):
        """Тест ML импутации для числовых колонок."""
        df = pd.DataFrame(
            {"A": [1, None, 3, 4], "B": [2, 3, 4, 5], "C": [1, 2, None, 4]}
        )
        result = _fill_missing_with_ml(df, ["A", "C"])
        # Проверяем что пропуски заполнены
        assert not result[["A", "C"]].isna().any().any()
        assert len(result) == len(df)

    def test_fill_missing_with_ml_categorical(self):
        """Тест ML импутации для категориальных колонок."""
        df = pd.DataFrame(
            {"A": ["a", None, "c", "a"], "B": [1, 2, 3, 4], "C": ["x", "y", None, "x"]}
        )
        result = _fill_missing_with_ml(df, ["A", "C"])
        # Проверяем что пропуски заполнены
        assert not result[["A", "C"]].isna().any().any()
        assert len(result) == len(df)

    def test_fill_missing_with_ml_empty_features(self):
        """Тест ML импутации когда нет других признаков."""
        df = pd.DataFrame({"A": [1, None, 3]})
        result = _fill_missing_with_ml(df, ["A"])
        # Должно остаться без изменений (нет признаков для обучения)
        assert result["A"].isna().sum() == 1

    def test_fill_missing_with_ml_no_missing(self):
        """Тест ML импутации без пропущенных значений."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = _fill_missing_with_ml(df, ["A", "B"])
        # Должно остаться без изменений
        pd.testing.assert_frame_equal(result, df)


class TestIntegration:
    """Интеграционные тесты для модуля очистки данных."""

    def test_full_cleaning_pipeline(self):
        """Тест полного pipeline очистки данных."""
        # Создаем "грязные" данные
        df = pd.DataFrame(
            {
                "A": [1, 1, None, 4, 100],  # дубликаты, пропуски, выбросы
                "B": [2, 2, 3, None, 5],  # дубликаты, пропуски
                "C": ["a", "a", "b", "c", "d"],  # дубликаты
            }
        )

        # Шаг 1: Удаление дубликатов
        step1 = remove_duplicates(df)
        assert len(step1) < len(df)

        # Шаг 2: Заполнение пропусков
        step2 = handle_missing_values(step1, method="mean")
        assert step2.isna().sum().sum() < step1.isna().sum().sum()

        # Шаг 3: Детекция выбросов
        outliers = detect_outliers(step2, method="iqr")
        assert isinstance(outliers, pd.DataFrame)

        # Проверяем что pipeline работает
        assert len(step2) <= len(step1) <= len(df)

    def test_robust_error_handling(self):
        """Тест устойчивости к ошибкам."""
        # Тест с различными проблемными случаями
        test_cases = [
            pd.DataFrame(),  # пустой
            pd.DataFrame({"A": []}),  # пустые колонки
            pd.DataFrame({"A": [None, None, None]}),  # только NaN
            pd.DataFrame({"A": ["a", "b", "c"]}),  # только строки
        ]

        for df in test_cases:
            # Все функции должны работать без исключений
            result1 = remove_duplicates(df)
            result2 = handle_missing_values(df, method="mean")
            result3 = detect_outliers(df, method="iqr")

            # Проверяем базовые инварианты
            assert isinstance(result1, pd.DataFrame)
            assert isinstance(result2, pd.DataFrame)
            assert isinstance(result3, (pd.DataFrame, pd.Series))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
