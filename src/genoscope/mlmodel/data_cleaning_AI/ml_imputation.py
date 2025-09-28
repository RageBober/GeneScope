from __future__ import annotations

import logging
from collections.abc import Iterable
from collections.abc import Sequence

import pandas as pd
from genoscope.data_analysis._typing import DF
from sklearn.base import ClassifierMixin
from sklearn.base import RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression

logger = logging.getLogger(__name__)




def fill_missing_with_ml(
    df: DF,
    columns: Sequence[str],
) -> DF:
    """
    Заполняет пропуски машинным обучением:
    * числовые — `RandomForestRegressor`
    * категориальные — `LogisticRegression`
    """
    df_copy: DF = df.copy()

    for col in columns:
        if col not in df_copy.columns:
            continue

        missing_mask: pd.Series[bool] = df_copy[col].isna()
        if not missing_mask.any():
            continue

        known_data: DF = df_copy.loc[~missing_mask]
        missing_data: DF = df_copy.loc[missing_mask]

        if known_data[col].dropna().empty:
            logger.warning(
                "[fill_missing_with_ml] Столбец «%s» полностью пуст, пропуски не заполним.",
                col,
            )
            continue

        # выбор модели
        model: RegressorMixin | ClassifierMixin
        if pd.api.types.is_numeric_dtype(df_copy[col]):
            model = RandomForestRegressor(n_estimators=10, random_state=42)
        else:
            if known_data[col].nunique(dropna=True) <= 1:
                logger.warning(
                    "[fill_missing_with_ml] «%s» имеет ≤1 класс, "
                    "LogisticRegression не обучим.",
                    col,
                )
                continue
            model = LogisticRegression(max_iter=1_000, n_jobs=1)

        # обучающая выборка без NaN-ов в фичах
        features: Iterable[str] = (c for c in df_copy.columns if c != col)
        feature_cols: list[str] = list(features)

        sub_known: DF = known_data.dropna(subset=feature_cols)
        if sub_known.empty:
            logger.warning(
                "[fill_missing_with_ml] Нет строк без NaN в features для «%s».",
                col,
            )
            continue

        X_train: DF = sub_known[feature_cols]
        y_train = sub_known[col]

        try:
            model.fit(X_train, y_train)
        except Exception as fit_err:
            logger.error(
                "[fill_missing_with_ml] Ошибка обучения модели для «%s»: %s",
                col,
                fit_err,
                exc_info=True,
            )
            continue

        sub_missing: DF = missing_data.dropna(subset=feature_cols)
        if sub_missing.empty:
            logger.warning(
                "[fill_missing_with_ml] Нет строк без NaN в features для предикта «%s».",
                col,
            )
            continue

        X_test: DF = sub_missing[feature_cols]

        # cинхронизируем порядок колонок
        common = X_train.columns.intersection(X_test.columns)
        X_train = X_train[common]
        X_test = X_test[common]

        try:
            preds = model.predict(X_test)
        except NotFittedError:
            logger.error(
                "[fill_missing_with_ml] Модель для «%s» не обучена (NotFittedError).",
                col,
                exc_info=True,
            )
            continue
        except Exception as pred_err:
            logger.error(
                "[fill_missing_with_ml] Ошибка предсказания для «%s»: %s",
                col,
                pred_err,
                exc_info=True,
            )
            continue

        df_copy.loc[X_test.index, col] = preds
        logger.info(
            "[fill_missing_with_ml] Пропуски в «%s» заполнены (%s).",
            col,
            type(model).__name__,
        )

    return df_copy
