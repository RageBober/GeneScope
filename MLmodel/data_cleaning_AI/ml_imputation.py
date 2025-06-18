# ml_imputation.py
import logging

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression

logger = logging.getLogger(__name__)


def fill_missing_with_ml(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """
    Заполняет пропущенные значения ML-моделью:
    - RandomForestRegressor (числовые)
    - LogisticRegression (категориальные)
    """
    df_copy = df.copy()

    for col in columns:
        # 1) Проверка, есть ли столбец
        if col not in df_copy.columns:
            continue

        # 2) Проверка, есть ли пропуски
        missing_mask = df_copy[col].isnull()
        if not missing_mask.any():
            continue  # нет пропусков

        # 3) Делим на известные/пропущенные
        known_data = df_copy.loc[~missing_mask]
        missing_data = df_copy.loc[missing_mask]

        # Если полностью пустой столбец, модель не обучить
        if known_data[col].dropna().empty:
            logger.warning(
                f"[fill_missing_with_ml] Столбец '{col}' полностью пуст, пропуски не заполним."
            )
            continue

        # 4) Выбираем модель в зависимости от типа
        if pd.api.types.is_numeric_dtype(df_copy[col]):
            model = RandomForestRegressor(n_estimators=10, random_state=42)
        else:
            # Если всего один уникальный класс, логрег не обучить
            if known_data[col].nunique() <= 1:
                logger.warning(
                    f"[fill_missing_with_ml] '{col}' имеет <=1 класс, LogisticRegression не обучим."
                )
                continue
            model = LogisticRegression(max_iter=1000)

        # 5) Обучающая выборка (убираем строки, где NaN в других фичах)
        features = df_copy.columns.drop(col)
        sub_known = known_data.dropna(subset=features)
        if sub_known.empty:
            logger.warning(
                f"[fill_missing_with_ml] Нет строк без пропусков в features для '{col}'."
            )
            continue

        X_train = sub_known[features]
        y_train = sub_known[col]

        # 6) Обучаем
        try:
            model.fit(X_train, y_train)
        except Exception as fit_err:
            logger.error(
                f"[fill_missing_with_ml] Ошибка обучения модели для '{col}': {fit_err}",
                exc_info=True,
            )
            continue

        # 7) Предикт (снова убираем NaN в фичах)
        sub_missing = missing_data.dropna(subset=features)
        if sub_missing.empty:
            logger.warning(
                f"[fill_missing_with_ml] Нет строк без NaN в features для предикта '{col}'."
            )
            continue

        X_test = sub_missing[features]
        # Синхронизируем набор столбцов (если что-то отвалилось)
        common_cols = X_train.columns.intersection(X_test.columns)
        X_train = X_train[common_cols]
        X_test = X_test[common_cols]

        # 8) Предсказываем и заполняем
        try:
            preds = model.predict(X_test)
        except NotFittedError:
            logger.error(
                f"[fill_missing_with_ml] Модель для '{col}' не обучена (NotFittedError).",
                exc_info=True,
            )
            continue
        except Exception as pred_err:
            logger.error(
                f"[fill_missing_with_ml] Ошибка предсказания для '{col}': {pred_err}", exc_info=True
            )
            continue

        df_copy.loc[X_test.index, col] = preds
        logger.info(
            f"[fill_missing_with_ml] Пропуски в '{col}' заполнены моделью {type(model).__name__}."
        )

    return df_copy
