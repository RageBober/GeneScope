# filter_variants.py


def filter_by_qual(df, qual_col="QUAL", min_qual=30.0):
    """
    Удаляем варианты (строки), где значение QUAL < min_qual.

    Параметры:
    - df (pd.DataFrame): DataFrame, где каждая строка — вариант
    - qual_col (str): столбец с качеством варианта (обычно QUAL в VCF)
    - min_qual (float): минимальное допустимое качество

    Возвращает:
    - pd.DataFrame (отфильтрованный)
    """
    if qual_col not in df.columns:
        print(
            f"[filter_by_qual] Колонка '{qual_col}' не найдена, пропускаем фильтрацию."
        )
        return df

    # Отфильтровываем строки, где QUAL >= min_qual
    return df[df[qual_col] >= min_qual]


def filter_by_depth(df, depth_col="DP", min_depth=10):
    """
    Фильтрация по глубине покрытия (DP). Удаляем варианты, где DP < min_depth.

    Параметры:
    - df: DataFrame
    - depth_col (str): столбец с глубиной (Coverage Depth). В VCF обычно DP.
    - min_depth (int): минимальная глубина

    Возвращает:
    - pd.DataFrame
    """
    if depth_col not in df.columns:
        print(f"[filter_by_depth] Колонка '{depth_col}' не найдена, пропускаем.")
        return df

    return df[df[depth_col] >= min_depth]


def filter_by_strand_bias(df, sb_col="SB", max_bias=0.01):
    """
    Удаляем варианты с сильным strand bias (артефакт секвенирования).
    Предполагается, что sb_col хранит числовую метрику (p-value или ratio).

    Параметры:
    - df
    - sb_col (str): название столбца со значением strand bias (p-значение или ratio).
    - max_bias (float): порог. Если sb_col > max_bias (или <?), фильтруем.
      В зависимости от того, как измеряется bias.

    Возвращает:
    - pd.DataFrame
    """
    if sb_col not in df.columns:
        print(f"[filter_by_strand_bias] Колонка '{sb_col}' не найдена, пропускаем.")
        return df

    # Тут всё зависит от того, как измеряется strand bias.
    # Если это p-value, то нужно df[sb_col] >= max_bias или наоборот.
    # Ниже для примера считаем: чем МЕНЬШЕ sb_col, тем сильнее bias.
    # Или наоборот. Это пример:

    return df[df[sb_col] > max_bias]


def filter_by_allele_frequency(df, af_col="AF", min_af=0.01, max_af=0.99):
    """
    Фильтрация по частоте аллеля. Удаляем слишком редкие (AF < min_af)
    или слишком частые (AF > max_af) варианты.

    Параметры:
    - df
    - af_col (str): столбец с аллельной частотой
    - min_af, max_af (float): допустимый диапазон

    Возвращает:
    - pd.DataFrame
    """
    if af_col not in df.columns:
        print(
            f"[filter_by_allele_frequency] Колонка '{af_col}' не найдена, пропускаем."
        )
        return df

    return df[(df[af_col] >= min_af) & (df[af_col] <= max_af)]


def filter_by_annotation(df, ann_col="ANN", allowed_types=None):
    """
    Фильтрация по аннотации варианта (ANN).
    Например, удалять синонимичные/некодирующие/неизвестные мутации.

    Параметры:
    - df
    - ann_col (str): столбец с аннотацией (в VCF это может быть INFO.ANN,
      где часто строка формата: 'Allele|Annotation|GeneName|...').
    - allowed_types (list|None): список аннотаций, которые хотим оставить.
      Если None, пропускаем.

    Возвращает:
    - pd.DataFrame
    """
    if allowed_types is None or len(allowed_types) == 0:
        # если не заданы allowed_types, ничего не фильтруем
        return df

    if ann_col not in df.columns:
        print(f"[filter_by_annotation] Колонка '{ann_col}' не найдена, пропускаем.")
        return df

    # Пример: 'ANN' может содержать CSV-список аннотаций,
    # или единую строку 'synonymous_variant|...' и т.д.
    # Ниже — очень упрощённый пример.

    def check_annotation(row_val):
        # row_val — строка, которая может содержать ключевое слово
        # допустим, мы ищем, есть ли AnnotationType,
        # и сравниваем с allowed_types (например, ['missense', 'stop_gained'])
        for a_type in allowed_types:
            if a_type in str(row_val):
                return True
        return False

    mask = df[ann_col].apply(check_annotation)
    return df[mask]


def filter_by_db_snp_clinvar(
    df,
    db_col="dbSNP_ID",
    clinvar_col="ClinVar_ID",
    require_dbsnp=True,
    require_clinvar=False,
):
    """
    Фильтрация по базам данных:
    - Если require_dbsnp=True, оставляем только варианты, у которых есть dbSNP_ID != None.
    - Если require_clinvar=True, оставляем только варианты, у которых ClinVar_ID != None.

    Параметры:
    - df
    - db_col (str): столбец с dbSNP ID
    - clinvar_col (str): столбец с ClinVar ID
    - require_dbsnp (bool): нужно ли оставлять только варианты, присутствующие в dbSNP
    - require_clinvar (bool): аналогично, для ClinVar

    Возвращает:
    - pd.DataFrame
    """
    filtered_df = df
    if require_dbsnp:
        if db_col in filtered_df.columns:
            filtered_df = filtered_df[filtered_df[db_col].notna()]
        else:
            print(
                f"[filter_by_db_snp_clinvar] Колонка '{db_col}' не найдена, пропускаем DB SNP фильтр."
            )

    if require_clinvar:
        if clinvar_col in filtered_df.columns:
            filtered_df = filtered_df[filtered_df[clinvar_col].notna()]
        else:
            print(
                f"[filter_by_db_snp_clinvar] Колонка '{clinvar_col}' не найдена, пропускаем ClinVar фильтр."
            )

    return filtered_df


def apply_variant_filter_pipeline(df, pipeline):
    """
    Аналогично apply_filter_pipeline, но для variant calling:
    Применяет список шагов, где каждый шаг — словарь вида:

    {
      "type": "qual",  # или "depth", "strand_bias", "allele_freq", "annotation", "db_snp_clinvar"
      # и другие параметры
    }

    Пример:
      [
        {"type": "qual", "qual_col": "QUAL", "min_qual": 30},
        {"type": "depth", "depth_col": "DP", "min_depth": 10},
        {"type": "allele_freq", "af_col": "AF", "min_af": 0.01, "max_af": 0.95},
      ]

    Возвращает:
    - pd.DataFrame (после всех шагов)
    """
    for step in pipeline:
        step_type = step.get("type")
        try:
            if step_type == "qual":
                df = filter_by_qual(
                    df,
                    qual_col=step.get("qual_col", "QUAL"),
                    min_qual=step.get("min_qual", 30.0),
                )

            elif step_type == "depth":
                df = filter_by_depth(
                    df,
                    depth_col=step.get("depth_col", "DP"),
                    min_depth=step.get("min_depth", 10),
                )

            elif step_type == "strand_bias":
                df = filter_by_strand_bias(
                    df,
                    sb_col=step.get("sb_col", "SB"),
                    max_bias=step.get("max_bias", 0.01),
                )

            elif step_type == "allele_freq":
                df = filter_by_allele_frequency(
                    df,
                    af_col=step.get("af_col", "AF"),
                    min_af=step.get("min_af", 0.01),
                    max_af=step.get("max_af", 0.99),
                )

            elif step_type == "annotation":
                df = filter_by_annotation(
                    df,
                    ann_col=step.get("ann_col", "ANN"),
                    allowed_types=step.get("allowed_types", None),
                )

            elif step_type == "db_snp_clinvar":
                df = filter_by_db_snp_clinvar(
                    df,
                    db_col=step.get("db_col", "dbSNP_ID"),
                    clinvar_col=step.get("clinvar_col", "ClinVar_ID"),
                    require_dbsnp=step.get("require_dbsnp", True),
                    require_clinvar=step.get("require_clinvar", False),
                )

            else:
                print(
                    f"[apply_variant_filter_pipeline] Неизвестный тип фильтра: {step_type}"
                )

        except Exception as e:
            print(f"[apply_variant_filter_pipeline] Ошибка при шаге '{step_type}': {e}")
            # Можно traceback.print_exc() для подробностей
            # По желанию — прерывать pipeline, либо продолжать

    return df


def filter_variants(df, column="IMPACT", impact_type="HIGH"):
    """
    Фильтрует варианты (например, мутации) по значению в столбце column.
    По умолчанию оставляет только те строки, где IMPACT == "HIGH".
    """
    if column not in df.columns:
        print(f"[filter_variants] Колонка '{column}' не найдена, пропускаем фильтр.")
        return df
    return df[df[column] == impact_type]
