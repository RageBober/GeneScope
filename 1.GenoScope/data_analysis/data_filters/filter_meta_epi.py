# filter_meta_epi.py

import pandas as pd

def filter_by_otu_coverage(df, coverage_col="OTU_COUNT", min_coverage=100):
    """
    Фильтрация по количеству ридов (или счётов) на таксон (OTU/ASV).
    Удаляем таксоны, у которых coverage_col < min_coverage.
    
    Параметры:
    - df (pd.DataFrame): DataFrame, где каждая строка — OTU/ASV либо сэмпл-таксон.
    - coverage_col (str): столбец с числом ридов (или подсчётом).
    - min_coverage (int): минимально допустимый уровень для OTU/ASV.
    
    Возвращает:
    - pd.DataFrame (фильтрованный)
    """
    if coverage_col not in df.columns:
        print(f"[filter_by_otu_coverage] Колонка '{coverage_col}' не найдена, пропускаем.")
        return df
    
    return df[df[coverage_col] >= min_coverage]


def filter_rare_taxa(df, taxon_col="TAXON", sample_col="SAMPLE_COUNT", min_samples=2):
    """
    Фильтрация "редких" таксонов: удаляем таксоны, которые встречаются 
    в количестве образцов < min_samples (или любая другая логика).
    
    Предполагается, что каждый таксон в df может встречаться 
    несколько раз (по разным сэмплам), или есть агрегированный столбец (sample_col) 
    — сколько сэмплов содержит данный таксон.
    
    Параметры:
    - df
    - taxon_col (str): столбец с названием таксона (например, 'ASV_ID' или 'OTU_ID')
    - sample_col (str): столбец, где указано, в скольких сэмплах встречается таксон 
                        (или «число образцов, где coverage > 0»).
    - min_samples (int): минимальное число образцов, в которых таксон должен встречаться
    
    Возвращает:
    - pd.DataFrame
    """
    # Если у вас другой формат (каждая строка = (sample, taxon, count)), 
    # то для фильтрации “редких” нужно сгруппировать:
    #   counts = df.groupby(taxon_col)[sample_col].count()  
    # И потом убрать тех, у кого count < min_samples.
    # Ниже - упрощённый вариант, если sample_col прямо хранит число образцов.
    
    if sample_col not in df.columns:
        print(f"[filter_rare_taxa] Колонка '{sample_col}' не найдена, пропускаем.")
        return df
    
    return df[df[sample_col] >= min_samples]


def filter_by_maf_epigenetics(df, maf_col="MAF", min_maf=0.05):
    """
    Фильтрация CpG/вариантов с низкой частотой аллеля (MAF).
    Предположим, что в столбце maf_col — частота минорного аллеля для CpG.
    
    Параметры:
    - df
    - maf_col (str): столбец с MAF (minor allele frequency).
    - min_maf (float): минимальная MAF, чтобы вариант считался релевантным.
    
    Возвращает:
    - pd.DataFrame
    """
    if maf_col not in df.columns:
        print(f"[filter_by_maf_epigenetics] Колонка '{maf_col}' не найдена, пропускаем.")
        return df
    
    return df[df[maf_col] >= min_maf]


def filter_by_cpg_confidence(df, coverage_col="CpG_COV", pvalue_col="CpG_PVALUE",
                             min_coverage=10, max_pvalue=0.05):
    """
    Фильтрация CpG по «уверенности сигнала»:
    - coverage_col: минимальное покрытие CpG (min_coverage)
    - pvalue_col: например, detection p-value (если оно > max_pvalue, 
                  считаем сигнал недостоверным).
    
    Параметры:
    - df
    - coverage_col (str): столбец с покрытием CpG (количество ридов?).
    - pvalue_col (str): столбец с p-значением (detection p-value).
    - min_coverage (int): минимальное покрытие для CpG-сайта
    - max_pvalue (float): если pvalue > max_pvalue, считаем сигнал некачественным
    
    Возвращает:
    - pd.DataFrame
    """
    filtered_df = df
    if coverage_col in filtered_df.columns:
        filtered_df = filtered_df[filtered_df[coverage_col] >= min_coverage]
    else:
        print(f"[filter_by_cpg_confidence] Колонка '{coverage_col}' не найдена, пропускаем coverage-фильтр.")

    if pvalue_col in filtered_df.columns:
        filtered_df = filtered_df[filtered_df[pvalue_col] <= max_pvalue]
    else:
        print(f"[filter_by_cpg_confidence] Колонка '{pvalue_col}' не найдена, пропускаем p-value-фильтр.")
    
    return filtered_df


def apply_meta_epi_filter_pipeline(df, pipeline):
    """
    Применяет цепочку шагов (pipeline) к метагеномике/эпигенетике.
    
    Пример pipeline:
    [
      {"type": "otu_coverage", "coverage_col": "OTU_COUNT", "min_coverage": 100},
      {"type": "rare_taxa", "taxon_col": "TAXON", "sample_col": "SAMPLE_COUNT", "min_samples": 3},
      {"type": "maf_epigenetics", "maf_col": "CpG_MAF", "min_maf": 0.1},
      {"type": "cpg_confidence", "coverage_col": "CpG_COV", "pvalue_col": "CpG_PVALUE",
       "min_coverage": 10, "max_pvalue": 0.05}
    ]
    
    Возвращает:
    - pd.DataFrame
    """
    for step in pipeline:
        step_type = step.get("type")
        try:
            if step_type == "otu_coverage":
                df = filter_by_otu_coverage(
                    df,
                    coverage_col=step.get("coverage_col", "OTU_COUNT"),
                    min_coverage=step.get("min_coverage", 100)
                )
            elif step_type == "rare_taxa":
                df = filter_rare_taxa(
                    df,
                    taxon_col=step.get("taxon_col", "TAXON"),
                    sample_col=step.get("sample_col", "SAMPLE_COUNT"),
                    min_samples=step.get("min_samples", 2)
                )
            elif step_type == "maf_epigenetics":
                df = filter_by_maf_epigenetics(
                    df,
                    maf_col=step.get("maf_col", "MAF"),
                    min_maf=step.get("min_maf", 0.05)
                )
            elif step_type == "cpg_confidence":
                df = filter_by_cpg_confidence(
                    df,
                    coverage_col=step.get("coverage_col", "CpG_COV"),
                    pvalue_col=step.get("pvalue_col", "CpG_PVALUE"),
                    min_coverage=step.get("min_coverage", 10),
                    max_pvalue=step.get("max_pvalue", 0.05)
                )
            else:
                print(f"[apply_meta_epi_filter_pipeline] Неизвестный тип фильтра: {step_type}")
        except Exception as e:
            print(f"[apply_meta_epi_filter_pipeline] Ошибка при шаге '{step_type}': {e}")
            # Можно вывести traceback
            # По желанию прерывать или продолжать
            
    return df
