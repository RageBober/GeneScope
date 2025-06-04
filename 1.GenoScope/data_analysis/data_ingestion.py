# data_ingestion.py

import os
import pandas as pd
from cyvcf2 import VCF
import pysam
import h5py
import gffutils
import json

# ВАЖНО: убедитесь, что DataIngestionFunc/auto_gff_loader.py 
# использует относительные импорты, если в том же пакете.
from .DataIngestionFunc.auto_gff_loader import auto_load_gff

############################
# Существующие базовые функции
############################

def load_csv(file_path):
    if not os.path.exists(file_path):
        print(f"[load_csv] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_csv] Error: файл '{file_path}' пустой.")
        return None

    try:
        df = pd.read_csv(file_path)
    except pd.errors.EmptyDataError:
        print(f"[load_csv] Error: файл '{file_path}' не содержит данных (EmptyDataError).")
        return None
    except pd.errors.ParserError as e:
        print(f"[load_csv] Parsing error CSV '{file_path}': {e}")
        return None
    except Exception as e:
        print(f"[load_csv] Непредвиденная ошибка чтения CSV '{file_path}': {e}")
        return None

    if df.empty:
        print(f"[load_csv] Предупреждение: файл '{file_path}' прочитан, но DataFrame пуст.")
        return None

    return df


def load_excel(file_path):
    if not os.path.exists(file_path):
        print(f"[load_excel] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_excel] Error: файл '{file_path}' пустой.")
        return None

    try:
        df = pd.read_excel(file_path)
    except ValueError as e:
        # Иногда при чтении Excel возникает ValueError, если формат не подходит
        print(f"[load_excel] ValueError при чтении Excel '{file_path}': {e}")
        return None
    except Exception as e:
        print(f"[load_excel] Непредвиденная ошибка чтения Excel '{file_path}': {e}")
        return None

    if df.empty:
        print(f"[load_excel] Предупреждение: файл '{file_path}' прочитан, но DataFrame пуст.")
        return None

    return df


def load_json(file_path):
    if not os.path.exists(file_path):
        print(f"[load_json] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_json] Error: файл '{file_path}' пустой.")
        return None

    try:
        df = pd.read_json(file_path)
    except ValueError as e:
        # Если JSON некорректен
        print(f"[load_json] Некорректный JSON в '{file_path}': {e}")
        return None
    except Exception as e:
        print(f"[load_json] Непредвиденная ошибка чтения JSON '{file_path}': {e}")
        return None

    if df.empty:
        print(f"[load_json] Предупреждение: JSON '{file_path}' загружен, но DataFrame пуст.")
        return None

    return df


def load_fasta(file_path):
    if not os.path.exists(file_path):
        print(f"[load_fasta] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_fasta] Error: файл '{file_path}' пустой.")
        return None

    sequences = []
    try:
        with open(file_path, 'r') as file:
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                sequences.append(sequence)

        if not sequences:
            print(f"[load_fasta] Предупреждение: файл '{file_path}' прочитан, но не найдено ни одной последовательности.")
            return None

        return pd.DataFrame({"sequence": sequences})
    except Exception as e:
        print(f"[load_fasta] Error loading FASTA '{file_path}': {e}")
        return None


def load_vcf(file_path):
    if not os.path.exists(file_path):
        print(f"[load_vcf] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_vcf] Error: файл '{file_path}' пустой.")
        return None

    try:
        vcf_reader = VCF(file_path)
        records = []
        for record in vcf_reader:
            records.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "ID": record.ID,
                "REF": record.REF,
                "ALT": [str(alt) for alt in record.ALT],
                "QUAL": record.QUAL,
                "FILTER": record.FILTER
            })

        if not records:
            print(f"[load_vcf] Предупреждение: VCF '{file_path}' прочитан, но записей не найдено.")
            return None

        return pd.DataFrame(records)
    except Exception as e:
        print(f"[load_vcf] Error loading VCF '{file_path}': {e}")
        return None


def load_bam(file_path):
    if not os.path.exists(file_path):
        print(f"[load_bam] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_bam] Error: файл '{file_path}' пустой.")
        return None

    try:
        bam_file = pysam.AlignmentFile(file_path, "rb")
        records = []
        for read in bam_file.fetch():
            records.append({
                "QNAME": read.query_name,
                "FLAG": read.flag,
                "RNAME": bam_file.get_reference_name(read.reference_id),
                "POS": read.reference_start,
                "MAPQ": read.mapping_quality,
                "CIGAR": read.cigarstring,
                "SEQ": read.query_sequence,
                "QUAL": read.query_qualities
            })
        bam_file.close()

        if not records:
            print(f"[load_bam] Предупреждение: BAM '{file_path}' прочитан, но записей не найдено.")
            return None

        return pd.DataFrame(records)
    except Exception as e:
        print(f"[load_bam] Error loading BAM '{file_path}': {e}")
        return None


def load_hdf5(file_path):
    if not os.path.exists(file_path):
        print(f"[load_hdf5] Error: файл '{file_path}' не существует.")
        return None

    if os.path.getsize(file_path) == 0:
        print(f"[load_hdf5] Error: файл '{file_path}' пустой.")
        return None

    try:
        with h5py.File(file_path, 'r') as hdf:
            if not hdf.keys():
                print(f"[load_hdf5] Предупреждение: HDF5 '{file_path}' не содержит групп/данных.")
                return None

            data = {key: hdf[key][:] for key in hdf.keys()}
        if not data:
            print(f"[load_hdf5] Предупреждение: HDF5 '{file_path}' прочитан, но данных не найдено.")
            return None

        return pd.DataFrame(data)
    except Exception as e:
        print(f"[load_hdf5] Error loading HDF5 '{file_path}': {e}")
        return None

############################
# Основная функция load_data
############################

def load_data(file_path, file_type):
    """
    Универсальная функция для загрузки данных из файла.
    
    Параметры:
        file_path (str): Путь к файлу.
        file_type (str): Тип файла ('csv', 'excel', 'json', 'fasta', 'vcf', 'bam', 'gff', 'hdf5').
                        Для Excel также поддерживаются 'xls' и 'xlsx'.
    
    Возвращает:
        pd.DataFrame или генератор DataFrame (в случае большого GFF), или None при ошибке.
    """
    file_type = file_type.lower()  # делаем проверку нечувствительной к регистру

    if file_type == 'csv':
        return load_csv(file_path)
    elif file_type in ['excel', 'xls', 'xlsx']:
        return load_excel(file_path)
    elif file_type == 'json':
        return load_json(file_path)
    elif file_type == 'fasta':
        return load_fasta(file_path)
    elif file_type == 'vcf':
        return load_vcf(file_path)
    elif file_type == 'bam':
        return load_bam(file_path)
    elif file_type == 'gff':
        # Используем auto_load_gff для гибкого определения:
        # если файл маленький => DataFrame, если большой => генератор чанков
        return auto_load_gff(file_path)
    elif file_type == 'hdf5':
        return load_hdf5(file_path)
    else:
        print(f"[load_data] Unsupported file type: '{file_type}'")
        return None
