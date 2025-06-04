#   
import os
import json
import gffutils
import pandas as pd

def validate_gff_file(file_path, lines_to_check=10):
    """
    Быстрая проверка GFF (или GTF) на формальность:
    - Минимум 9 столбцов
    - # комментарии игнорируются
    Возвращает True, если файл похож на GFF, иначе False и причину.
    """
    try:
        with open(file_path, 'r') as f:
            checked = 0
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 9:
                    return False, f"Строка имеет меньше 9 полей: {line}"
                checked += 1
                if checked >= lines_to_check:
                    break
        if checked == 0:
            return False, "Файл не содержит данных (все строки пустые или комментарии)."
        return True, "OK"
    except Exception as e:
        return False, f"Ошибка при чтении файла: {e}"


def load_gff_advanced(
    file_path,
    dbfn=':memory:',
    force=True,
    filter_types=None,
    parse_attrs=False,
    disk_mode=False
):
    """
    Расширенная версия загрузки GFF с дополнительными функциями:
    1. Валидация на корректность формата (минимальная).
    2. Возможность сохранять gffutils-базу на диск (dbfn != ':memory:').
    3. Фильтрация по типам фич (filter_types) – напр. ['gene','mRNA'].
    4. Распаковка атрибутов (parse_attrs=True) – переносит ключи атрибутов в отдельные колонки.

    Параметры:
    - file_path (str): путь к GFF/GTФ файлу.
    - dbfn (str): путь к базе для gffutils; по умолчанию ':memory:'.
    - force (bool): пересоздавать базу, если уже существует.
    - filter_types (list|None): список типов фич (featuretype) для отбора. Если None – берём все.
    - parse_attrs (bool): если True, распаковываем словарь feature.attributes в колонки DataFrame.
    - disk_mode (bool): если True, после обработки не удаляем базу на диске (dbfn).
    Если False и dbfn != ':memory:', временный файл будет удалён.

    Возвращает:
        pd.DataFrame или None при ошибке.
    """
    # Шаг 1. Проверка GFF
    ok, msg = validate_gff_file(file_path)
    if not ok:
        print(f"[load_gff_advanced] Предупреждение: {msg}")
        # можно вернуть None или продолжить
        # вернём None, если не прошли проверку
        return None

    # Если пользователь не хочет хранить в памяти, dbfn != ':memory:'
    if dbfn == ':memory:' and disk_mode:
        # чтобы не путать пользователя: если disk_mode=True, но dbfn не задан, создадим временный
        dbfn = 'temp_gff_db.db'

    try:
        db = gffutils.create_db(
            file_path,
            dbfn=dbfn,
            force=force,
            keep_order=True
        )
    except Exception as e:
        print(f"[load_gff_advanced] Ошибка при создании gffutils DB: {e}")
        return None

    features = []
    try:
        # Итерируем по всем фичам
        for feature in db.all_features():
            if filter_types and feature.featuretype not in filter_types:
                # Пропускаем ненужные типы
                continue

            row = {
                "seqid": feature.seqid,
                "source": feature.source,
                "type": feature.featuretype,
                "start": feature.start,
                "end": feature.end,
                "score": feature.score,
                "strand": feature.strand,
            }
            if not parse_attrs:
                # Просто складываем всё в JSON
                row["attributes"] = json.dumps(feature.attributes)
            else:
                # Распаковываем наиболее частые ключи
                for key, val in feature.attributes.items():
                    # в GFF атрибуты обычно списки
                    if isinstance(val, list):
                        if len(val) == 1:
                            row[key] = val[0]
                        else:
                            row[key] = ",".join(val)
                    else:
                        row[key] = val
            features.append(row)

        df = pd.DataFrame(features)
        return df

    except Exception as e:
        print(f"[load_gff_advanced] Ошибка при извлечении feature: {e}")
        return None
    finally:
        # Если хотим удалить базу, когда закончили, и user не хочет хранить её
        if dbfn != ':memory:' and not disk_mode:
            try:
                os.remove(dbfn)
            except Exception:
                pass


def chunk_read_gff(file_path, chunk_size=200000, filter_types=None):
    """
    (Примерный) Chunk-based парсер GFF «вручную» без gffutils, для ОЧЕНЬ больших файлов.
    *Внимание*: gffutils не предназначен для "частичных" баз, поэтому тут мы
    просто читаем строки, фильтруем их и выдаём кусками. Разметка и атрибуты - на усмотрение.

    Параметры:
    - file_path (str): путь к GFF
    - chunk_size (int): сколько строк (не считая комментариев) за раз возвращать
    - filter_types (list|None): типы фич, которые нас интересуют.
    Yields:
    pd.DataFrame (каждый чанк)
    """
    # Структура GFF (9 полей): seqid, source, type, start, end, score, strand, phase, attributes
    # Упрощённая реализация
    columns = ["seqid","source","type","start","end","score","strand","phase","attributes"]
    buffer = []
    count = 0
    with open(file_path, "r") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields)<9:
                continue  # пропускаем некорректную строку
            # фильтр по type
            ftype = fields[2]
            if filter_types and ftype not in filter_types:
                continue
            # Сохраняем
            buffer.append(fields)
            count +=1
            if count>=chunk_size:
                # Выдаём DataFrame
                df_chunk = pd.DataFrame(buffer, columns=columns)
                # Преобразуем типы start/end к int, score к float (если не '.')
                df_chunk["start"] = df_chunk["start"].astype(int, errors="ignore")
                df_chunk["end"]   = df_chunk["end"].astype(int, errors="ignore")
                # score может быть '.', если неизвестно
                df_chunk["score"] = df_chunk["score"].apply(lambda x: float(x) if x!="." else None)
                yield df_chunk
                buffer.clear()
                count=0
    # остались данные?
    if buffer:
        df_chunk = pd.DataFrame(buffer, columns=columns)
        df_chunk["start"] = df_chunk["start"].astype(int, errors="ignore")
        df_chunk["end"]   = df_chunk["end"].astype(int, errors="ignore")
        df_chunk["score"] = df_chunk["score"].apply(lambda x: float(x) if x!="." else None)
        yield df_chunk
