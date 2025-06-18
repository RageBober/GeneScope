# auto_gff_loader.py
import os

from .gff_enhanced import chunk_read_gff, load_gff_advanced


def auto_load_gff(
    file_path,
    size_threshold_mb=200,
    chunk_size=200000,
    filter_types=None,
    parse_attrs=True,
    dbfn=":memory:",
    disk_mode=True,
):
    """
    Автоматически определяет, как загружать GFF:
    1) Если файл больше size_threshold_mb (по умолчанию 500 МБ),
    используем chunk_read_gff (чанковый парсер).
    2) Иначе используем load_gff_advanced, чтобы получить DataFrame целиком.

    Параметры:
    - file_path (str): Путь к GFF-файлу.
    - size_threshold_mb (int): Порог размера (в МБ), выше которого используем чанковое чтение.
    - chunk_size (int): Сколько строк на один чанк для chunk_read_gff.
    - filter_types (list | None): Список типов (featuretype), которые хотим оставить (например, ['gene','mRNA']).
    - parse_attrs (bool): Если True, раскладываем атрибуты в отдельные столбцы при обычной загрузке.
    - dbfn (str): Путь к базе gffutils при load_gff_advanced, если не ':memory:'.
    - disk_mode (bool): Если True, не удаляем файл базы после чтения; иначе удалим.

    Возвращает:
    - Если файл не превышает порог: pd.DataFrame.
    - Если превышает порог: генератор (yield) DataFrame чанками.
    """
    file_size = os.path.getsize(file_path)
    size_mb = file_size / (1024 * 1024)
    print(f"[auto_load_gff] Файл {file_path} ~ {size_mb:.2f} МБ")

    if size_mb > size_threshold_mb:
        # используем чанковый режим
        print("[auto_load_gff] Файл большой, переходим к chunk_read_gff...")
        return chunk_read_gff(file_path, chunk_size=chunk_size, filter_types=filter_types)
    else:
        # используем load_gff_advanced
        print("[auto_load_gff] Файл относительно небольшой, используем load_gff_advanced...")
        df = load_gff_advanced(
            file_path,
            dbfn=dbfn,
            force=True,
            filter_types=filter_types,
            parse_attrs=parse_attrs,
            disk_mode=disk_mode,
        )
        return df
