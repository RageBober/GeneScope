import tkinter as tk
import types
from tkinter import filedialog, messagebox, ttk

import pandas as pd

from data_analysis.data_cleaning import handle_missing_values, remove_duplicates
from data_analysis.data_ingestion import load_data


class ProgressWindow:
    """
    Окно прогресса для отображения этапов обработки данных в GenoScope.

    Показывает процент выполнения и обновляется при каждом шаге обработки.
    """

    def __init__(self, total_steps: int, title: str = "Обработка данных") -> None:
        """
        Инициализация окна прогресса.

        Аргументы:
            total_steps (int): Всего этапов в процессе.
            title (str): Заголовок окна.
        """
        self.total_steps = total_steps
        self.current_step = 0

        self.root = tk.Toplevel()
        self.root.title(title)
        self.root.geometry("400x100")

        self.label = ttk.Label(self.root, text="Прогресс: 0%")
        self.label.pack(pady=10)

        self.progressbar = ttk.Progressbar(self.root, maximum=100, length=300, mode="determinate")
        self.progressbar.pack(pady=10)
        self.root.update()

    def update(self, step=1):
        """
        Увеличивает текущий прогресс и обновляет отображение в окне.

        Аргументы:
            step (int): На сколько шагов увеличить прогресс (по умолчанию 1).
        """
        self.current_step += step
        percent = int(100 * self.current_step / self.total_steps)
        percent = min(percent, 100)
        self.label.config(text=f"Прогресс: {percent}%")
        self.progressbar["value"] = percent
        self.root.update()

    def close(self):
        self.root.destroy()


class GenoScopeApp:
    """
    Главное окно приложения GenoScope.

    Обеспечивает графический интерфейс для загрузки данных, удаления дубликатов, фильтрации и других операций над таблицей данных.
    """

    def __init__(self, root):
        self.root = root
        self.root.title("GenoScope")
        self.root.geometry("800x600")

        self.data = None
        self._build_ui()

    def _build_ui(self):
        """
        Формирует и размещает все виджеты интерфейса:
        кнопки загрузки, фильтрации, настройки удаления дубликатов и др.
        """
        # Верхняя панель
        top_frame = tk.Frame(self.root)
        top_frame.pack(pady=5)

        tk.Button(top_frame, text="Загрузить файл", command=self._upload_file).pack(
            side=tk.LEFT, padx=5
        )

        # Блок настроек удаления дубликатов
        dup_frame = tk.LabelFrame(self.root, text="Удаление дубликатов")
        dup_frame.pack(fill="x", padx=5, pady=5)

        tk.Label(dup_frame, text="subset (через запятую):").grid(row=0, column=0)
        self.subset_entry = tk.Entry(dup_frame, width=20)
        self.subset_entry.insert(0, "sample_id")
        self.subset_entry.grid(row=0, column=1, padx=4)

        tk.Label(dup_frame, text="priority_col:").grid(row=0, column=2)
        self.priority_entry = tk.Entry(dup_frame, width=15)
        self.priority_entry.insert(0, "score")
        self.priority_entry.grid(row=0, column=3, padx=4)

        tk.Label(dup_frame, text="priority:").grid(row=0, column=4)
        self.priority_var = tk.StringVar(value="max")
        ttk.Combobox(
            dup_frame,
            textvariable=self.priority_var,
            values=["max", "min"],
            width=5,
            state="readonly",
        ).grid(row=0, column=5, padx=4)

        # Блок заполнения пропусков
        miss_frame = tk.LabelFrame(self.root, text="Заполнение пропусков")
        miss_frame.pack(fill="x", padx=5, pady=5)

        tk.Label(miss_frame, text="Метод:").grid(row=0, column=0)
        self.miss_var = tk.StringVar(value="mean")
        ttk.Combobox(
            miss_frame,
            textvariable=self.miss_var,
            values=["ffill", "bfill", "mean", "median", "mode", "interpolate", "ml"],
            width=15,
            state="readonly",
        ).grid(row=0, column=1, padx=4)
        tk.Button(miss_frame, text="Заполнить пропуски", command=self._fill_missing).grid(
            row=0, column=2, padx=4
        )

        # Поле вывода
        self.output = tk.Text(self.root, height=20, width=100)
        self.output.pack(padx=5, pady=5)

    def _upload_file(self):
        """
        Открывает диалог выбора файла и загружает данные в приложение.
        """

        path = filedialog.askopenfilename(
            filetypes=[
                ("CSV файлы", "*.csv"),
                ("FASTA файлы", "*.fasta"),
                ("JSON файлы", "*.json"),
                ("Excel файлы", "*.xlsx;*.xls"),
                ("VCF файлы", "*.vcf"),
                ("BAM файлы", "*.bam"),
                ("GFF файлы", "*.gff"),
                ("HDF5 файлы", "*.hdf5"),
                ("Все файлы", "*.*"),
            ]
        )
        if not path:
            return

        self.output.insert(tk.END, f"\nЗагружаем файл: {path}\n")

        file_type = path.split(".")[-1].lower()
        if file_type in ["xls", "xlsx"]:
            file_type = "excel"
        elif file_type not in [
            "csv",
            "fasta",
            "json",
            "vcf",
            "bam",
            "gff",
            "hdf5",
            "excel",
        ]:
            self.output.insert(tk.END, "Неизвестный формат, предполагаем CSV.\n")
            file_type = "csv"

        try:
            data = load_data(path, file_type)
            if data is None:
                messagebox.showerror("Ошибка", "Не удалось загрузить файл.")
                return

            if isinstance(data, pd.DataFrame):
                self._process_dataframe(data)
            elif isinstance(data, types.GeneratorType):
                self._process_generator(data)
            else:
                messagebox.showerror("Ошибка", "Неверный формат данных.")

        except Exception as e:
            messagebox.showerror("Ошибка", str(e))

    def _process_dataframe(self, df):
        self.output.insert(tk.END, "\nДанные успешно загружены.\n")
        df = self._remove_duplicates(df)
        self.data = df
        self.output.insert(tk.END, f"{df.head()}\n")

    def _process_generator(self, generator):
        self.output.insert(tk.END, "Большой файл: начинаем чанковую обработку...\n")

        progress = ProgressWindow(total_steps=10)  # Эмулируем 10 чанков
        chunk_num = 0

        for chunk in generator:
            chunk_num += 1
            self.output.insert(tk.END, f"Чанк {chunk_num}, строк: {len(chunk)}\n")
            chunk = self._remove_duplicates(chunk)
            progress.update()

        progress.close()
        self.output.insert(tk.END, f"Обработка чанков завершена. Всего чанков: {chunk_num}\n")

    def _remove_duplicates(self, df):
        subset = [s.strip() for s in self.subset_entry.get().split(",") if s.strip()]
        priority_col = self.priority_entry.get().strip() or None
        priority = self.priority_var.get()

        return remove_duplicates(
            df,
            subset=subset if subset else None,
            priority_col=priority_col,
            priority=priority,
        )

    def _fill_missing(self):
        if self.data is None:
            messagebox.showwarning("Нет данных", "Сначала загрузите файл.")
            return

        method = self.miss_var.get()
        self.data = handle_missing_values(self.data, method=method)
        self.output.insert(
            tk.END,
            f"Пропуски заполнены ({method}). Остаток NaN: {self.data.isna().sum().sum()}\n",
        )


if __name__ == "__main__":
    root = tk.Tk()
    app = GenoScopeApp(root)
    root.mainloop()
