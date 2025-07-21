# find_todos.py
import os
import re

KEYWORDS = [r"TODO", r"FIXME", r"NotImplementedError"]


def search_todos(root="."):
    regex = re.compile("|".join(KEYWORDS))
    for dirpath, _, files in os.walk(root):
        for fname in files:
            if fname.endswith(".py"):
                path = os.path.join(dirpath, fname)
                with open(path, encoding="utf-8", errors="ignore") as f:
                    for idx, line in enumerate(f, 1):
                        if regex.search(line):
                            print(f"{path}:{idx}: {line.strip()}")


if __name__ == "__main__":
    search_todos(".")
