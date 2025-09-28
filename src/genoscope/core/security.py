import os
import re
from collections.abc import Iterable
from pathlib import Path

import aiofiles
from fastapi import HTTPException
from fastapi import UploadFile

try:
    import magic  # python-magic (не обязателен)
except Exception:
    magic = None

SAFE_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def secure_filename(name: str) -> str:
    base = os.path.basename(name or "file")
    base = base.replace(" ", "_")
    base = SAFE_CHARS.sub("_", base)
    if not base or base in (".", ".."):
        base = "file"
    return base[:200]


def sniff_mime(path: Path) -> str | None:
    if magic is None:
        return None
    try:
        m = magic.Magic(mime=True)
        return m.from_file(str(path))
    except Exception:
        return None


def validate_saved_file(
    path: Path, allowed_exts: Iterable[str], max_size: int
) -> tuple[bool, str]:
    if not path.exists():
        return False, "File not found"
    st = path.stat()
    if st.st_size == 0:
        return False, "Empty file"
    if st.st_size > max_size:
        return False, "File too large"
    ext = path.suffix.lower()
    if allowed_exts and ext not in set(allowed_exts):
        return False, f"Extension {ext} is not allowed"
    # MIME-проверка мягкая (научные форматы часто octet-stream)
    _ = sniff_mime(path)
    return True, "ok"


async def async_save_upload(file: UploadFile, dest: Path, max_size: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    bytes_read = 0
    async with aiofiles.open(dest, "wb") as out:
        while True:
            chunk = await file.read(1024 * 1024)
            if not chunk:
                break
            bytes_read += len(chunk)
            if bytes_read > max_size:
                try:
                    await out.flush()
                except Exception:
                    pass
                dest.unlink(missing_ok=True)
                raise HTTPException(status_code=413, detail="File too large")
            await out.write(chunk)
