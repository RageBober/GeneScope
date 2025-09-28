from pathlib import Path

from fastapi import APIRouter
from fastapi.responses import HTMLResponse

router = APIRouter()
_UI_FILE = Path(__file__).with_name("ui.html")


@router.get("/ui", response_class=HTMLResponse, include_in_schema=False)
def ui_page():
    return HTMLResponse(_UI_FILE.read_text(encoding="utf-8"))
