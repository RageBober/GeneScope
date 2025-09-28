from __future__ import annotations

import uuid
from pathlib import Path

import pandas as pd
from fastapi import APIRouter
from fastapi import Depends
from fastapi import FastAPI
from fastapi import File
from fastapi import HTTPException
from fastapi import Query
from fastapi import UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.responses import FileResponse
from fastapi.responses import ORJSONResponse
from fastapi.staticfiles import StaticFiles
from starlette.middleware.gzip import GZipMiddleware
from starlette.responses import RedirectResponse

from genoscope.core.exceptions import install_exception_handlers
from genoscope.core.logging_config import setup_logging
from genoscope.core.security import async_save_upload
from genoscope.core.security import secure_filename
from genoscope.core.security import validate_saved_file
from genoscope.core.settings import settings

# Import metrics if available
try:
    from genoscope.monitoring.metrics import MetricsMiddleware, metrics_endpoint
    METRICS_OK = True
except ImportError:
    METRICS_OK = False
    MetricsMiddleware = None
    metrics_endpoint = None

# ─────────────────────────────
#  Наши модели/схемы/сервисы (новые)
# ─────────────────────────────
from .models_data import Dataset
from .models_data import FilterRun
from .models_data import Report
from .schemas import EvidenceCard
from .schemas import FilterResponse
from .schemas import UploadResponse
from .services import build_effective_params
from .services import filter_df
from .services import load_any_to_df
from .services import make_html_report
from .services import maybe_render_pdf
from .services import normalize_columns
from .ui import router as ui_router

# ─────────────────────────────
#  База данных (опционально)
# ─────────────────────────────
DB_OK = True
try:
    from sqlmodel import Session  # type: ignore
    from sqlmodel import select  # type: ignore

    from .db import get_session  # твой модуль
    from .db import init_db  # твой модуль

    _ = (Dataset, FilterRun, Report)
    try:
        init_db()
    except Exception:
        DB_OK = False
except Exception:
    DB_OK = False
    Session = object  # только для type hints

    def get_session():
        raise HTTPException(503, "DB not configured")


# ─────────────────────────────
#  Старые сущности (Job/Celery) — опционально
# ─────────────────────────────
JOB_OK = False
CELERY_OK = False
try:
    from .models import Job  # таблица Job из старого кода

    JOB_OK = True
except Exception:
    Job = None  # type: ignore[assignment]

try:
    from .celery_app import run_pipeline
    from .celery_app import sleepy

    CELERY_OK = True
except Exception:

    def sleepy(seconds: int):
        return None

    def run_pipeline(*a, **k):
        return None


# ─────────────────────────────
#  Genomics router — опционально
# ─────────────────────────────
GENOMICS_OK = False
genomics_router = None
try:
    from .genomics_router import genomics_router

    GENOMICS_OK = True
except Exception:
    genomics_router = None

# ─────────────────────────────
#  Billing router — опционально
# ─────────────────────────────
BILLING_ROUTER_OK = False
billing_router = None
try:
    from .billing_router import billing_router

    BILLING_ROUTER_OK = True
except Exception:
    billing_router = None

# ─────────────────────────────
#  Pipeline router — опционально
# ─────────────────────────────
PIPELINE_ROUTER_OK = False
pipeline_router = None
try:
    from .pipeline_router import pipeline_router

    PIPELINE_ROUTER_OK = True
except Exception:
    pipeline_router = None


# ─────────────────────────────
#  Приложение
# ─────────────────────────────
# Настройка логирования
setup_logging(log_level="DEBUG" if settings.debug else "INFO")

app = FastAPI(
    title="GenoScope API",
    default_response_class=ORJSONResponse,
    docs_url=None,  # отключаем дефолтные /docs
    openapi_url="/openapi.json",  # оставляем OpenAPI по этому пути
)
install_exception_handlers(app, include_trace=(settings.ENV != "prod"))
router = APIRouter(prefix="")

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.ALLOW_ORIGINS,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(GZipMiddleware, minimum_size=1000)

# Add metrics middleware if available
if METRICS_OK and MetricsMiddleware:
    app.add_middleware(MetricsMiddleware)

# каталоги
UPLOAD_DIR = Path(settings.UPLOAD_DIR)
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = Path("data")
DATA_DIR.mkdir(parents=True, exist_ok=True)
REPORTS_DIR = DATA_DIR / "reports"
REPORTS_DIR.mkdir(parents=True, exist_ok=True)
DATASETS_DIR = DATA_DIR / "datasets"
DATASETS_DIR.mkdir(parents=True, exist_ok=True)


@app.get("/docs", include_in_schema=False)
def custom_docs():
    return get_swagger_ui_html(
        openapi_url=app.openapi_url,
        title=f"{app.title} — Docs",
        swagger_favicon_url="/favicon.svg",  # берём нашу иконку
    )


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Статика / фавиконки / манифест
STATIC_DIR = Path(__file__).with_name("static")
STATIC_DIR.mkdir(parents=True, exist_ok=True)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/favicon.svg", include_in_schema=False)
def favicon_svg():
    return FileResponse(STATIC_DIR / "favicon.svg", media_type="image/svg+xml")


@app.get("/favicon.ico", include_in_schema=False)
def favicon_ico():
    p = STATIC_DIR / "favicon.ico"
    return FileResponse(p if p.exists() else STATIC_DIR / "favicon.svg")


@app.get("/favicon-16x16.png", include_in_schema=False)
def favicon_16():
    p = STATIC_DIR / "favicon-16.png"
    return FileResponse(
        p if p.exists() else STATIC_DIR / "favicon.svg", media_type="image/png"
    )


@app.get("/favicon-32x32.png", include_in_schema=False)
def favicon_32():
    p = STATIC_DIR / "favicon-32.png"
    return FileResponse(
        p if p.exists() else STATIC_DIR / "favicon.svg", media_type="image/png"
    )


@app.get("/site.webmanifest", include_in_schema=False)
def webmanifest():
    return FileResponse(
        STATIC_DIR / "site.webmanifest", media_type="application/manifest+json"
    )


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# отдача загруженных файлов как статику
app.mount("/files", StaticFiles(directory=str(UPLOAD_DIR)), name="uploads")


# ─────────────────────────────
#  helpers
# ─────────────────────────────
def df_to_records(df: pd.DataFrame, limit: int | None = None) -> list[dict]:
    """Безопасная сериализация в JSON: NaN/pd.NA -> None, numpy -> python types."""
    if limit is not None:
        df = df.head(limit)
    df2 = df.where(pd.notnull(df), None)
    return df2.to_dict(orient="records")  # type: ignore[no-any-return]


# ─────────────────────────────
#  HEALTH + DEMO (старые)
# ─────────────────────────────
@app.get("/", include_in_schema=False)
async def root():
    return RedirectResponse(url="/ui")


@app.get("/health")
async def health_check():
    return {"status": "ok", "db": DB_OK, "celery": CELERY_OK}


# Metrics endpoint
if METRICS_OK and metrics_endpoint:
    @app.get("/metrics")
    async def get_metrics(request):
        return await metrics_endpoint(request)


@app.post("/demo")
async def demo_task(seconds: int = 10):
    if not CELERY_OK:
        raise HTTPException(503, "Celery not configured")
    task = sleepy.delay(seconds)
    return {"task_id": task.id}


# ─────────────────────────────
#  СТАРЫЕ: /files, /jobs/{id}, /reports/{id}
# ─────────────────────────────
@router.post("/files")
async def upload_file(
    file: UploadFile = File(...),
    session: Session = Depends(get_session) if DB_OK else None,
):
    safe = secure_filename(file.filename)
    dest = UPLOAD_DIR / f"{uuid.uuid4()}_{safe}"
    await async_save_upload(file, dest, settings.MAX_FILE_SIZE)

    ok, msg = validate_saved_file(
        dest, set(settings.ALLOWED_EXTENSIONS), settings.MAX_FILE_SIZE
    )
    if not ok:
        dest.unlink(missing_ok=True)
        raise HTTPException(status_code=400, detail=msg)
    if DB_OK and JOB_OK:
        job = Job(filename=dest.name, status="PENDING")  # type: ignore[call-arg]
        session.add(job)  # type: ignore[arg-type]
        session.commit()  # type: ignore[call-arg]
        session.refresh(job)  # type: ignore[call-arg]
        status = "SAVED"
        if CELERY_OK:
            run_pipeline.delay(job.id, str(dest))  # type: ignore[union-attr]
            status = "PENDING"
        return {"job_id": job.id, "status": status}  # type: ignore[union-attr]
    return {"saved": True, "path": str(dest), "tip": "DB disabled, file saved only"}


@router.get("/jobs/{job_id}")
def job_status(
    job_id: int,
    session: Session = Depends(get_session) if DB_OK else None,
):
    if not (DB_OK and JOB_OK):
        raise HTTPException(503, "DB not configured")
    job = session.get(Job, job_id)  # type: ignore[arg-type]
    if not job:
        raise HTTPException(404, "Job not found")
    return job


@router.get("/reports/{job_id}")
def download_report(
    job_id: int,
    session: Session = Depends(get_session) if DB_OK else None,
):
    if not (DB_OK and JOB_OK):
        raise HTTPException(503, "DB not configured")
    job = session.get(Job, job_id)  # type: ignore[arg-type]
    if not job:
        raise HTTPException(404, "Job not found")
    if getattr(job, "status", None) != "DONE" or not getattr(job, "result_path", None):
        raise HTTPException(404, "Report not ready")
    report_path = Path(job.result_path)  # type: ignore[arg-type]
    if not report_path.exists():
        raise HTTPException(404, "Report file missing on disk")
    return FileResponse(report_path, media_type="text/html")


# ─────────────────────────────
#  НОВЫЕ: /datasets/upload, /variants/filter, /report/{rid}
# ─────────────────────────────
@router.post("/datasets/upload", response_model=UploadResponse)
async def datasets_upload(
    file: UploadFile = File(...),
    session: Session = Depends(get_session) if DB_OK else None,
):
    if not DB_OK:
        raise HTTPException(503, "DB not configured")

    # Проверка размера файла
    if file.size and file.size > settings.MAX_FILE_SIZE:
        raise HTTPException(
            413, f"File too large: {file.size} > {settings.MAX_FILE_SIZE}"
        )

    # Проверка типа файла
    if not file.filename:
        raise HTTPException(400, "Filename is required")

    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in settings.ALLOWED_EXTENSIONS:
        raise HTTPException(400, f"File type {file_ext} not allowed")

    tmp = None
    try:
        # сохранить временно
        safe = secure_filename(file.filename)
        tmp = UPLOAD_DIR / f"raw_{uuid.uuid4()}_{safe}"
        await async_save_upload(file, tmp, settings.MAX_FILE_SIZE)

        ok, msg = validate_saved_file(
            tmp, set(settings.ALLOWED_EXTENSIONS), settings.MAX_FILE_SIZE
        )
        if not ok:
            raise HTTPException(status_code=400, detail=msg)

        # прочитать → нормализовать
        try:
            df = load_any_to_df(tmp)
            if df is None or df.empty:
                raise HTTPException(400, "File contains no valid data")
        except ValueError as e:
            logger.warning(f"Invalid file format for {safe}: {e}")
            raise HTTPException(400, f"Invalid file format: {str(e)[:100]}")
        except pd.errors.EmptyDataError:
            logger.warning(f"Empty data file: {safe}")
            raise HTTPException(400, "File contains no data")
        except pd.errors.ParserError as e:
            logger.warning(f"File parsing error for {safe}: {e}")
            raise HTTPException(400, f"File format error: {str(e)[:100]}")
        except MemoryError:
            logger.error(f"File too large to process in memory: {safe}")
            raise HTTPException(413, "File too large to process")
        except PermissionError:
            logger.error(f"Permission denied accessing file: {safe}")
            raise HTTPException(500, "File access error")
        except FileNotFoundError:
            logger.error(f"File disappeared during processing: {safe}")
            raise HTTPException(500, "File processing error")
        except Exception as e:
            logger.exception(f"Unexpected error processing file {safe}: {e}")
            raise HTTPException(500, f"File processing failed: {type(e).__name__}")

    except HTTPException:
        if tmp and tmp.exists():
            tmp.unlink(missing_ok=True)
        raise
    except Exception as e:
        if tmp and tmp.exists():
            tmp.unlink(missing_ok=True)
        raise HTTPException(500, f"Upload failed: {e!s}")
    finally:
        if tmp and tmp.exists():
            tmp.unlink(missing_ok=True)

    try:
        df, mapping = normalize_columns(df)
    except Exception as e:
        raise HTTPException(500, f"Data normalization failed: {e!s}")

    # создать Dataset и сохранить Parquet
    ds = Dataset(
        original_name=file.filename,
        path_parquet="",
        n_rows=len(df),
        n_cols=len(df.columns),
        columns=list(map(str, df.columns)),
        mapping=mapping,
    )
    session.add(ds)  # type: ignore[arg-type]
    session.commit()  # type: ignore[call-arg]
    session.refresh(ds)  # type: ignore[call-arg]

    pq_path = DATASETS_DIR / f"dataset_{ds.id}.parquet"
    df.to_parquet(pq_path, index=False)  # требует pyarrow
    ds.path_parquet = str(pq_path)
    session.add(ds)  # type: ignore[arg-type]
    session.commit()  # type: ignore[call-arg]

    head = df_to_records(df, limit=10)
    return UploadResponse(
        dataset_id=ds.id,  # type: ignore[arg-type]
        rows=len(df),
        cols=len(df.columns),
        columns=list(map(str, df.columns)),
        head=head,
        mapping=mapping,
    )


@router.post("/variants/filter", response_model=FilterResponse)
async def variants_filter(
    dataset_id: int,
    preset: str | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    min_qual: float | None = None,
    genes: list[str] | None = Query(None),
    chroms: list[str] | None = Query(None),
    limit: int = 50,
    create_report: bool = False,
    locale: str = "ru",
    clinic: str | None = None,
    logo_url: str | None = None,
    session: Session = Depends(get_session) if DB_OK else None,
):
    # Валидация параметров
    if min_af is not None and (min_af < 0 or min_af > 1):
        raise HTTPException(400, "min_af must be between 0 and 1")

    if max_af is not None and (max_af < 0 or max_af > 1):
        raise HTTPException(400, "max_af must be between 0 and 1")

    if min_af is not None and max_af is not None and min_af >= max_af:
        raise HTTPException(400, "min_af must be less than max_af")

    if min_qual is not None and min_qual < 0:
        raise HTTPException(400, "min_qual must be non-negative")

    if limit <= 0 or limit > 10000:
        raise HTTPException(400, "limit must be between 1 and 10000")

    if locale not in ["ru", "en"]:
        raise HTTPException(400, "locale must be 'ru' or 'en'")

    if dataset_id <= 0:
        raise HTTPException(400, "dataset_id must be positive")
    if not DB_OK:
        raise HTTPException(503, "DB not configured")

    ds = session.get(Dataset, dataset_id)  # type: ignore[arg-type]
    if not ds or not Path(ds.path_parquet).exists():  # type: ignore[union-attr]
        raise HTTPException(404, "Dataset not found")

    df = pd.read_parquet(ds.path_parquet)  # type: ignore[arg-type]
    user = dict(
        dataset_id=dataset_id,
        preset=preset,
        min_af=min_af,
        max_af=max_af,
        min_qual=min_qual,
        genes=genes,
        chroms=chroms,
    )
    eff, cards_dicts = build_effective_params(user)
    dff = filter_df(df, eff)
    total = len(dff)
    preview = df_to_records(dff, limit=limit)

    fr = FilterRun(dataset_id=ds.id, params=eff, preset=preset, total=total)  # type: ignore[arg-type]
    session.add(fr)  # type: ignore[arg-type]
    session.commit()  # type: ignore[call-arg]
    session.refresh(fr)  # type: ignore[call-arg]

    report_id = None
    if create_report:
        html = make_html_report(
            preview_rows=df_to_records(dff, limit=200),
            total=total,
            cards=cards_dicts,
            clinic=clinic,
            logo_url=logo_url,
            locale=locale,
        )
        rpath = REPORTS_DIR / f"report_{fr.id}.html"
        out = maybe_render_pdf(html, rpath, fmt="html")
        rep = Report(filter_run_id=fr.id, path=str(out), fmt="html")
        session.add(rep)  # type: ignore[arg-type]
        session.commit()  # type: ignore[call-arg]
        report_id = rep.id  # type: ignore[assignment]

    cards = [EvidenceCard(**c) for c in cards_dicts]
    return FilterResponse(
        filter_run_id=fr.id,
        total=total,
        preview=preview,
        applied=cards,
        report_id=report_id,
    )


@app.get("/report/{rid}")
def get_report(
    rid: int,
    session: Session = Depends(get_session) if DB_OK else None,
):
    if not DB_OK:
        raise HTTPException(503, "DB not configured")
    rep = session.get(Report, rid)  # type: ignore[arg-type]
    if rep is None:
        rep = session.exec(select(Report).where(Report.filter_run_id == rid)).first()  # type: ignore[call-arg]
        if rep is None:
            raise HTTPException(404, "Report not ready")
    p = Path(rep.path)  # type: ignore[arg-type]
    if not p.exists():
        raise HTTPException(404, "Report file missing")
    media = "application/pdf" if rep.fmt == "pdf" else "text/html"  # type: ignore[attr-defined]
    return FileResponse(p, media_type=media, filename=p.name)


# ─────────────────────────────
#  УТИЛИТА ЗАПУСКА
# ─────────────────────────────
def run() -> None:
    import uvicorn

    uvicorn.run("genoscope.api.main:app", host="0.0.0.0", port=8000, reload=False)


# регистрируем все роуты из router
app.include_router(router)
app.include_router(ui_router)

# Include genomics router if available
if GENOMICS_OK and genomics_router:
    app.include_router(genomics_router)
    print("✓ Genomics router loaded")
else:
    print("✗ Genomics router not available")

# Include billing router if available
if BILLING_ROUTER_OK and billing_router:
    app.include_router(billing_router)
    print("✓ Billing router loaded")
else:
    print("✗ Billing router not available")

# Include pipeline router if available
if PIPELINE_ROUTER_OK and pipeline_router:
    app.include_router(pipeline_router)
    print("✓ Pipeline router loaded")
else:
    print("✗ Pipeline router not available")
