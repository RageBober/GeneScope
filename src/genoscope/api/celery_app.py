import os
from datetime import datetime  # ← добавить
from pathlib import Path  # ← для расширения, опционально

from celery import Celery
from sqlmodel import Session

from genoscope.data_analysis.data_cleaning import handle_missing_values
from genoscope.data_analysis.data_ingestion import load_data

from .db import engine
from .models import Job

broker_url = os.getenv("CELERY_BROKER", "redis://redis:6379/0")

celery_app = Celery("genoscope", broker=broker_url, backend=broker_url)


@celery_app.task(name="genoscope.sleep_demo")
def sleepy(seconds: int = 10):
    import logging
    import time

    logging.info(f"⏳ starting sleepy({seconds})")
    for i in range(seconds):
        logging.info(f"tick {i+1}/{seconds}")
        time.sleep(1)
    logging.info("✅ sleepy done")
    return {"slept": seconds}


@celery_app.task(name="genoscope.pipeline")
def run_pipeline(job_id: int, filepath: str):
    # 1. помечаем job RUNNING
    with Session(engine) as session:
        job = session.get(Job, job_id)
        job.status = "RUNNING"
        session.commit()

    try:
        # 2. pipeline
        file_type = Path(filepath).suffix.lstrip(".")
        df = load_data(filepath, file_type)
        df = handle_missing_values(df)
        # TODO: ещё clean / filter / analyse …

        # 3. сохраняем отчёт-заглушку
        report_path = f"{filepath}.html"
        with open(report_path, "w") as f:
            f.write("<h1>Pipeline finished OK</h1>")

        # 4. DONE
        with Session(engine) as session:
            job = session.get(Job, job_id)
            job.status = "DONE"
            job.result_path = report_path
            job.updated_at = datetime.utcnow()
            session.commit()

    except Exception as exc:
        # 5. ERROR
        with Session(engine) as session:
            job = session.get(Job, job_id)
            job.status = "ERROR"
            job.updated_at = datetime.utcnow()
            session.commit()
        raise exc
