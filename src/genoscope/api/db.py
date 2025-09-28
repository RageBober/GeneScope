# src/genoscope/api/db.py

from sqlmodel import Session
from sqlmodel import SQLModel
from sqlmodel import create_engine

DB_URL = "sqlite:///./genoscope.db"
engine = create_engine(DB_URL, echo=False, connect_args={"check_same_thread": False})


def init_db() -> None:
    SQLModel.metadata.create_all(engine)


def get_session():
    with Session(engine) as session:
        yield session
