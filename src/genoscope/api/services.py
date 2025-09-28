from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

# ─────────────────────────────
#  VCF → DataFrame
# ─────────────────────────────


def _read_vcf_to_df(p: Path) -> pd.DataFrame:
    """
    Читает VCF/BCF (включая .vcf.gz) в DataFrame с базовыми колонками.
    Требует cyvcf2.
    """
    try:
        from cyvcf2 import VCF  # type: ignore
    except Exception as e:
        raise RuntimeError("VCF: требуется cyvcf2 (poetry add cyvcf2)") from e

    rows: list[dict[str, Any]] = []
    for v in VCF(str(p)):
        rows.append(
            {
                "CHROM": str(v.CHROM),
                "POS": int(v.POS),
                "REF": v.REF,
                "ALT": (
                    ",".join(v.ALT) if isinstance(v.ALT, (list, tuple)) else str(v.ALT)
                ),
                "QUAL": float(v.QUAL) if v.QUAL is not None else None,
                "AF": _extract_af(getattr(v, "INFO", {})),
            }
        )
    return pd.DataFrame(rows)


def _extract_af(info: dict) -> float | None:
    for key in ("AF", "gnomAD_AF", "ALLELE_FREQ"):
        if key in info:
            val = info[key]
            if isinstance(val, (list, tuple)):
                val = val[0]
            try:
                return float(val)
            except Exception:
                return None
    return None


def load_any_to_df(path: Path) -> pd.DataFrame:
    """
    Универсальная загрузка: CSV/TSV/XLSX/VCF(.gz)/BCF → DataFrame.
    """
    name = path.name.lower()
    ext = path.suffix.lower()

    # VCF/BCF проверяем по имени (важно для .vcf.gz)
    if name.endswith(".vcf") or name.endswith(".vcf.gz") or name.endswith(".bcf"):
        return _read_vcf_to_df(path)

    if ext == ".csv":
        return pd.read_csv(path)
    if ext in (".tsv", ".tab"):
        return pd.read_csv(path, sep="\t")
    if ext in (".xlsx", ".xls"):
        try:
            return pd.read_excel(path)
        except Exception as e:
            raise RuntimeError("XLSX: требуется openpyxl (poetry add openpyxl)") from e

    raise RuntimeError(f"Неизвестный формат файла: {path.suffix}")


# ─────────────────────────────
#  Нормализация колонок
# ─────────────────────────────

STANDARD = ["CHROM", "POS", "REF", "ALT", "GENE", "AF", "QUAL"]

SYNONYMS: dict[str, str] = {
    "chrom": "CHROM",
    "chr": "CHROM",
    "chromosome": "CHROM",
    "position": "POS",
    "pos": "POS",
    "start": "POS",
    "reference": "REF",
    "ref": "REF",
    "alternate": "ALT",
    "alt": "ALT",
    "alt_allele": "ALT",
    "gene": "GENE",
    "gene_symbol": "GENE",
    "symbol": "GENE",
    "af": "AF",
    "allele_frequency": "AF",
    "freq": "AF",
    "qual": "QUAL",
    "quality": "QUAL",
    "q": "QUAL",
}


def normalize_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Приводит имена колонок к стандарту и даёт mapping: исходное → нормализованное.
    Гарантирует наличие STANDARD-колонок (заполняет None, если отсутствуют).
    """
    mapping: dict[str, str] = {}
    new_cols: list[str] = []
    for c in df.columns:
        key = re.sub(r"[^a-z0-9]+", "", str(c).lower())
        std = SYNONYMS.get(key, str(c))
        mapping[str(c)] = std
        new_cols.append(std)

    df = df.copy()
    df.columns = new_cols

    for s in STANDARD:
        if s not in df.columns:
            df[s] = None

    # типизация
    if "POS" in df.columns:
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    if "AF" in df.columns:
        df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    if "QUAL" in df.columns:
        df["QUAL"] = pd.to_numeric(df["QUAL"], errors="coerce")

    return df, mapping


# ─────────────────────────────
#  Пресеты и «карточки доказательств»
# ─────────────────────────────

PRESETS: dict[str, dict[str, Any]] = {
    "rare_monogenic": {
        "name": "Rare monogenic",
        "source": "gnomAD",
        "version": "v4.0",
        "date": "2025-08-11",
        "min_af": None,
        "max_af": 0.01,
        "min_qual": 30.0,
        "genes": None,
        "chroms": None,
    },
    "oncology": {
        "name": "Oncology panel",
        "source": "OncoKB",
        "version": "2025-07",
        "date": "2025-08-11",
        "min_af": 0.01,
        "max_af": None,
        "min_qual": 20.0,
        "genes": None,
        "chroms": None,
    },
    "pharmacogenomics": {
        "name": "Pharmacogenomics",
        "source": "PharmGKB",
        "version": "2025-06",
        "date": "2025-08-11",
        "min_af": None,
        "max_af": 0.05,
        "min_qual": 20.0,
        "genes": None,
        "chroms": None,
    },
}


def build_effective_params(
    user: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """
    Смешивает пользовательские параметры с пресетом (если задан),
    возвращает эффективные параметры + список «карточек» о применённых правилах.
    """
    applied: list[dict[str, Any]] = []
    eff: dict[str, Any] = {k: v for k, v in user.items() if v is not None}

    preset_key = user.get("preset")
    if preset_key:
        p = PRESETS.get(preset_key)
        if p:
            applied.append(
                {
                    "name": p["name"],
                    "source": p["source"],
                    "version": p["version"],
                    "date": p["date"],
                    "details": {
                        k: p[k]
                        for k in ("min_af", "max_af", "min_qual", "genes", "chroms")
                    },
                }
            )
            for k in ("min_af", "max_af", "min_qual", "genes", "chroms"):
                if eff.get(k) is None and p.get(k) is not None:
                    eff[k] = p[k]

    return eff, applied


# ─────────────────────────────
#  Фильтрация
# ─────────────────────────────


def filter_df(df: pd.DataFrame, eff: dict[str, Any]) -> pd.DataFrame:
    """
    Применяет пороги/списки к DataFrame.
    """
    if df.empty:
        return df

    m = pd.Series(True, index=df.index)
    if eff.get("min_af") is not None:
        m &= df["AF"].fillna(-1) >= float(eff["min_af"])
    if eff.get("max_af") is not None:
        m &= df["AF"].fillna(1e9) <= float(eff["max_af"])
    if eff.get("min_qual") is not None:
        m &= df["QUAL"].fillna(-1) >= float(eff["min_qual"])
    if eff.get("genes"):
        genes = set(map(str.upper, eff["genes"]))
        m &= df["GENE"].fillna("").str.upper().isin(genes)
    if eff.get("chroms"):
        chroms = set(map(str, eff["chroms"]))
        m &= df["CHROM"].fillna("").astype(str).isin(chroms)
    return df[m]


# ─────────────────────────────
#  Отчёт HTML → (PDF опционально)
# ─────────────────────────────

HTML_TMPL = """<!doctype html>
<html lang="{{ lang }}">
<head>
  <meta charset="utf-8">
  <title>{{ title }}</title>
  <style>
    body { font-family: system-ui, Arial, sans-serif; margin: 24px; }
    h1  { margin-bottom: 4px; }
    .muted { color: #666; font-size: 12px; }
    table { border-collapse: collapse; width: 100%; margin-top: 16px; }
    th, td { border: 1px solid #ddd; padding: 6px 8px; font-size: 12px; }
    th { background: #fafafa; text-align: left; }
    .cards { margin-top: 14px; }
    .card { border: 1px solid #ddd; padding: 8px; margin: 8px 0; border-radius: 6px; }
  </style>
</head>
<body>
  <div style="display:flex; align-items:center; gap:12px;">
    {% if logo_url %}<img src="{{ logo_url }}" alt="logo" height="36">{% endif %}
    <h1>{{ title }}</h1>
  </div>
  <div class="muted">{{ clinic }} • {{ date }}</div>

  <div class="cards">
    <h3>{{ t['rules'] }}</h3>
    {% for c in cards %}
      <div class="card">
        <b>{{ c.name }}</b> — {{ c.source }} ({{ c.version }}, {{ c.date }})
        {% if c.details %}<pre>{{ c.details }}</pre>{% endif %}
      </div>
    {% endfor %}
  </div>

  <h3>{{ t['summary'] }}: {{ total }}</h3>
  <table>
    <thead><tr>
      {% for col in cols %}<th>{{ col }}</th>{% endfor %}
    </tr></thead>
    <tbody>
      {% for row in rows %}
        <tr>{% for col in cols %}<td>{{ row.get(col) }}</td>{% endfor %}</tr>
      {% endfor %}
    </tbody>
  </table>

  <p class="muted">{{ t['disclaimer'] }}</p>
</body>
</html>"""

I18N = {
    "ru": {
        "rules": "Карточки доказательств",
        "summary": "Найдено вариантов",
        "disclaimer": "Отчёт предназначен для исследовательских целей. Не является диагнозом.",
    },
    "en": {
        "rules": "Evidence cards",
        "summary": "Variants found",
        "disclaimer": "This report is for research purposes only and not a diagnosis.",
    },
}


def make_html_report(
    preview_rows: list[dict[str, Any]],
    total: int,
    cards: list[dict[str, Any]],
    clinic: str | None,
    logo_url: str | None,
    locale: str,
) -> str:
    """
    Рендерит HTML-отчёт. PDF делается в maybe_render_pdf().
    """
    from jinja2 import Template  # локальный импорт, чтобы не тянуть зависимость зря

    t = I18N.get(locale, I18N["en"])
    cols = ["CHROM", "POS", "REF", "ALT", "GENE", "AF", "QUAL"]
    return Template(HTML_TMPL).render(
        title="GenoScope Report",
        clinic=clinic or "Clinic",
        date=datetime.now().strftime("%Y-%m-%d"),
        logo_url=logo_url,
        t=t,
        cols=cols,
        rows=preview_rows,
        total=total,
        cards=cards,
        lang=locale,
    )


def maybe_render_pdf(html: str, out_path: Path, fmt: str = "html") -> Path:
    """
    Сохраняет отчёт: HTML (всегда) или PDF (если установлен weasyprint).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "pdf":
        try:
            from weasyprint import HTML as WHTML  # type: ignore

            WHTML(string=html).write_pdf(str(out_path))
            return out_path
        except Exception:
            # Если нет weasyprint или ошибка рендера — откат к HTML
            fmt = "html"

    # HTML по умолчанию
    out_path.write_text(html, encoding="utf-8")
    return out_path
