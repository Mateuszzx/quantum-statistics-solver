from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd


def export_figure_dataframe(
    df: pd.DataFrame,
    stem: str,
    figure_id: str,
    source_notebook: str,
    params: dict | None = None,
):
    stem_path = Path(stem)
    stem_path.parent.mkdir(parents=True, exist_ok=True)

    data = df.reset_index() if isinstance(df.index, pd.MultiIndex) or df.index.name else df.copy()

    csv_path = stem_path.with_suffix(".csv")
    json_path = stem_path.with_suffix(".json")

    data.to_csv(csv_path, index=False)

    meta = {
        "figure_id": figure_id,
        "source_notebook": source_notebook,
        "data_file": csv_path.name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rows": int(len(data)),
        "columns": list(map(str, data.columns)),
        "params": params or {},
    }

    json_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return csv_path, json_path

