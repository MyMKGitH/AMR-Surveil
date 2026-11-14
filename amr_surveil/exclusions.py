"""A single, explicit record of every record excluded from an analysis and why.

Surveillance practice requires a record of excluded samples and the reason for each. Every
function that drops, ignores or defers a record accepts an optional :class:`ExclusionLog`
so that nothing disappears silently.
"""

from __future__ import annotations

import pandas as pd

EXCLUSION_COLUMNS = ["stage", "table", "record_id", "reason", "detail"]


class ExclusionLog:
    """Append-only collector of excluded records."""

    def __init__(self) -> None:
        self._rows: list[dict[str, str]] = []

    def add(self, stage: str, table: str, record_id: object, reason: str, detail: str = "") -> None:
        self._rows.append(
            {
                "stage": stage,
                "table": table,
                "record_id": "" if record_id is None or pd.isna(record_id) else str(record_id),
                "reason": reason,
                "detail": detail,
            }
        )

    def add_many(
        self, stage: str, table: str, record_ids: object, reason: str, detail: str = ""
    ) -> None:
        for record_id in record_ids:  # type: ignore[attr-defined]
            self.add(stage, table, record_id, reason, detail)

    def __len__(self) -> int:
        return len(self._rows)

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame(self._rows, columns=EXCLUSION_COLUMNS)
        if frame.empty:
            return frame
        return (
            frame.drop_duplicates()
            .sort_values(["stage", "table", "reason", "record_id"], kind="stable")
            .reset_index(drop=True)
        )

    def summary(self) -> pd.DataFrame:
        frame = self.to_frame()
        if frame.empty:
            return pd.DataFrame(columns=["stage", "table", "reason", "n_records"])
        return (
            frame.groupby(["stage", "table", "reason"], sort=True)
            .size()
            .reset_index(name="n_records")
        )
