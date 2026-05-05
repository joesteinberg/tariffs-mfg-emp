#!/usr/bin/env python3
"""Construct Trump-2025 sector tariffs for the C model.

The Trade War Tracker import files contain monthly Census import values and
calculated duties at the HS10 level.  This script measures the 2025 tariff shock
as the December 2025 HS10 duty rate less the 2024 HS10 duty/import baseline,
then aggregates that shock to the paper's four goods sectors using 2024 import
weights.
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path
from urllib.request import urlopen
import zipfile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TWT_DIR = ROOT / "external" / "how-restrictive-us-trade"
DEFAULT_OUTPUT = ROOT / "programs" / "python" / "output" / "trump2025_tariffs.csv"
DEFAULT_DIAGNOSTICS = ROOT / "programs" / "python" / "output" / "trump2025_tariff_diagnostics.csv"
DEFAULT_HEADER = ROOT / "programs" / "c" / "src" / "trump2025_tariffs.h"
WITS_H6_GTAP_URL = "https://wits.worldbank.org/data/public/concordance/Concordance_H6_to_GP.zip"

MODEL_SECTORS = ["UPH", "UPL", "DNH", "DNL"]
MODEL_SECTOR_NAMES = {
    "UPH": "Oil",
    "UPL": "Steel",
    "DNH": "Toys",
    "DNL": "Cars",
}

# GTAP goods sectors mapped to the four goods sectors used by icio.py.
GTAP_TO_MODEL = {
    # "Oil": agriculture, mining, wood, paper, refined petroleum, fabricated metals.
    "B_T": "DNL",
    "C_B": "UPH",
    "CHM": "UPL",
    "CMT": "DNL",
    "COA": "UPH",
    "CTL": "UPH",
    "EEQ": "DNH",
    "ELE": "DNH",
    "ELY": "UPH",
    "FMP": "UPH",
    "FRS": "UPH",
    "FSH": "DNH",
    "GAS": "UPH",
    "GDT": "UPH",
    "GRO": "UPH",
    "I_S": "UPL",
    "LEA": "DNH",
    "LUM": "UPH",
    "MIL": "DNL",
    "MVH": "DNL",
    "NFM": "UPL",
    "NMM": "UPL",
    "OAP": "UPH",
    "OCR": "UPH",
    "OFD": "DNL",
    "OIL": "UPH",
    "OME": "DNL",
    "OMF": "DNL",
    "OMT": "DNL",
    "OSD": "UPH",
    "OTN": "DNL",
    "OXT": "UPH",
    "PCR": "DNL",
    "PDR": "UPH",
    "PFB": "UPH",
    "PPP": "UPH",
    "P_C": "UPH",
    "RPP": "UPL",
    "SGR": "DNL",
    "TEX": "DNH",
    "VOL": "DNL",
    "V_F": "UPH",
    "WAP": "DNH",
    "WHT": "UPH",
    "WOL": "UPH",
    "BPH": "DNL",
}


def read_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trade-war-tracker-dir", type=Path, default=DEFAULT_TWT_DIR)
    parser.add_argument("--gtap-concordance", type=Path, default=None,
                        help="Optional WITS H6-to-GTAP CSV or ZIP. Downloaded if omitted.")
    parser.add_argument("--target-date", default="2025-12")
    parser.add_argument("--baseline-year", default="2024")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--diagnostics", type=Path, default=DEFAULT_DIAGNOSTICS)
    parser.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    parser.add_argument("--allow-negative", action="store_true",
                        help="Keep negative HS10 tariff changes instead of clipping them to zero.")
    return parser.parse_args()


def read_wits_concordance(path: Path | None) -> pd.DataFrame:
    if path is None:
        with urlopen(WITS_H6_GTAP_URL) as response:
            raw = response.read()
        zf = zipfile.ZipFile(io.BytesIO(raw))
        member = next(name for name in zf.namelist() if name.lower().endswith(".csv"))
        df = pd.read_csv(zf.open(member))
    elif path.suffix.lower() == ".zip":
        zf = zipfile.ZipFile(path)
        member = next(name for name in zf.namelist() if name.lower().endswith(".csv"))
        df = pd.read_csv(zf.open(member))
    else:
        df = pd.read_csv(path)

    df = df.rename(columns={
        "HS 2022 Product Code": "hs6",
        "HS 2022 Product Description": "hs6_description",
        "GTAP Product Code": "gtap",
        "GTAP Product Description": "gtap_description",
    })
    df["hs6"] = df["hs6"].astype(str).str.zfill(6)
    df["gtap"] = df["gtap_description"].astype(str).str.split(" - ", n=1).str[0]
    df["model_sector"] = df["gtap"].map(GTAP_TO_MODEL)
    return df[["hs6", "hs6_description", "gtap", "gtap_description", "model_sector"]]


def read_import_file(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df = df[["I_COMMODITY", "time", "CON_VAL_MO", "CAL_DUT_MO"]].copy()
    df["hs10"] = df["I_COMMODITY"].astype(str).str.zfill(10)
    df["CON_VAL_MO"] = pd.to_numeric(df["CON_VAL_MO"], errors="coerce").fillna(0.0)
    df["CAL_DUT_MO"] = pd.to_numeric(df["CAL_DUT_MO"], errors="coerce").fillna(0.0)
    return df[["hs10", "time", "CON_VAL_MO", "CAL_DUT_MO"]]


def aggregate_period(df: pd.DataFrame, period: str, is_year: bool) -> pd.DataFrame:
    if is_year:
        mask = df["time"].astype(str).str.startswith(period + "-")
    else:
        mask = df["time"].astype(str).eq(period)
    out = df.loc[mask].groupby("hs10", as_index=False)[["CON_VAL_MO", "CAL_DUT_MO"]].sum()
    return out.rename(columns={"CON_VAL_MO": "imports", "CAL_DUT_MO": "duties"})


def subtract_imports(total: pd.DataFrame, china: pd.DataFrame) -> pd.DataFrame:
    merged = pd.merge(total, china, how="left", on="hs10", suffixes=("_total", "_china")).fillna(0.0)
    out = pd.DataFrame({
        "hs10": merged["hs10"],
        "imports": merged["imports_total"] - merged["imports_china"],
        "duties": merged["duties_total"] - merged["duties_china"],
    })
    out.loc[out["imports"] < 0.0, "imports"] = 0.0
    out.loc[out["duties"] < 0.0, "duties"] = 0.0
    return out


def tariff_rate(imports: pd.Series, duties: pd.Series) -> pd.Series:
    return duties.div(imports).replace([np.inf, -np.inf], np.nan).fillna(0.0)


def build_source_panel(
    source: str,
    current: pd.DataFrame,
    baseline: pd.DataFrame,
    concordance: pd.DataFrame,
    allow_negative: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    current = current.rename(columns={"imports": "target_imports", "duties": "target_duties"})
    baseline = baseline.rename(columns={"imports": "weight_imports", "duties": "baseline_duties"})
    df = pd.merge(baseline, current, how="left", on="hs10").fillna(0.0)
    df = df[df["weight_imports"] > 0.0].copy()
    df["baseline_tariff"] = tariff_rate(df["weight_imports"], df["baseline_duties"])
    df["target_tariff"] = tariff_rate(df["target_imports"], df["target_duties"])
    df["raw_tariff_shock"] = df["target_tariff"] - df["baseline_tariff"]
    df["tariff_shock"] = df["raw_tariff_shock"] if allow_negative else df["raw_tariff_shock"].clip(lower=0.0)
    df["hs6"] = df["hs10"].str[:6]
    df = pd.merge(df, concordance, how="left", on="hs6")
    df["source"] = source

    mapped = df[df["model_sector"].isin(MODEL_SECTORS)].copy()
    grouped_rows = []
    for model_sector, group in mapped.groupby("model_sector"):
        grouped_rows.append({
            "model_sector": model_sector,
            "tariff_rate": np.average(group["tariff_shock"], weights=group["weight_imports"]),
            "weight_imports": group["weight_imports"].sum(),
            "target_imports": group["target_imports"].sum(),
            "hs10_count": group["hs10"].nunique(),
        })
    grouped = pd.DataFrame(grouped_rows)
    grouped["source"] = source
    grouped["sector_name"] = grouped["model_sector"].map(MODEL_SECTOR_NAMES)

    diagnostics = pd.DataFrame([
        {
            "source": source,
            "metric": "baseline_imports_total",
            "value": df["weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "mapped_imports_total",
            "value": mapped["weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "unmatched_imports_total",
            "value": df.loc[~df["model_sector"].isin(MODEL_SECTORS), "weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "unmatched_import_share",
            "value": (
                df.loc[~df["model_sector"].isin(MODEL_SECTORS), "weight_imports"].sum()
                / df["weight_imports"].sum()
            ),
        },
        {
            "source": source,
            "metric": "negative_shock_import_share",
            "value": df.loc[df["raw_tariff_shock"] < 0.0, "weight_imports"].sum() / df["weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "baseline_weighted_tariff",
            "value": np.average(df["baseline_tariff"], weights=df["weight_imports"]),
        },
        {
            "source": source,
            "metric": "target_weighted_tariff",
            "value": np.average(df["target_tariff"], weights=df["weight_imports"]),
        },
        {
            "source": source,
            "metric": "shock_weighted_tariff",
            "value": np.average(df["tariff_shock"], weights=df["weight_imports"]),
        },
    ])

    return grouped, diagnostics


def write_c_header(path: Path, tariffs: pd.DataFrame, metadata: dict[str, str]) -> None:
    rates = {(row.model_sector, row.source): float(row.tariff_rate) for row in tariffs.itertuples()}
    lines = [
        "#ifndef __TRUMP2025_TARIFFS_H__",
        "#define __TRUMP2025_TARIFFS_H__",
        "",
        "/* Generated by programs/python/trump2025_tariffs.py.",
        f"   Source repo commit: {metadata['twt_commit']}",
        f"   Target date: {metadata['target_date']}; baseline year: {metadata['baseline_year']}. */",
        "",
        "static const double trump2025_tau[NS][NC] = {",
    ]
    for sector in MODEL_SECTORS:
        chn = rates.get((sector, "CHN"), 0.0)
        row = rates.get((sector, "ROW"), 0.0)
        lines.append(f"  {{0.0, {chn:.12f}, {row:.12f}}}, /* {sector} */")
    lines.extend([
        "  {0.0, 0.0, 0.0}, /* SVC */",
        "  {0.0, 0.0, 0.0}  /* CNS */",
        "};",
        "",
        "#endif",
        "",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def git_commit(repo: Path) -> str:
    head = repo / ".git" / "HEAD"
    if not head.exists():
        return "unknown"
    ref = head.read_text(encoding="utf-8").strip()
    if ref.startswith("ref: "):
        ref_path = repo / ".git" / ref.split(" ", 1)[1]
        if ref_path.exists():
            return ref_path.read_text(encoding="utf-8").strip()
    return ref


def main() -> None:
    args = read_args()
    imports_dir = args.trade_war_tracker_dir / "data" / "imports-hs10"
    china_file = imports_dir / "5700data-current.parquet"
    total_file = imports_dir / "TOTALdata-current.parquet"
    if not china_file.exists() or not total_file.exists():
        raise FileNotFoundError(
            "Trade War Tracker parquet files not found. Expected "
            f"{china_file} and {total_file}."
        )

    concordance = read_wits_concordance(args.gtap_concordance)
    china = read_import_file(china_file)
    total = read_import_file(total_file)

    china_current = aggregate_period(china, args.target_date, is_year=False)
    china_baseline = aggregate_period(china, args.baseline_year, is_year=True)
    total_current = aggregate_period(total, args.target_date, is_year=False)
    total_baseline = aggregate_period(total, args.baseline_year, is_year=True)
    row_current = subtract_imports(total_current, china_current)
    row_baseline = subtract_imports(total_baseline, china_baseline)

    china_rates, china_diag = build_source_panel(
        "CHN", china_current, china_baseline, concordance, args.allow_negative
    )
    row_rates, row_diag = build_source_panel(
        "ROW", row_current, row_baseline, concordance, args.allow_negative
    )

    tariffs = pd.concat([china_rates, row_rates], ignore_index=True)
    tariffs["tariff_percent"] = 100.0 * tariffs["tariff_rate"]
    tariffs = tariffs[[
        "source", "model_sector", "sector_name", "tariff_rate", "tariff_percent",
        "weight_imports", "target_imports", "hs10_count",
    ]].sort_values(["source", "model_sector"])
    diagnostics = pd.concat([china_diag, row_diag], ignore_index=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tariffs.to_csv(args.output, index=False)
    diagnostics.to_csv(args.diagnostics, index=False)
    write_c_header(args.header, tariffs, {
        "target_date": args.target_date,
        "baseline_year": args.baseline_year,
        "twt_commit": git_commit(args.trade_war_tracker_dir),
    })

    print(f"Wrote {args.output}")
    print(f"Wrote {args.diagnostics}")
    print(f"Wrote {args.header}")
    print(tariffs.to_string(index=False))


if __name__ == "__main__":
    main()
