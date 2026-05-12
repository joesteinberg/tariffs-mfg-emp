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
import re
from urllib.request import urlopen
import zipfile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TWT_DIR = ROOT / "external" / "how-restrictive-us-trade"
DEFAULT_OUTPUT = ROOT / "programs" / "python" / "output" / "trump2025_tariffs.csv"
DEFAULT_DIAGNOSTICS = ROOT / "programs" / "python" / "output" / "trump2025_tariff_diagnostics.csv"
DEFAULT_HEADER = ROOT / "programs" / "c" / "src" / "trump2025_tariffs.h"
DEFAULT_CONCORDANCE_DIR = ROOT / "data" / "concordances"
DEFAULT_H6_GTAP_CONCORDANCE = DEFAULT_CONCORDANCE_DIR / "Concordance_H6_to_GP.zip"
DEFAULT_H6_H3_CONCORDANCE = DEFAULT_CONCORDANCE_DIR / "Concordance_H6_to_H3.zip"
DEFAULT_HS2007_ISIC4_DO = DEFAULT_CONCORDANCE_DIR / "hs2007_to_isic4_roy_zenodo.do"
WITS_H6_GTAP_URL = "https://wits.worldbank.org/data/public/concordance/Concordance_H6_to_GP.zip"
WITS_H6_H3_URL = "https://wits.worldbank.org/data/public/concordance/Concordance_H6_to_H3.zip"
HS2007_ISIC4_URL = (
    "https://zenodo.org/records/16901880/files/"
    "6%20digit%20HS%202007%20to%202%20digit%20ISIC%20Rev%204%20coding.do?download=1"
)

MODEL_SECTORS = ["UPH", "UPL", "DNH", "DNL"]
MODEL_SECTOR_NAMES = {
    "UPH": "Oil",
    "UPL": "Steel",
    "DNH": "Toys",
    "DNL": "Cars",
}

# Constructed from the final OECD ICIO industry clustering in icio.py. The
# industry membership is reported in programs/python/output/sectors.tex and uses
# the ICIO codes from data/industry_names_elasts.csv.
ICIO_TO_MODEL = {
    "A01_02": "UPH",
    "A03": "DNH",
    "B05_06": "UPH",
    "B07_08": "UPH",
    "B09": "UPH",
    "C10T12": "DNL",
    "C13T15": "DNH",
    "C16": "UPH",
    "C17_18": "UPH",
    "C19": "UPH",
    "C20": "UPL",
    "C21": "DNL",
    "C22": "UPL",
    "C23": "UPL",
    "C24": "UPL",
    "C25": "UPH",
    "C26": "DNH",
    "C27": "DNH",
    "C28": "DNL",
    "C29": "DNL",
    "C30": "DNL",
    "C31T33": "DNL",
}

# Mechanical OECD ICIO aggregation of ISIC Rev. 4 two-digit divisions. For
# example, ICIO C10T12 pools ISIC 10, 11, and 12; C31T33 pools 31, 32, and 33.
# These ICIO codes are the same row_sector values consumed by icio.py.
ISIC2_TO_ICIO = {
    1: "A01_02",
    2: "A01_02",
    3: "A03",
    5: "B05_06",
    6: "B05_06",
    7: "B07_08",
    8: "B07_08",
    9: "B09",
    10: "C10T12",
    11: "C10T12",
    12: "C10T12",
    13: "C13T15",
    14: "C13T15",
    15: "C13T15",
    16: "C16",
    17: "C17_18",
    18: "C17_18",
    19: "C19",
    20: "C20",
    21: "C21",
    22: "C22",
    23: "C23",
    24: "C24",
    25: "C25",
    26: "C26",
    27: "C27",
    28: "C28",
    29: "C29",
    30: "C30",
    31: "C31T33",
    32: "C31T33",
    33: "C31T33",
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
    parser.add_argument("--concordance-method", choices=["isic", "gtap"], default="isic",
                        help="Sector bridge to use. 'isic' follows the ICIO industries from icio.py.")
    parser.add_argument("--gtap-concordance", type=Path, default=None,
                        help="Optional WITS H6-to-GTAP CSV or ZIP. Downloaded if omitted.")
    parser.add_argument("--hs2022-hs2007-concordance", type=Path, default=None,
                        help="Optional WITS H6-to-H3 CSV or ZIP. Downloaded if omitted.")
    parser.add_argument("--hs2007-isic4-do", type=Path, default=None,
                        help="Optional HS2007-to-ISIC Rev.4 Stata .do file. Downloaded if omitted.")
    parser.add_argument("--target-date", default="2025-12")
    parser.add_argument("--baseline-year", default="2024")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--diagnostics", type=Path, default=DEFAULT_DIAGNOSTICS)
    parser.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    parser.add_argument("--allow-negative", action="store_true",
                        help="Keep negative HS10 tariff changes instead of clipping them to zero.")
    return parser.parse_args()


def read_csv_or_zipped_csv(path: Path | None, url: str) -> pd.DataFrame:
    def read_csv_fallback(handle) -> pd.DataFrame:
        try:
            return pd.read_csv(handle, encoding="utf-8-sig")
        except UnicodeDecodeError:
            if hasattr(handle, "seek"):
                handle.seek(0)
            return pd.read_csv(handle, encoding="latin1")

    if path is None:
        with urlopen(url) as response:
            raw = response.read()
        zf = zipfile.ZipFile(io.BytesIO(raw))
        member = next(name for name in zf.namelist() if name.lower().endswith(".csv"))
        df = read_csv_fallback(zf.open(member))
    elif path.suffix.lower() == ".zip":
        zf = zipfile.ZipFile(path)
        member = next(name for name in zf.namelist() if name.lower().endswith(".csv"))
        df = read_csv_fallback(zf.open(member))
    else:
        df = read_csv_fallback(path)
    return df


def read_text_file_or_url(path: Path | None, url: str) -> str:
    if path is None:
        with urlopen(url) as response:
            return response.read().decode("utf-8-sig")
    return path.read_text(encoding="utf-8-sig")


def add_equal_concordance_weights(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    df = df.drop_duplicates(keys + ["model_sector"]).copy()
    df["concordance_weight"] = 1.0 / df.groupby(keys)["model_sector"].transform("nunique")
    return df


def use_cached_file(path: Path | None, default_path: Path) -> Path | None:
    if path is not None:
        return path
    if default_path.exists():
        return default_path
    return None


def read_wits_gtap_concordance(path: Path | None) -> pd.DataFrame:
    path = use_cached_file(path, DEFAULT_H6_GTAP_CONCORDANCE)
    df = read_csv_or_zipped_csv(path, WITS_H6_GTAP_URL)
    df = df.rename(columns={
        "HS 2022 Product Code": "hs6",
        "HS 2022 Product Description": "hs6_description",
        "GTAP Product Code": "gtap",
        "GTAP Product Description": "gtap_description",
    })
    df["hs6"] = df["hs6"].astype(str).str.zfill(6)
    df["gtap"] = df["gtap_description"].astype(str).str.split(" - ", n=1).str[0]
    df["model_sector"] = df["gtap"].map(GTAP_TO_MODEL)
    df = df[df["model_sector"].isin(MODEL_SECTORS)]
    return add_equal_concordance_weights(
        df[["hs6", "hs6_description", "gtap", "gtap_description", "model_sector"]],
        ["hs6"],
    )


def clean_stata_do_commands(text: str) -> list[str]:
    commands = []
    current = ""
    for raw_line in text.splitlines():
        line = raw_line.strip().replace("{", "").replace("}", "")
        if not line or line.startswith("*"):
            continue
        if current and re.match(r"^(gen|replace|drop|codebook)\b", line):
            commands.append(current)
            current = ""
        if "///" in line:
            line, trailing = line.split("///", 1)
            line = line.strip()
            current = f"{current} {line}".strip()
            if trailing.strip() and not line.endswith(("|", "&")):
                if current:
                    commands.append(current)
                current = ""
            continue
        if "//" in line:
            line = line.split("//", 1)[0].strip()
        current = f"{current} {line}".strip()
        if current:
            commands.append(current)
        current = ""
    if current:
        commands.append(current)
    return commands


def apply_hs2007_isic4_stata_rules(hs6: pd.Series, do_text: str) -> pd.DataFrame:
    out = pd.DataFrame({"hs2007": hs6.astype(str).str.zfill(6)})
    out["commoditycode"] = pd.to_numeric(out["hs2007"], errors="coerce").astype("Int64")
    out["twodigit"] = (out["commoditycode"] // 10000).astype("Int64")
    out["isic2"] = np.nan

    for command in clean_stata_do_commands(do_text):
        assignment = re.match(r"^(?:gen|replace)\s+isic4\s*=\s*(\d+)\s+if\s+(.+)$", command)
        drop = re.match(r"^drop\s+if\s+(.+)$", command)
        if assignment:
            isic2 = int(assignment.group(1))
            condition = assignment.group(2).replace("isic4", "isic2")
            mask = out.eval(condition, engine="python")
            out.loc[mask.fillna(False), "isic2"] = isic2
        elif drop:
            condition = drop.group(1).replace("isic4", "isic2")
            mask = out.eval(condition, engine="python")
            out.loc[mask.fillna(False), "isic2"] = np.nan

    return out[["hs2007", "isic2"]]


def primary_isic_from_hs2007(hs2007: str) -> int | None:
    code = int(hs2007)
    chapter = code // 10000
    if chapter == 1:
        return 1
    if chapter == 3:
        return 3
    if 60110 <= code <= 60499 or 130110 <= code <= 130239:
        return 1
    if 400110 <= code <= 400130 or 410320 <= code <= 410390:
        return 1
    if 270111 <= code <= 270300 or code in {270500, 270900, 271410, 271490}:
        return 5
    return None


def read_isic_concordance(
    h6_h3_path: Path | None,
    hs2007_isic4_do_path: Path | None,
) -> pd.DataFrame:
    h6_h3_path = use_cached_file(h6_h3_path, DEFAULT_H6_H3_CONCORDANCE)
    hs2007_isic4_do_path = use_cached_file(hs2007_isic4_do_path, DEFAULT_HS2007_ISIC4_DO)
    h6_h3 = read_csv_or_zipped_csv(h6_h3_path, WITS_H6_H3_URL).rename(columns={
        "HS 2022 Product Code": "hs6",
        "HS 2022 Product Description": "hs6_description",
        "HS 2007 Product Code": "hs2007",
        "HS 2007 Product Description": "hs2007_description",
    })
    h6_h3["hs6"] = h6_h3["hs6"].astype(str).str.zfill(6)
    h6_h3["hs2007"] = h6_h3["hs2007"].astype(str).str.zfill(6)

    do_text = read_text_file_or_url(hs2007_isic4_do_path, HS2007_ISIC4_URL)
    hs2007_isic = apply_hs2007_isic4_stata_rules(h6_h3["hs2007"].drop_duplicates(), do_text)
    hs2007_isic["fallback_isic2"] = hs2007_isic["hs2007"].apply(primary_isic_from_hs2007)
    hs2007_isic["isic2"] = hs2007_isic["isic2"].fillna(hs2007_isic["fallback_isic2"])
    hs2007_isic["isic2"] = hs2007_isic["isic2"].astype("Int64")
    hs2007_isic["icio_sector"] = hs2007_isic["isic2"].map(ISIC2_TO_ICIO)
    hs2007_isic["model_sector"] = hs2007_isic["icio_sector"].map(ICIO_TO_MODEL)

    df = pd.merge(h6_h3, hs2007_isic, how="left", on="hs2007")
    df = df[df["model_sector"].isin(MODEL_SECTORS)].copy()
    return add_equal_concordance_weights(
        df[[
            "hs6", "hs6_description", "hs2007", "hs2007_description",
            "isic2", "icio_sector", "model_sector",
        ]],
        ["hs6"],
    )


def read_concordance(args: argparse.Namespace) -> pd.DataFrame:
    if args.concordance_method == "gtap":
        return read_wits_gtap_concordance(args.gtap_concordance)
    return read_isic_concordance(args.hs2022_hs2007_concordance, args.hs2007_isic4_do)


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
    baseline_imports_total = df["weight_imports"].sum()
    negative_shock_import_share = (
        df.loc[df["raw_tariff_shock"] < 0.0, "weight_imports"].sum() / baseline_imports_total
    )
    baseline_weighted_tariff = np.average(df["baseline_tariff"], weights=df["weight_imports"])
    target_weighted_tariff = np.average(df["target_tariff"], weights=df["weight_imports"])
    shock_weighted_tariff = np.average(df["tariff_shock"], weights=df["weight_imports"])
    df["hs6"] = df["hs10"].str[:6]
    df = pd.merge(df, concordance, how="left", on="hs6")
    df["source"] = source
    df["allocated_weight_imports"] = df["weight_imports"] * df["concordance_weight"].fillna(0.0)
    df["allocated_target_imports"] = df["target_imports"] * df["concordance_weight"].fillna(0.0)

    mapped = df[df["model_sector"].isin(MODEL_SECTORS)].copy()
    grouped_rows = []
    for model_sector, group in mapped.groupby("model_sector"):
        grouped_rows.append({
            "model_sector": model_sector,
            "tariff_rate": np.average(group["tariff_shock"], weights=group["allocated_weight_imports"]),
            "weight_imports": group["allocated_weight_imports"].sum(),
            "target_imports": group["allocated_target_imports"].sum(),
            "hs10_count": group["hs10"].nunique(),
        })
    grouped = pd.DataFrame(grouped_rows)
    grouped["source"] = source
    grouped["sector_name"] = grouped["model_sector"].map(MODEL_SECTOR_NAMES)

    diagnostics = pd.DataFrame([
        {
            "source": source,
            "metric": "baseline_imports_total",
            "value": baseline_imports_total,
        },
        {
            "source": source,
            "metric": "mapped_imports_total",
            "value": mapped["allocated_weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "unmatched_imports_total",
            "value": baseline_imports_total - mapped["allocated_weight_imports"].sum(),
        },
        {
            "source": source,
            "metric": "unmatched_import_share",
            "value": (
                (baseline_imports_total - mapped["allocated_weight_imports"].sum())
                / baseline_imports_total
            ),
        },
        {
            "source": source,
            "metric": "negative_shock_import_share",
            "value": negative_shock_import_share,
        },
        {
            "source": source,
            "metric": "baseline_weighted_tariff",
            "value": baseline_weighted_tariff,
        },
        {
            "source": source,
            "metric": "target_weighted_tariff",
            "value": target_weighted_tariff,
        },
        {
            "source": source,
            "metric": "shock_weighted_tariff",
            "value": shock_weighted_tariff,
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
        f"   Target date: {metadata['target_date']}; baseline year: {metadata['baseline_year']}.",
        f"   Concordance method: {metadata['concordance_method']}. */",
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

    concordance = read_concordance(args)
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
        "concordance_method": args.concordance_method,
        "twt_commit": git_commit(args.trade_war_tracker_dir),
    })

    print(f"Wrote {args.output}")
    print(f"Wrote {args.diagnostics}")
    print(f"Wrote {args.header}")
    print(tariffs.to_string(index=False))


if __name__ == "__main__":
    main()
