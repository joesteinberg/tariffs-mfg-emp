# Tariff Sector Concordance

This note records the preferred concordance for aggregating HS tariff shocks into
the four goods sectors used by the ICIO/C model.

## Problem

The Trade War Tracker import panels are HS10 monthly Census data. The model,
however, is calibrated on OECD ICIO industries that are then clustered in
`programs/python/icio.py` into four goods sectors:

- `UPH` / "Oil": agriculture, mining, wood, paper, refined petroleum, fabricated metals
- `UPL` / "Steel": chemicals, rubber/plastics, minerals, basic metals
- `DNH` / "Toys": fishing, textiles, electronics, electrical equipment
- `DNL` / "Cars": food/beverages, pharmaceuticals, machinery, motor vehicles, other transport, other manufacturing

The original `trump2025_tariffs.py` bridge used WITS HS2022-to-GTAP and then a
hand-coded GTAP-to-model dictionary. That is transparent enough to reproduce, but
it is not the same taxonomy used to build the model sectors.

## Preferred Approach

The tariff aggregation should route through the ICIO industries:

```text
HS10 tariff line
  -> HS6 product code
  -> HS2007 HS6 code
  -> ISIC Rev. 4 two-digit industry
  -> OECD ICIO industry group
  -> model sector
```

The implemented default does this as follows:

1. Truncate the Census HS10 commodity code to HS6.
2. Use the WITS HS2022-to-HS2007 concordance.
3. Apply the published HS2007-to-ISIC Rev. 4 Stata rules from Roy's Zenodo
   concordance (`10.5281/zenodo.16901880`).
4. Add explicit primary-sector fallbacks for HS codes that the manufacturing
   Stata rules drop but that correspond to ICIO goods sectors, mainly live
   animals/agricultural products, fish, and raw fuels.
5. Map ISIC Rev. 4 two-digit sectors into the exact OECD ICIO groups used by
   `icio.py`, then map those groups into `UPH`, `UPL`, `DNH`, and `DNL`.
6. If an HS2022 code maps to more than one model sector through the concordance,
   split the HS10 observation equally across the distinct model sectors. The
   tariff-rate shock is unchanged; only the import weights are allocated.

This is still not perfect: HS product codes identify traded products, while ICIO
industries identify production activities. Many-to-many product-industry links
would ideally be split with concordance weights or import-use weights. In the
absence of such weights in the public concordances used here, equal splitting is
the explicit rule.

## Mapping Tables

The script has two explicit dictionaries for the ISIC route.

`ISIC2_TO_ICIO` is mechanical. It maps ISIC Rev. 4 two-digit industries into the
OECD ICIO industry codes used in `data/industry_names_elasts.csv` and `icio.py`.
For example, ISIC 10, 11, and 12 become ICIO `C10T12`; ISIC 31, 32, and 33
become `C31T33`; ISIC 05 and 06 become `B05_06`.

`ICIO_TO_MODEL` is constructed from the sector aggregation produced by
`programs/python/icio.py` and reported in `programs/python/output/sectors.tex`:

| Model sector | Label | ICIO industries |
|---|---|---|
| `UPH` | Oil | `A01_02`, `B05_06`, `B07_08`, `B09`, `C16`, `C17_18`, `C19`, `C25` |
| `UPL` | Steel | `C20`, `C22`, `C23`, `C24` |
| `DNH` | Toys | `A03`, `C13T15`, `C26`, `C27` |
| `DNL` | Cars | `C10T12`, `C21`, `C28`, `C29`, `C30`, `C31T33` |

These dictionaries are explicit in code because they are small classification
tables, but their inputs are now documented here and in
`data/concordances/README.md`.

## Preserved Inputs

The concordance files accessed for this workflow are stored in
`data/concordances/`:

- `Concordance_H6_to_H3.zip`
- `Concordance_H6_to_GP.zip`
- `hs2007_to_isic4_roy_zenodo.do`

The local copies make the pipeline reproducible without relying on a live
download. The script uses these local files by default when they exist, and only
downloads from the public URLs if a file is missing or a different path is not
provided.

## Why This Is Better Than GTAP

The GTAP bridge introduced an extra taxonomy that was not used in the model
calibration. The ISIC bridge instead anchors the tariff aggregation to the same
industry names and groupings reported by `icio.py` and `industry_names_elasts.csv`.
That makes the tariff shock, upstreamness measure, trade elasticity, and sector
labels internally consistent.

The old GTAP route is still available with:

```bash
python3 programs/python/trump2025_tariffs.py --concordance-method gtap
```

The default route is now:

```bash
python3 programs/python/trump2025_tariffs.py --concordance-method isic
```

The script accepts local concordance files via:

- `--hs2022-hs2007-concordance`
- `--hs2007-isic4-do`
- `--gtap-concordance`

If omitted, the script downloads the public concordance files.
