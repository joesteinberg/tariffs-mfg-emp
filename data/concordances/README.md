# Concordance Inputs

These files preserve the public concordance inputs used by
`programs/python/trump2025_tariffs.py`.

| File | Source | Role | SHA-256 |
|---|---|---|---|
| `Concordance_H6_to_H3.zip` | WITS product concordance, HS 2022 to HS 2007 | Bridges Trade War Tracker HS2022-style HS6 codes to the HS2007 vintage used by the ISIC rules. | `28a1add2f2cb57d2300a83abb5390218a1a4dcd7d1abedee374d29ff0eef1daf` |
| `Concordance_H6_to_GP.zip` | WITS product concordance, HS 2022 to GTAP | Preserved so the old GTAP-based aggregation remains reproducible with `--concordance-method gtap`. | `1539560a86330e417767816d574da4bd1500dda4e82392d803977250f8f2050f` |
| `hs2007_to_isic4_roy_zenodo.do` | Roy, "Harmonized System (HS) 6 Digit to ISIC Rev 4 Concordance", Zenodo, DOI `10.5281/zenodo.16901880` | Stata rules mapping HS2007 six-digit products to ISIC Rev. 4 two-digit industries. | `657cfef1a215e515d65c037343834ce4e4e821732d7ef90ef6c7ed0b8aaf9998` |

The tariff script first looks for these local files. If they are absent, it
falls back to downloading from the public source URLs embedded in the script.
