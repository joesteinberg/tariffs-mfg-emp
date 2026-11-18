# README for "Tariffs, Manufacturing Employment, and Supply Chains" (Joseph B. Steinberg)

## Overview

The code in this replication package builds a four-region, six-sector input-output representation of the world economy from public data (OECD, 2023; BEA, n.d.), solves the dynamic general-equilibrium model described in Section 2 of the paper under each tariff scenario, and produces every figure and generated table in the paper and its online appendix. It uses Python and C. One main script, `run_replication.sh`, runs all of the steps in order. The model solves are memory-intensive: each needs about 23 GB of RAM, and the 18 solves with a 200-period horizon (Appendix H) need about 90 GB each. Running the 92 standard solves eight at a time and the 200-period solves two at a time (`JOBS=8 OMP_NUM_THREADS=8 JOBS_NT200=2 THREADS_NT200=32 ./run_replication.sh`), a replicator should expect the code to run for about 11 hours on a 64-core workstation with 256 GB of RAM.

## Data Availability and Provenance Statements

### Statement about Rights

- [x] I certify that the author of the manuscript has legitimate access to and permission to use the data used in this manuscript.
- [x] I certify that the author of the manuscript has documented permission to redistribute/publish the data contained within this replication package. Appropriate permissions are documented in the LICENSE.txt file.

### License for Data

The package is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (see `LICENSE.txt`). The third-party data it contains are public data from the OECD, the U.S. Bureau of Economic Analysis, the World Bank, and Zenodo, which remain subject to their providers' terms, including attribution.

### Summary of Availability

- [x] All data **are** publicly available.

### Summary of Data Availability

| Data name | Data files | Location | Provided | Citation |
| --- | --- | --- | --- | --- |
| OECD Inter-Country Input-Output Tables, 2020 table | `2020_SML.csv` | `data/` | Yes | OECD (2023) |
| BEA Components of Value Added by Industry | `ComponentsOfVa.xlsx` | `data/` | Yes | BEA (n.d.) |
| Industry names and Caliendo-Parro trade elasticities | `industry_names_elasts.csv` | `data/` | Yes | Caliendo and Parro (2015) |
| WITS HS 2022-to-HS 2007 concordance | `Concordance_H6_to_H3.zip` | `data/concordances/` | Yes | World Bank (n.d.) |
| HS 2007-to-ISIC Rev. 4 concordance | `hs2007_to_isic4_roy_zenodo.do` | `data/concordances/` | Yes | Roy (2025) |
| Trade War Tracker HS10 U.S. imports and duties | `data/imports-hs10/*.parquet` in the Trade War Tracker repository | `external/how-restrictive-us-trade/` | No (download instructions below) | Trade War Tracker (2026) |

### Details on each Data Source

**OECD Inter-Country Input-Output (ICIO) Tables.** The model is calibrated to the 2020 table of the 2023 release of the OECD ICIO database (OECD, 2023), in the "SML" format, which reports intermediate and final use by country and industry, value added, and gross output in millions of U.S. dollars. It was downloaded from https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html (accessed May 5, 2026).

Datafile: `data/2020_SML.csv`

**BEA Components of Value Added by Industry.** U.S. compensation of employees by industry in 2023 (sheet `UVCT2-A` of the BEA GDP-by-Industry table "Components of Value Added by Industry") gives each model sector's share of U.S. labor compensation, which sets the U.S. capital shares in the model and the sector employment shares in Table 1. It was downloaded from https://apps.bea.gov/iTable/ (accessed May 5, 2026). I added a column, `ICIO Industry Code`, that maps the BEA industries to the ICIO industries. Where one BEA industry covers two ICIO industries that belong to different model sectors (chemicals and pharmaceuticals; agriculture and fishing), `icio.py` splits its compensation in proportion to the two industries' U.S. value added in the ICIO table. BEA data are in the public domain.

Datafile: `data/ComponentsOfVa.xlsx`

**Industry names and trade elasticities.** An author-created file that gives each ICIO goods industry a short name and its trade elasticity from Table 1 of Caliendo and Parro (2015) (99 percent sample estimates). These elasticities and the industries' upstreamness (computed from the ICIO table following Antràs et al., 2012) determine the sector clustering described in Section 3.1 and Appendix B.

Datafile: `data/industry_names_elasts.csv`

**Concordances.** Two public concordances map HS products to ICIO industries for the 2025 tariff schedule: the World Bank WITS HS 2022-to-HS 2007 concordance (World Bank, n.d., https://wits.worldbank.org/product_concordance.html, accessed May 12, 2026) and the HS 2007-to-ISIC Rev. 4 rules of Roy (2025), DOI 10.5281/zenodo.16901880 (accessed May 12, 2026). See `data/concordances/README.md` for file hashes.

Datafiles: `data/concordances/Concordance_H6_to_H3.zip`, `data/concordances/hs2007_to_isic4_roy_zenodo.do`

**Trade War Tracker.** The 2025 tariff schedule (Section 4.8 and Appendix B) is measured from monthly U.S. Census Bureau HS10 import values and calculated duties, compiled by the Trade War Tracker project (Trade War Tracker, 2026) and distributed under the MIT license at https://github.com/tradewartracker/how-restrictive-us-trade. These files are about 350 MB and are not included in this package. The resulting tariff rates by sector and partner region are included in `programs/c/src/trump2025_tariffs.h`, so the download is needed only to rebuild that file. To obtain the data used in the paper (repository commit `6f2bf4af885ca2e0c5a37c6f868d487559228f32`), run from the root of this package:

```bash
git clone https://github.com/tradewartracker/how-restrictive-us-trade external/how-restrictive-us-trade
git -C external/how-restrictive-us-trade checkout 6f2bf4af885ca2e0c5a37c6f868d487559228f32
```

Datafiles: `external/how-restrictive-us-trade/data/imports-hs10/*data-current.parquet` (not provided)

## Dataset list

| Data file | Source | Notes | Provided |
| --- | --- | --- | --- |
| `data/2020_SML.csv` | OECD (2023) | Raw 2020 ICIO table | Yes |
| `data/ComponentsOfVa.xlsx` | BEA (n.d.) | Raw table plus author-added ICIO industry codes | Yes |
| `data/industry_names_elasts.csv` | Caliendo and Parro (2015); author | Industry names and trade elasticities | Yes |
| `data/concordances/Concordance_H6_to_H3.zip` | World Bank (n.d.) | Raw concordance | Yes |
| `data/concordances/hs2007_to_isic4_roy_zenodo.do` | Roy (2025) | Raw concordance (Stata syntax, read as text) | Yes |
| `programs/c/src/trump2025_tariffs.h` | Trade War Tracker (2026); concordances above | 2025 tariff increase by model sector and partner region, generated by `programs/python/trump2025_tariffs.py` | Yes |
| `external/how-restrictive-us-trade/` | Trade War Tracker (2026) | HS10 imports and duties | No |

## Computational requirements

### Software Requirements

- [x] The replication package contains a program that sets up the necessary directory structure (`run_replication.sh`; the makefile and scripts create their output directories). Python dependencies are listed in `requirements.txt`.

- Python 3.10.12
  - `numpy` 1.26.4
  - `pandas` 2.2.3
  - `matplotlib` 3.5.1
  - `seaborn` 0.13.2
  - `scikit-learn` 1.7.1
  - `scipy` 1.15.3
  - `openpyxl` 3.1.5
  - `pyarrow` 23.0.1 (only to rebuild the 2025 tariff schedule)
  - Run `pip install -r requirements.txt` once before running the code.
- C compiler and libraries (Ubuntu packages in parentheses):
  - GCC 11.4.0 with OpenMP (`build-essential`)
  - GNU Make 4.3
  - GSL 2.7.1 (`libgsl-dev`)
  - LAPACKE 3.10.0 (`liblapacke-dev`)
  - CBLAS from OpenBLAS 0.3.20 (`libopenblas-dev`)
- GNU bash 5.1

The scripts use bash and were run on Linux (Ubuntu 22.04). The Python scripts format numbers with the `en_US.utf8` locale, which must be installed, and request the Roboto font for figures (without it, Matplotlib falls back to its default font, which changes only the figures' appearance).

### Controlled Randomness

- [x] No pseudo random generator is used in the analysis described here. (The hierarchical clustering in `icio.py` is deterministic.)

### Memory, Runtime, Storage Requirements

#### Summary time to reproduce

Approximate time needed to reproduce the analyses on a 64-core workstation with 256 GB of RAM:

- [x] 8-24 hours

#### Summary of required storage space

Approximate storage space needed:

- [x] 250 MB - 2 GB (the package is about 70 MB; the model output is about 420 MB and the figures and tables about 25 MB)

#### Computational Details

The code was last run on a workstation with an AMD Ryzen Threadripper 3990X (64 cores, 128 threads), 251 GB of RAM, and Ubuntu 22.04.4 LTS (Linux 6.8).

- `icio.py` takes under a minute.
- Each of the 92 standard model solves (100-period horizon) needs about 23 GB of RAM and about 130 CPU-minutes; with 8 OpenMP threads and up to eight solves running at once, each took about 17 minutes (2 to 40 minutes), about 3.5 hours for all 92. Set `JOBS` to the number of solves the machine's memory can hold at once (for example, `JOBS=8` with 256 GB) and `OMP_NUM_THREADS` to the number of cores divided by `JOBS`.
- Each of the 18 solves with a 200-period horizon (`bin/model_nt200.bin`, Appendix H) needs about 90 GB of RAM; with 32 OpenMP threads and two solves running at once, each took about 45 minutes (40 to 100 minutes), about 7 hours for all 18. They run after the standard solves, `JOBS_NT200` at a time (default 1; `JOBS_NT200=2` fits in 256 GB), with `THREADS_NT200` threads each if set (otherwise `OMP_NUM_THREADS`).
- Creating the figures and tables takes about 5 minutes.

## Description of programs/code

- `run_replication.sh` runs every step below in order.
- `programs/python/icio.py` aggregates the ICIO table to the model's four regions and six sectors, clusters the goods industries into the four goods sectors (Section 3.1 and Appendix B), computes the sectors' shares of U.S. labor compensation, and writes the model's input-output matrix (`programs/python/output/iomat.txt`), the data figures (Figure 1), and the sector tables (Table 1, Table B.1). It prints the labor-compensation shares in the form hard-coded in `programs/c/src/calibrate.c` (`lshare`).
- `programs/python/trump2025_tariffs.py` measures the 2025 tariff increase on each HS10 product (December 2025 duties over imports minus the 2024 ratio, floored at zero), maps products to model sectors through the concordances, aggregates with 2024 import weights by partner region, and writes `programs/c/src/trump2025_tariffs.h`. It runs only with `run_replication.sh --tariffs`, since it needs the Trade War Tracker data.
- `programs/c/` contains the model. `src/main.c` reads the scenario options, `src/calibrate.c` calibrates the model to `iomat.txt` (writing `output/params.txt`), `src/eqm.c` defines the equilibrium conditions and writes the solutions, and `src/solver.c` and `src/gnewton.c` solve them. `src/globals.h` holds the model's dimensions and the alternative Armington elasticities and supply-chain frictions of Appendices D and E. `makefile` builds `bin/model.bin` (100-period horizon) and `bin/model_nt200.bin` (200 periods).
- `programs/c/run_all.sh` runs every model solve used in the paper. Each solve writes `output/vars0_*` (no-tariff equilibrium) and `output/vars1_*` (equilibrium with tariffs) by region, and a log. File names encode the scenario: `t25` 25% tariff; `s0`, `s1`, `s3`, `s4`, `s6` tariffs on "oil", "steel", "toys", "cars", all goods; `c2` all partners, `c0` China, `c3` low-tariff partners, `c4` high-tariff partners; `r1` retaliation, `r2` symmetric tariffs; `a4` all adjustment frictions, `a0` none, `a1`-`a3` one at a time; `trump2025` the 2025 schedule; suffixes `_ghh`, `_phiZ_lo/hi`, `_sig_NAME`, and `_psibX_nt200` for the appendix variants.
- `programs/c/calibrate_phiz.py` finds the supply-chain friction $\phi_Z$ that delivers a given short-run trade elasticity (Section 3, Appendices D and E). Its results are hard-coded in `programs/c/src/globals.h` (`PHIZ_LO`, `PHIZ_HI`, `SIG_PHIZ`) and `programs/python/model_io.py` (`PHIZ`), so it is not part of the main run. To repeat the calibration, run from `programs/c`: `python3 calibrate_phiz.py --target 1.0 --target 0.2` and, for each elasticity vector NAME, `python3 calibrate_phiz.py --sig NAME`.
- `programs/python/make_results.sh` runs `results.py` (figures) and `tables.py` (tables, and the CSV files behind the numbers in the text) on the model output, for the baseline or, with an option, for one of the appendix variants. `model_io.py` holds the code they share.
- `programs/python/sig_table.py`, `phiz_table.py`, and `psib_table.py` make the tables that compare the appendix variants, and `static_benchmark.py` solves the static model of Appendix I.

All generated files are written to `programs/python/output/` (figures, tables, CSV files) and `programs/c/output/` (model solutions and logs).

### License for Code

The code is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. See `LICENSE.txt`.

## Instructions to Replicators

- Install the software listed above, then run `pip install -r requirements.txt`.
- (Optional) To rebuild the 2025 tariff schedule, download the Trade War Tracker data as described above.
- From the root of the package, run `JOBS=8 ./run_replication.sh` (choose `JOBS` as described above), or `JOBS=8 ./run_replication.sh --tariffs` to also rebuild the 2025 tariff schedule.

### Details on various programs

- If running programs individually, ORDER IS IMPORTANT: `icio.py` (and `trump2025_tariffs.py`) before compiling the model; all model solves before `make_results.sh`; `psib_table.py` before `static_benchmark.py`. The Python scripts run from `programs/python`; `run_all.sh` and `calibrate_phiz.py` run from `programs/c`.
- The model solves in `run_all.sh` are independent of each other and can be split across machines, as long as all output files end up in `programs/c/output/`.
- `make_results.sh` stops before writing anything if a model solve it needs is missing.

## List of tables and programs

The provided code reproduces:

- [x] All numbers provided in text in the paper
- [x] All tables and figures in the paper

Output files are in `programs/python/output/`. Tables are written as LaTeX (`.tex`) and, for the generated tables, as `.csv` files, which also hold the numbers cited in the text.

| Figure/Table | Program | Output file | Note |
| --- | --- | --- | --- |
| Table 1 | `icio.py`; `programs/c/src/calibrate.c` | `sectors.tex`; `programs/c/output/params.txt` | Assembled by hand from these outputs and the parameter values in Section 3 |
| Figure 1(a)-(b) | `icio.py` | `fig_data_linkages_downstream_0.pdf`, `fig_data_linkages_upstream_0.pdf` | |
| Figure 1(c)-(f) | `icio.py` | `fig_data_sectoral_trade_by_use_{im,nx,im2,nx2}.pdf` | |
| Figure 2(a) | `results.py` | `fig_lr_emp_across_scenarios_all_vs_sectors.pdf` | |
| Figure 2(b)-(f) | `results.py` | `fig_dyn_labor_goods_s{2,1,4,0,3}.pdf` | In figure names, `s0`-`s4` are "oil", "steel", "toys", "cars", all goods |
| Figure 3(a) | `results.py` | `fig_lr_emp_vs_gdp.pdf` | |
| Figure 3(b) | `results.py` | `fig_dyn_tb_across.pdf` | |
| Figure 3(c) | `results.py` | `fig_lr_emp_across_scenarios_compare.pdf` | |
| Figure 3(d) | `results.py` | `fig_dyn_labor_goods_atb_noadj.pdf` | |
| Figure 3(e) | `results.py` | `fig_dyn_labor_goods_atb_retaliation.pdf` | |
| Figure 3(f) | `results.py` | `fig_dyn_labor_goods_trump2025.pdf` | |
| Table B.1 | `icio.py` | `industry_elasts.tex` | |
| Table B.2 | `icio.py`; `programs/c/src/globals.h` | `sectors.tex` | Assembled by hand |
| Table C.1 | `tables.py` | `tab_numbers.tex` | |
| Table C.2 | `tables.py` | `tab_numbers_sectors.tex` | |
| Table C.3 | `tables.py` | `tab_goods_cns.tex` | |
| Table D.1 | `programs/c/src/globals.h`; `sig_table.py` | `tab_sig.csv` | Assembled by hand from the elasticity vectors and the calibrated $\phi_Z$ |
| Table D.2 | `sig_table.py` | `tab_sig.tex` | |
| Table D.3 | `sig_table.py` | `tab_sig_dyn.tex` | |
| Figure D.1 | `results.py --sig NAME` | `fig_dyn_labor_goods_s4_sig_{cp,krs,half,u2,umean,u5}.pdf` | |
| Figure D.2 | `results.py --sig NAME` | `fig_dyn_labor_goods_trump2025_sig_{cp,krs,half,u2,umean,u5}.pdf` | |
| Figure D.3 | `results.py` | `fig_lr_emp_sig.pdf` | |
| Table E.1 | `phiz_table.py` | `tab_phiz.tex` | |
| Figure E.1 | `results.py`, `results.py --pz lo/hi` | `fig_dyn_labor_goods_{s4,trump2025}{_phiZ_lo,,_phiZ_hi}.pdf` | |
| Figure E.2 | `results.py --pz lo/hi` | `fig_lr_emp_across_scenarios_all_vs_sectors_phiZ_{lo,hi}.pdf` | |
| Table E.2 | `tables.py` | `tab_frictions.tex` | |
| Table F.1 | `tables.py` | `tab_blocs.tex` | |
| Table F.2 | `tables.py` | `tab_diversion.tex` | |
| Figure F.1 | `results.py` | `fig_div_imports_s{0,1,2,3,4}.pdf` | |
| Figure F.2 | `results.py` | `fig_dyn_labor_goods_s4{,_chn,_ltp,_htp}.pdf` | |
| Table G.1 | `tables.py` | `tab_ghh.tex` | |
| Figure G.1 | `results.py --ghh` | `fig_lr_emp_across_scenarios_all_vs_sectors_ghh.pdf` | |
| Figure G.2 | `results.py --ghh` | `fig_dyn_labor_decomp_{s4,trump2025}_ghh.pdf` | |
| Figure H.1 | `results.py` | `fig_dyn_tb_{s4,trump2025}.pdf` | |
| Table H.1 | `psib_table.py` | `tab_psib.tex` | |
| Figure H.2 | `results.py`, `results.py --psib 0.01 --nt 200` | `fig_dyn_nfa_across.pdf`, `fig_dyn_{nfa,tb}_across_psib0.01_nt200.pdf` | |
| Figure H.3 | `results.py --psib X --nt 200` | `fig_dyn_labor_goods_{s4,trump2025}_psib{0.005,0.01,0.02}_nt200.pdf` | |
| Table I.1 | `static_benchmark.py` | `tab_static.tex` | |

Numbers in the text come from these tables and from the CSV files written alongside them, mainly `tab_headline.csv` (goods employment, aggregates, and the trade balance by scenario), `tab_numbers*.csv`, `paths_external.csv` (the trade balance and its composition, net foreign assets, and real exports and imports by year), `paths_labor_decomp.csv`, `tab_diversion.csv`, `tab_frictions.csv`, `tab_ghh.csv`, `tab_psib.csv`, `tab_sig.csv`, `tab_phiz.csv`, `tab_static.csv`, and `tab_simple_benchmark.csv` (Appendix I). `icio.py` prints the sector shares of U.S. labor compensation and the clustering.

## References

Antràs, Pol, Davin Chor, Thibault Fally, and Russell Hillberry. 2012. "Measuring the Upstreamness of Production and Trade Flows." *American Economic Review* 102 (3): 412-16.

Bureau of Economic Analysis (BEA). n.d. "Components of Value Added by Industry [dataset]." GDP by Industry, sheet UVCT2-A. https://apps.bea.gov/iTable/ (accessed May 5, 2026).

Caliendo, Lorenzo, and Fernando Parro. 2015. "Estimates of the Trade and Welfare Effects of NAFTA." *Review of Economic Studies* 82 (1): 1-44.

OECD. 2023. "OECD Inter-Country Input-Output Tables, 2023 edition [dataset]." https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html (accessed May 5, 2026).

Roy, Jayjit. 2025. "Harmonized System (HS) 6 Digit to ISIC Rev 4 Concordance [dataset]." Zenodo. https://doi.org/10.5281/zenodo.16901880 (accessed May 12, 2026).

Trade War Tracker. 2026. "How Restrictive Is U.S. Trade Policy? [dataset and code]." GitHub repository, commit 6f2bf4af885ca2e0c5a37c6f868d487559228f32. https://github.com/tradewartracker/how-restrictive-us-trade.

World Bank. n.d. "WITS Product Concordance [dataset]." https://wits.worldbank.org/product_concordance.html (accessed May 12, 2026).
