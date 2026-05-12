# README for "Tariffs, Manufacturing Employment, and Supply Chains"

This repository contains the data, source code, generated intermediate files, and generated figures for:

Joseph B. Steinberg, "Tariffs, Manufacturing Employment, and Supply Chains," *NBER Working Paper 34236* (2025), https://www.nber.org/papers/w34236.

The code constructs a three-region, six-sector input-output representation of the world economy from public data, solves the dynamic general-equilibrium model described in Section 2 of the paper, and creates the figures used in Sections 3 and 4. The main replication path is:

1. Build the calibrated input-output data with `programs/python/icio.py`.
2. Compile and run the C model in `programs/c/`.
3. Create figures from the C model outputs with `programs/python/results.py`.
4. Optionally sync generated PDFs and generated table files to the external writing folder with `scripts/export_overleaf.sh`.

The code reproduces all manuscript figures and generated tables except the manually assembled calibration table in Table 1. The C code writes the internally calibrated parameter values to `programs/c/output/params.txt`; Table 1 summarizes those values and externally assigned parameters rather than listing every large matrix.

## Data Availability and Provenance Statements

All data used by this project are publicly available and are included in this repository under `data/`. The author has permission to redistribute the data contained in this replication package. See `LICENSE.txt` for license terms.

### Summary of Data Availability

| Data source | Repository file | Provided | Source and access |
| --- | --- | --- | --- |
| OECD 2023 Inter-Country Input-Output Tables, 2020 ICIO matrix | `data/2020_SML.csv` | Yes | Downloaded from the OECD Inter-Country Input-Output Tables page, https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html. Accessed May 5, 2026. |
| U.S. Bureau of Economic Analysis, Components of Value Added by Industry | `data/ComponentsOfVa.xlsx` | Yes | Downloaded from BEA GDP by Industry table TVA113, "Components of Value Added by Industry." The code uses sheet `UVCT2-A`, compensation of employees. The `ICIO Industry Code` column was manually added for merging with the ICIO industry codes. Accessed May 5, 2026. |
| Industry names and trade elasticities | `data/industry_names_elasts.csv` | Yes | Author-created file mapping ICIO industries to short names and trade elasticities. The `elast_CP` values are from the 99 percent sample estimates in Table 1 of Caliendo and Parro (2015). Accessed May 5, 2026. |
| Trade War Tracker HS10 import and duty data | `external/how-restrictive-us-trade/` | No | Cloned from https://github.com/tradewartracker/how-restrictive-us-trade. The Trump 2025 tariff script uses Census HS10 import values and calculated duties from this repository. |
| WITS HS2022-to-HS2007 concordance | `data/concordances/Concordance_H6_to_H3.zip` | Yes | World Bank WITS concordance file `Concordance_H6_to_H3.zip`, https://wits.worldbank.org/product_concordance.html. Used to bridge HS2022-style HS6 products to HS2007 before applying the ISIC concordance. Accessed May 12, 2026. |
| HS2007-to-ISIC Rev. 4 concordance | `data/concordances/hs2007_to_isic4_roy_zenodo.do` | Yes | Roy, "Harmonized System (HS) 6 Digit to ISIC Rev 4 Concordance," Zenodo, DOI `10.5281/zenodo.16901880`. Used to map HS2007 products to ISIC Rev. 4 two-digit industries. Accessed May 12, 2026. |
| WITS HS2022-to-GTAP concordance | `data/concordances/Concordance_H6_to_GP.zip` | Yes | World Bank WITS concordance file `Concordance_H6_to_GP.zip`, https://wits.worldbank.org/product_concordance.html. Preserved to reproduce the earlier GTAP-based tariff aggregation with `--concordance-method gtap`. Accessed May 12, 2026. |

### Dataset List

| Data file | Format | Used by | Notes |
| --- | --- | --- | --- |
| `data/2020_SML.csv` | CSV | `programs/python/icio.py` | Raw 2020 OECD ICIO matrix used to construct the three-region, six-sector input-output table. The paper discusses this data construction in Section 3.1, "Sectoral aggregation." |
| `data/ComponentsOfVa.xlsx` | Excel workbook | `programs/python/icio.py` | BEA compensation-of-employees data used to compute sectoral employment compensation shares for the sector classification summary. The paper references this source in the notes to Table 1. |
| `data/industry_names_elasts.csv` | CSV | `programs/python/icio.py` | Author-created industry-code crosswalk containing industry names and Caliendo-Parro trade elasticities. These elasticities enter the clustering described in Section 3.1 and summarized in Table 1. |
| `data/concordances/` | ZIP and Stata `.do` files | `programs/python/trump2025_tariffs.py` | Preserved public concordance files for aggregating HS tariff shocks into model sectors. See `data/concordances/README.md` for sources and SHA-256 hashes. |
| `programs/python/output/trump2025_tariffs.csv` | CSV | `programs/c/src/trump2025_tariffs.h` | Generated sector/source tariff rates for the Trump 2025 scenario. Rates are December 2025 HS10 duty rates minus 2024 baseline duty rates, clipped at zero and weighted by 2024 imports. |
| `programs/python/output/trump2025_tariff_diagnostics.csv` | CSV | Diagnostics | Generated diagnostics for the Trump 2025 tariff construction, including unmatched import shares and aggregate weighted tariff checks. |

## Computational Requirements

The code was last run on the following machine:

- OS: Ubuntu 22.04.4 LTS, Linux 6.8.0-94
- CPU: AMD Ryzen Threadripper 3990X, 64 cores / 128 threads
- RAM: 251 GiB
- Storage used by this repository: 177 MB total, including 66 MB in `data/`, 51 MB in `programs/c/output/`, and 1.6 MB in `programs/python/output/`

The full set of model runs is compute-intensive. On this workstation, `icio.py` takes less than 10 minutes, a single C model run takes about 151 seconds, `results.py` is much faster than the model solves, and `scripts/export_overleaf.sh` takes about 1 second. `programs/c/run_all.sh` runs 29 model exercises, so a complete clean rerun takes about 75 minutes if run serially.

### Software Requirements

The project was run without a Python virtual environment. The versions on the replication machine were:

- Python 3.10.12
- `numpy` 1.26.4
- `pandas` 2.2.3
- `matplotlib` 3.5.1
- `seaborn` 0.13.2
- `scikit-learn` 1.7.1
- `openpyxl` 3.1.5
- GCC 11.4.0
- GNU Make 4.3
- GNU bash 5.1.16
- GSL 2.7.1
- BLAS/LAPACK/LAPACKE 3.10.0
- OpenBLAS 0.3.20

On Ubuntu, the system-level requirements are provided by packages such as `build-essential`, `libgsl-dev`, `liblapacke-dev`, and `libopenblas-dev`. The Python scripts expect the `en_US.utf8` locale to be available. They also request the Roboto font for plot formatting; missing Roboto should mainly affect appearance, but missing locale support can cause the scripts to fail. If Matplotlib cannot write to its default cache directory, set a writable cache location before running the Python scripts:

```bash
export MPLCONFIGDIR=/tmp/matplotlib
```

Numerical results may differ slightly across machines because the C solver uses GSL, BLAS, LAPACK/LAPACKE, OpenMP, and compiler-dependent floating-point operations.

## Description of Programs and Outputs

### `programs/python/icio.py`

This script constructs the calibrated input-output data used by the C model and creates descriptive figures for the calibration section of the paper. It implements the data aggregation and sector clustering described in Section 3.1.

Inputs:

- `data/2020_SML.csv`
- `data/ComponentsOfVa.xlsx`
- `data/industry_names_elasts.csv`

Main outputs:

| Output | Description |
| --- | --- |
| `programs/python/output/iomat.txt` | Balanced three-region, six-sector input-output matrix, normalized by U.S. GDP. This is read by the C model in `programs/c/src/calibrate.c`. |
| `programs/python/output/sectors.tex` | Generated sector classification summary with industries, upstreamness, trade elasticity, and employment shares. The paper's Table 1 is manually assembled, but this file provides the generated sector-clustering information used in that table. |
| `programs/python/output/fig_data_upstream_vs_elast.pdf` | Auxiliary clustering plot showing industry upstreamness and trade elasticities. |
| `programs/python/output/fig_data_sectoral_trade_by_region_{ex,im,nx,ex2,im2,nx2}.pdf` | Auxiliary U.S. trade exposure figures by partner region. Suffixes without `2` are normalized by sectoral gross output; suffixes with `2` are normalized by U.S. GDP. |
| `programs/python/output/fig_data_sectoral_trade_by_use_{ex,im,nx,ex2,im2,nx2}.pdf` | U.S. trade exposure figures by use category: intermediate, investment, and consumption. These include several manuscript Figure 1 panels. |
| `programs/python/output/fig_data_linkages_downstream_{0,1,2,3,4}.pdf` | Downstream linkage figures. Index `0` uses all U.S. input sources, `1` domestic sources, `2` foreign sources, `3` China, and `4` rest of world. |
| `programs/python/output/fig_data_linkages_upstream_{0,1,2,3,4}.pdf` | Upstream linkage figures with the same index convention as downstream linkages. |

### `programs/python/trump2025_tariffs.py`

This script builds the Trump 2025 tariff matrix used by the C model. It expects Trade War Tracker data cloned to `external/how-restrictive-us-trade/` and writes sector/source tariff rates for China and rest of world.

The tariff shock is computed at the HS10 level as the December 2025 duty/import rate minus the 2024 duty/import baseline, clipped at zero, and then aggregated with 2024 import weights. HS10 products are truncated to HS6 and mapped to the model's four goods sectors through an ICIO-consistent concordance:

```text
HS10 tariff line -> HS6 product -> HS2007 HS6 -> ISIC Rev. 4 two-digit industry -> OECD ICIO industry -> model sector
```

The default concordance route uses the preserved files in `data/concordances/`: WITS HS2022-to-HS2007, Roy's HS2007-to-ISIC Rev. 4 rules, and the ICIO-to-model sector grouping generated by `programs/python/icio.py` and reported in `programs/python/output/sectors.tex`. The script first uses these local files; if they are missing, it downloads the public concordances. More detail, including the mapping tables, is in `docs/tariff_sector_concordance.md`.

The earlier WITS HS2022-to-GTAP route is preserved for comparison and can be run with:

```bash
python3 programs/python/trump2025_tariffs.py --concordance-method gtap
```

Main outputs:

| Output | Description |
| --- | --- |
| `programs/python/output/trump2025_tariffs.csv` | Sector/source tariff rates used by the model. |
| `programs/python/output/trump2025_tariff_diagnostics.csv` | Unmatched import shares and aggregate weighted-tariff checks. |
| `programs/c/src/trump2025_tariffs.h` | Generated C constants read by `set_tariffs()`. |

### `programs/c/`

The C code solves the dynamic general-equilibrium model described in Section 2 and the counterfactual tariff exercises described in Section 4.

Important files:

| File | Purpose |
| --- | --- |
| `programs/c/makefile` | Builds `programs/c/bin/model.bin`. |
| `programs/c/run_all.sh` | Runs exactly the C model exercises needed by `programs/python/results.py`. |
| `programs/c/src/main.c` | Parses model command-line options and runs the benchmark and tariff equilibrium solves. |
| `programs/c/src/calibrate.c` | Reads `programs/python/output/iomat.txt`, calibrates model parameters, and writes `programs/c/output/params.txt`. |
| `programs/c/src/trump2025_tariffs.h` | Generated Trump 2025 sector/source tariff matrix used when the model is run with `-p 1`. |
| `programs/c/src/eqm.c` | Defines equilibrium variables/equations, writes equilibrium time series, and writes solver seed files. |
| `programs/c/src/solver.c`, `programs/c/src/gnewton.c` | Nonlinear-equation solver routines. |

Main outputs:

| Output | Description |
| --- | --- |
| `programs/c/output/params.txt` | Internally calibrated scalar, vector, matrix, and 3D parameter values. This file underlies the internally calibrated parameter entries summarized in Table 1. |
| `programs/c/output/vars0_{usa,chn,row}_a*.csv` | Free-trade benchmark equilibrium time series by country/region and adjustment-cost setting. |
| `programs/c/output/vars1_{usa,chn,row}_t*_s*_c*_r*_d*_a*.csv` | Counterfactual tariff equilibrium time series by country/region and scenario. These are the inputs to `results.py`. |
| `programs/c/output/vars1_{usa,chn,row}_trump2025_r*_d*_a*.csv` | Trump 2025 tariff-matrix equilibrium time series. |
| `programs/c/output/log_*.txt` | Solver logs for the corresponding model runs, including command-line scenario information and convergence output. |
| `programs/c/output/seed0.bin`, `programs/c/output/seed1.bin` | Binary solver state files written after successful solves. `seed0.bin` is for the free-trade benchmark and `seed1.bin` is for the tariff equilibrium. If `read_seed` is enabled in `main.c`, `eqm.c` can read these files as initial guesses for the Newton solver in later counterfactuals. |

The model binary uses these command-line options:

| Option | Meaning | Values used in replication |
| --- | --- | --- |
| `-t` | Tariff rate in percent | `25` |
| `-s` | Targeted sector/use flag | `0` oil, `1` steel, `3` toys, `4` cars, `6` all goods; `10,11,13,14,16` intermediate-only analogues; `20,21,23,24,26` final-good-only analogues |
| `-c` | Targeted country/region | `0` China, `1` rest of world, `2` China + rest of world |
| `-r` | Retaliation flag | `0` no retaliation, `1` symmetric retaliation |
| `-d` | Duration flag | `0` permanent tariff, `1` temporary tariff revoked after four periods |
| `-a` | Adjustment-cost flag | `0` no adjustment costs, `4` labor + capital + supply-chain adjustment costs |
| `-p` | Tariff policy flag | Optional; `0` scalar command-line tariff, `1` Trump 2025 sector/source tariff matrix |

Examples:

```bash
# 25 percent permanent U.S. tariff on all goods from China and rest of world,
# with all adjustment costs.
./bin/model.bin -t 25 -s 6 -c 2 -r 0 -d 0 -a 4

# Same tariff, but with symmetric retaliation.
./bin/model.bin -t 25 -s 6 -c 2 -r 1 -d 0 -a 4

# 25 percent permanent U.S. tariff on final goods in the toys sector only.
./bin/model.bin -t 25 -s 23 -c 2 -r 0 -d 0 -a 4

# Trump 2025 sector/source tariff matrix. The scalar tariff, sector, and
# country flags are ignored for tariff assignment under -p 1.
./bin/model.bin -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1
```

Use `programs/c/run_all.sh` to comprehensively produce the C outputs used by the paper figures.

### `programs/python/results.py`

This script reads C model outputs and generates the model-result figures used in Section 4, plus auxiliary output used outside the manuscript.

Inputs:

- `programs/c/output/vars0_usa_a0.csv`
- `programs/c/output/vars0_usa_a4.csv`
- `programs/c/output/vars1_usa_t25_*.csv` for the scenarios in `programs/c/run_all.sh`

Main output groups:

| Output pattern | Description |
| --- | --- |
| `fig_dyn_labor_goods_s{0,1,2,3,4}.pdf` | Dynamic U.S. goods-sector employment effects. Scenario index `0` oil, `1` steel, `2` toys, `3` cars, `4` all goods. Several files are used in Figure 2. |
| `fig_dyn_labor_sectors_s{0,1,2,3,4}.pdf` | Dynamic employment effects for goods, services, and construction. |
| `fig_dyn_y_goods_s{0,1,2,3,4}.pdf` | Dynamic gross-output effects by goods sector. |
| `fig_dyn_inv_goods_s{0,1,2,3,4}.pdf` | Dynamic investment effects by goods sector. |
| `fig_dyn_rks_goods_s{0,1,2,3,4}.pdf` | Dynamic relative capital-cost effects by goods sector. |
| `fig_dyn_pm2_goods_s{0,1,2,3,4}.pdf` | Dynamic relative material-cost effects by goods sector. |
| `fig_dyn_pm_goods_s{0,1,2,3,4}.pdf` | Dynamic intermediate-price effects by goods sector. |
| `fig_dyn_pf_goods_s{0,1,2,3,4}.pdf` | Dynamic final-price effects by goods sector. |
| `fig_dyn_imports_goods_s{0,1,2,3,4}.pdf` | Dynamic import effects by goods sector. |
| `fig_dyn_exports_goods_s{0,1,2,3,4}.pdf` | Dynamic export effects by goods sector. |
| `fig_dyn_nx_goods_s{0,1,2,3,4}.pdf` | Dynamic net-export effects by goods sector. |
| `fig_dyn_rer_s{0,1,2,3,4}.pdf` | Dynamic real-exchange-rate effects. |
| `fig_dyn_macro.pdf` | Dynamic real GDP effects across scenarios. |
| `fig_dyn_consumption.pdf` | Dynamic consumption effects across scenarios. |
| `fig_dyn_labor_goods_atb_resize.pdf` | Resized all-goods employment-dynamics figure. |
| `fig_dyn_labor_goods_atb_noadj.pdf` | All-goods employment dynamics without adjustment costs, used in Figure 3. |
| `fig_dyn_labor_goods_atb_retaliation.pdf` | All-goods employment dynamics with symmetric retaliation, used in Figure 3. |
| `fig_dyn_labor_goods_atb_temp.pdf` | All-goods employment dynamics with temporary tariffs, used in Figure 3. |
| `fig_lr_emp_across_scenarios_all.pdf` | Long-run total goods employment effects. |
| `fig_lr_emp_across_scenarios_all_vs_sectors.pdf` | Long-run total and sectoral goods employment effects, used in Figure 2. |
| `fig_lr_emp_across_scenarios_chn.pdf` | Long-run goods employment effects for China-only tariffs. |
| `fig_lr_emp_across_scenarios_row.pdf` | Long-run goods employment effects for rest-of-world-only tariffs. |
| `fig_lr_emp_across_scenarios_compare.pdf` | Long-run goods employment effects for tariffs on one vs. both foreign regions, used in Figure 3. |
| `fig_lr_emp_vs_gdp.pdf` | Long-run goods employment vs. real GDP effects, used in Figure 3. |
| `fig_lr_emp_all_vs_m_vs_f.pdf` | Long-run employment effects of tariffs on intermediates, final goods, or both, used in Figure 3. |
| `fig_lr_goods_labor_s{0,1,2,3,4}.pdf` | Long-run sectoral goods employment effects by scenario. |

### `scripts/export_overleaf.sh`

This script syncs generated PDFs and generated `.tex` files from `programs/python/output/` to the folder where the paper is written. It is not needed to reproduce the outputs in this repository. It deletes old generated PDFs in the destination `figures/` directory and old generated `.tex` files in the destination `tables/` directory, but it preserves the manual `calibration.tex` table.

## Instructions to Replicators

Run commands from the repository root unless otherwise stated.

1. Ensure the Python packages and C libraries listed above are installed.

2. Optionally set a writable Matplotlib cache directory:

```bash
export MPLCONFIGDIR=/tmp/matplotlib
```

3. Construct the model input data and descriptive calibration figures:

```bash
cd programs/python
python3 icio.py
```

4. Optional: regenerate the Trump 2025 sector/source tariff matrix:

```bash
python3 trump2025_tariffs.py
```

The generated C header is already included in the repository, so this step is only needed if the tariff construction should be refreshed.

5. Compile the C model:

```bash
cd ../c
make model
```

6. Run all model exercises used by the paper figures:

```bash
./run_all.sh
```

7. Create figures from the model outputs:

```bash
cd ../python
python3 results.py
```

8. Optional: sync generated outputs to the external writing folder:

```bash
cd ../..
./scripts/export_overleaf.sh
```

The default destination for the export script is `../tariffs-mfg-emp-overleaf`. A different destination can be passed as the first argument.

## List of Tables/Figures and Programs

The public material in this repository reproduces all figures in the paper and all generated tables except the manually assembled calibration table.

| Paper output | Program | Output file(s) | Notes |
| --- | --- | --- | --- |
| Table 1, Calibration | Manual table in writing folder; inputs from `icio.py` and C calibration | `programs/python/output/sectors.tex`, `programs/c/output/params.txt` | The final table is manually assembled. Externally assigned parameter values come from the cited sources. Internally calibrated values are written by the C code to `params.txt`. |
| Figure 1(a), Intermediate inputs (% GO) | `programs/python/icio.py` | `programs/python/output/fig_data_linkages_downstream_0.pdf` | Described in Section 3.1. |
| Figure 1(b), Intermediate sales (% GO) | `programs/python/icio.py` | `programs/python/output/fig_data_linkages_upstream_0.pdf` | Described in Section 3.1. |
| Figure 1(c), Imports (% GO) | `programs/python/icio.py` | `programs/python/output/fig_data_sectoral_trade_by_use_im.pdf` | Described in Section 3.1. |
| Figure 1(d), Trade balance (% GO) | `programs/python/icio.py` | `programs/python/output/fig_data_sectoral_trade_by_use_nx.pdf` | Described in Section 3.1. |
| Figure 1(e), Imports (% GDP) | `programs/python/icio.py` | `programs/python/output/fig_data_sectoral_trade_by_use_im2.pdf` | Described in Section 3.1. |
| Figure 1(f), Trade balance (% GDP) | `programs/python/icio.py` | `programs/python/output/fig_data_sectoral_trade_by_use_nx2.pdf` | Described in Section 3.1. |
| Figure 2(a), Long-run goods employment | `programs/python/results.py` | `programs/python/output/fig_lr_emp_across_scenarios_all_vs_sectors.pdf` | Discussed in Section 4.1. |
| Figure 2(b), Employment dynamics: tariffs on Toys | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_s2.pdf` | Discussed in Section 4.2. |
| Figure 2(c), Employment dynamics: tariffs on all goods | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_s4.pdf` | Discussed in Section 4.2. |
| Figure 2(d), Employment dynamics: tariffs on Steel | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_s1.pdf` | Discussed in Section 4.2. |
| Figure 2(e), Employment dynamics: tariffs on Oil | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_s0.pdf` | Discussed in Section 4.2. |
| Figure 2(f), Employment dynamics: tariffs on Cars | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_s3.pdf` | Discussed in Section 4.2. |
| Figure 3(a), Macroeconomic consequences | `programs/python/results.py` | `programs/python/output/fig_lr_emp_vs_gdp.pdf` | Discussed in Section 4.3. |
| Figure 3(b), Tariffs on one vs. both countries | `programs/python/results.py` | `programs/python/output/fig_lr_emp_across_scenarios_compare.pdf` | Discussed in Section 4.4. |
| Figure 3(c), Tariffs on only intermediates or final goods | `programs/python/results.py` | `programs/python/output/fig_lr_emp_all_vs_m_vs_f.pdf` | Discussed in Section 4.5. |
| Figure 3(d), Tariffs on all goods without adjustment costs | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_atb_noadj.pdf` | Discussed in Section 4.6. |
| Figure 3(e), Tariffs on all goods with retaliation | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_atb_retaliation.pdf` | Discussed in Section 4.7. |
| Figure 3(f), Temporary tariffs | `programs/python/results.py` | `programs/python/output/fig_dyn_labor_goods_atb_temp.pdf` | Discussed in Section 4.8. |

## References

Antras, Pol, Davin Chor, Thibault Fally, and Russell Hillberry. 2012. "Measuring the Upstreamness of Production and Trade Flows." *American Economic Review* 102(3): 412-16.

Bureau of Economic Analysis. n.d. "Components of Value Added by Industry." GDP by Industry, table TVA113, sheet UVCT2-A. Accessed May 5, 2026. https://apps.bea.gov/iTable/

Caliendo, Lorenzo, and Fernando Parro. 2015. "Estimates of the Trade and Welfare Effects of NAFTA." *Review of Economic Studies* 82(1): 1-44. https://doi.org/10.1093/restud/rdu035

OECD. 2023. "OECD Inter-Country Input-Output Tables." Accessed May 5, 2026. https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html

Roy, Jayjit. 2025. "Harmonized System (HS) 6 Digit to ISIC Rev 4 Concordance." Zenodo. Accessed May 12, 2026. https://doi.org/10.5281/zenodo.16901880

Steinberg, Joseph B. 2025. "Tariffs, Manufacturing Employment, and Supply Chains." *NBER Working Paper* 34236. https://doi.org/10.3386/w34236

World Bank. n.d. "WITS Product Concordance." Accessed May 12, 2026. https://wits.worldbank.org/product_concordance.html
