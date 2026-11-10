# Three-region results (archived)

Model output and generated figures/tables from the original three-region version of the model
(United States, China, rest of world), as they stood at commit `5e69ab5` on branch `jie-revision`.
The paper now uses the four-region model (United States, China, low-tariff and high-tariff
partners); its output is in `programs/c/output` and `programs/python/output`.

- `c-output/`: the C model's output (`vars*.csv`, `ss1_*.csv`, logs, and the phi_Z calibration
  caches `phiz_calibration_ae*.csv`). The `*_pz*.csv` files are intermediate solves from the phi_Z
  search (`calibrate_phiz.py`) and are not under version control.
- `python-output/`: figures and tables generated from that output by `results.py`, `tables.py`,
  `phiz_table.py`, `sig_table.py`, and `psib_table.py`.

To reproduce: `make model_nc3` builds the three-region model (`bin/model_nc3.bin`, which reads
`programs/python/output/iomat_nc3.txt` and writes files with the suffix `_nc3`); the Python scripts
read them with `--nc 3`. The phi_Z constants in `src/globals.h` (`PHIZ_LO`, `PHIZ_HI`,
`SIG_PHIZ`) are now the four-region values; the archived `_phiZ_*` and `_sig_*` runs used
`PHIZ_LO 0.6585`, `PHIZ_HI 12.44`, and `SIG_PHIZ = {2.9, 4.74, 2.966, 3.465, 2.488, 3.71, 4.745}`,
so to reproduce those, pass the value directly (`-z 0.6585`, or `-e NAME -z VALUE`).
(The archived files themselves have no `_nc3` suffix, because the three-region model was the
baseline when they were written.) The data inputs `iomat_nc3.{txt,tex}` and `trump2025_tariffs_nc3.csv`
remain in `programs/python/output`.
