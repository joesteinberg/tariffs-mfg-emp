# Notation crosswalk: paper Appendix A ↔ C code

Appendix A ("Complete model", `app:model` in the Overleaf `paper.tex`) is written from
`programs/c/src/eqm.c`, `eqm.h`, `calibrate.c`, `calibrate.h`. Indices: country `i` (0 = USA,
1 = CHN, 2 = LTP, 3 = HTP), sector `s` (0–3 goods, `SVC` = 4, `CNS` = 5), source `j`, time `t`.
Traded set 𝒯 = sectors 0–4.

| Paper | Code | Notes |
|---|---|---|
| V^s_{i,t} | `va_t` | `a_ts`·`prod_va(k,l,A,alpha)`; `a_ts` = 1 |
| Ỹ, Y^s_{i,t} | `y_t` | `va/lam_va` minus labor adj. cost (Eq. A.3) |
| A^s_i, α^s_i | `A`, `alpha` | α = 1 − l0/va0 for USA; 0.34 elsewhere |
| λ^{s,v}_i, λ^{s,r}_i | `lam_va`, `lam` | |
| M^{s,r}_{i,t} | `md_t[i][s][r]` | |
| L^s_{i,t}, L̄_i | `l_t`, `ll0` | `fixl` = 1 baseline |
| K^s_{i,t}, X^s_{i,t} | `k_t`, `is_t` | |
| Φ(x), Φ'(x) | `phiK`, `dphiK` | φ_K = `etaK` = 0.52 |
| φ_L | `etaL` = 6.5 | cost in own output, × L_{t−1} |
| τ^k_i | `tauk` | US ≈ 0.073 |
| P^{s,y}, P^{s,m}, P^{s,f} | `py_t`, `pm_t`, `p_t` | |
| P^{s,v} | `py − DOT_PROD(lam, pm)` / `lam_va` | inside `mpk_rk`, `mpl_w` |
| P^x, X_{i,t} | `pi_t`, `ii_t` | `A^x` = `G`, ε^{s,x} = `eps[i][1][s]` |
| P^c, C, ε^{s,c}, ρ^c | `dcagg`-based, `cagg`, `eps[i][0][s]`, `rho` | code `rho` = (ρ^c−1)/ρ^c |
| numeraire (A.26) | `cpi_t[t][0]` = 1 (`price_norm`) | Laspeyres, weights `c0` |
| W_{i,t} | `w_t` | |
| Z^{s,m}_{i,j,t}, Z^{s,f}_{i,j,t} | `m2_t`, `q2_t` | |
| Ĝ, G^{s,m}, G^{s,f} | `prod_m`/`prod_q`, `m_t`, `q_t` | |
| μ^{s,m}, μ^{s,f} | `mu`·`M`^ρ, `theta`·`H`^ρ | scale absorbed into μ in the paper |
| ζ^{s,m}, ζ^{s,f} | `zeta`, `sig` | code stores (ζ−1)/ζ |
| φ_Z | `etaM` = `etaF` = `phiz` | |
| τ^{s,m}, τ^{s,f} | `tau_m_ts`, `tau_f_ts` | |
| Q_t (bond price) | `pb_t` | |
| Λ_{t+1}/Λ_t | `Q_t` | equals `pb_t` when ψ_B = 0 |
| B_{i,t} | `b_t` | last country by bond clearing |
| A.8 (labor FOC) | `mpl_w` | |
| A.9–A.10 (capital) | `mpk_rk` | evaluated at t+1 |
| A.13 (sourcing FOC) | `foc_m2`, `foc_q2` | |
| A.11 (composite) | `prod_m_chk`, `prod_q_chk` | |
| A.20 (Euler) | `euler` | |
| A.22 (goods clearing) | `mkt_clear_y` | |
| A.27 (current account) | `bop` | imposed for i = 0..2 |
| intratemporal consumption | `mucs_mucr` | goods vs. services |
| T = 100 | `NT` | |
