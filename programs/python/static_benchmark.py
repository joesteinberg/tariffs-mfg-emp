"""Static steady-state benchmark for the U.S. (referee 2, major comment 1).

A small open economy version of the model's steady state, calibrated to the same U.S. block of the
ICIO table the C code uses (output/iomat.txt), with:
  - aggregate labor fixed, mobile across sectors;
  - the trade balance fixed at its no-tariff value (so tariffs cannot change net foreign assets);
  - foreign prices fixed at one and foreign demand for U.S. exports A_s * p_s^(-zeta_s);
  - the same technologies as the model: Leontief between value added and intermediates, Cobb-Douglas
    value added, CES (rho_c) consumption, Cobb-Douglas investment, Armington composites by use;
  - capital supplied elastically at the steady-state user cost, proportional to the price of
    investment goods; steady-state investment is a fixed fraction of capital income.

Usage:  python static_benchmark.py   (writes output/tab_static.{csv,tex} and prints a check)
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output')
NS, NC, CNS = 6, 4, 5
T = range(5)                                   # traded sectors (4 goods + services)
GOODS = range(4)
NAMES = ['Oil', 'Steel', 'Toys', 'Cars']
ZETA = np.array([8.0, 2.0, 4.0, 1.05, 1.05])   # baseline Armington elasticities (calibrate.c)
RHO_C = 0.65
LSHARE = np.array([0.029574, 0.014399, 0.018609, 0.041832, 0.843225, 0.052360])  # calibrate.c
CAP_SHARE = 0.34


def load_data():
    M = np.loadtxt(os.path.join(OUT, 'iomat.txt')).reshape(NS * NC + 2, NS * NC + 2 * NC + 1)
    d = {}
    d['go'] = M[NS * NC + 1, :NS].copy()
    d['va'] = M[NS * NC, :NS].copy()
    # intermediate purchases by US sector j of sector r, by source (0 = US)
    zm = np.array([[[M[src * NS + r, j] for j in range(NS)] for src in range(NC)] for r in T])
    d['a'] = zm.sum(1) / d['go']                  # a[r, j]: value of r per unit of j's output
    d['mu_m'] = zm[:, 0, :].sum(1) / zm.sum((1, 2))
    zc = np.array([[M[src * NS + r, NS * NC] for src in range(NC)] for r in T])
    zx = np.array([[M[src * NS + r, NS * NC + NC] for src in range(NC)] for r in T])
    d['mu_f'] = (zc[:, 0] + zx[:, 0]) / (zc + zx).sum(1)
    d['cons'] = zc.sum(1)                          # by sector, all sources
    d['inv'] = np.append(zx.sum(1), M[CNS, NS * NC + NC])   # incl. construction
    d['x'] = np.array([M[r, NS:NS * NC].sum() + M[r, NS * NC + 1:NS * NC + NC].sum()
                       + M[r, NS * NC + NC + 1:NS * NC + 2 * NC].sum() for r in T])
    d['imp'] = (zm[:, 1:, :].sum((1, 2)) + zc[:, 1:].sum(1) + zx[:, 1:].sum(1))
    d['nx'] = d['x'].sum() - d['imp'].sum()
    return d


def calibrate(d, zeta=ZETA, alpha=None, io=None):
    """Parameters. alpha: capital shares (default: from LSHARE as in calibrate.c).
    io: optional (lam_v, a) to override each sector's input structure."""
    p = dict(d)
    p['zeta'] = np.asarray(zeta, float)
    if alpha is None:
        alpha = 1.0 - 100 * (1 - CAP_SHARE) * LSHARE / d['va']
    p['alpha'] = np.asarray(alpha, float)
    p['lam_v'] = d['va'] / d['go']
    if io is not None:
        p['lam_v'], p['a'] = io
    p['eps_c'] = d['cons'] / d['cons'].sum()
    p['eps_x'] = d['inv'] / d['inv'].sum()
    rk = (p['alpha'] * d['va']).sum()
    p['s_inv'] = d['inv'].sum() / rk               # steady-state investment / capital income
    p['L'] = ((1 - p['alpha']) * d['va']).sum()    # labor endowment (wage = 1 in the data)
    return p


def prices(w, p, tm, tf):
    """Home output prices, composite prices, and the user cost of capital, given the wage."""
    z = p['zeta']
    pr = np.ones(NS)
    for _ in range(500):
        pm = (p['mu_m'] * pr[:5] ** (1 - z) + (1 - p['mu_m']) * (1 + tm) ** (1 - z)) ** (1 / (1 - z))
        pf = (p['mu_f'] * pr[:5] ** (1 - z) + (1 - p['mu_f']) * (1 + tf) ** (1 - z)) ** (1 / (1 - z))
        px = np.prod(np.append(pf, pr[CNS]) ** p['eps_x'])
        R = px                                     # user cost proportional to the investment price
        cva = w ** (1 - p['alpha']) * R ** p['alpha']
        new = p['lam_v'] * cva + pm @ p['a']
        if np.max(np.abs(new - pr)) < 1e-14:
            break
        pr = new
    pc = (p['eps_c'] @ pf ** (1 - RHO_C)) ** (1 / (1 - RHO_C))
    return pr, pm, pf, px, R, cva, pc


def quantities(Y, w, p, tm, tf, pr_all):
    """Residuals of goods-market clearing and the labor-market residual, given output Y."""
    pr, pm, pf, px, R, cva, pc = pr_all
    z = p['zeta']
    va = p['lam_v'] * cva * Y
    gdp = va.sum()
    rk = (p['alpha'] * va).sum()
    inv = p['s_inv'] * rk
    gm = p['a'] @ Y
    # tariff revenue depends on consumption, which depends on revenue: solve the linear fixed point
    sh_m_f = (1 - p['mu_m']) * ((1 + tm) / pm) ** (-z)        # import quantity per unit composite
    sh_f_f = (1 - p['mu_f']) * ((1 + tf) / pf) ** (-z)
    cshare = p['eps_c'] * (pf / pc) ** (1 - RHO_C)           # consumption spending shares
    tm_rev = (tm * sh_m_f * gm).sum()
    gf_inv = p['eps_x'][:5] * inv / pf
    # T = tm_rev + sum tf*sh_f_f*(cshare*cons/pf + gf_inv), cons = gdp + T - nx - inv
    k = (tf * sh_f_f * cshare / pf).sum()
    base = tm_rev + (tf * sh_f_f * gf_inv).sum()
    cons = (gdp - p['nx'] - inv + base) / (1 - k)
    gf = cshare * cons / pf + gf_inv
    hm = p['mu_m'] * (pr[:5] / pm) ** (-z) * gm
    hf = p['mu_f'] * (pr[:5] / pf) ** (-z) * gf
    x = p['x'] * pr[:5] ** (-z)
    res = np.empty(NS + 1)
    res[:5] = Y[:5] - (hm + hf + x)
    res[CNS] = Y[CNS] - p['eps_x'][CNS] * inv / pr[CNS]
    lab = (1 - p['alpha']) * va / w
    res[NS] = lab.sum() - p['L']
    return res, lab


def solve(p, tm=np.zeros(5), tf=np.zeros(5), guess=None):
    def f(v):
        w, Y = v[0], v[1:]
        return quantities(Y, w, p, tm, tf, prices(w, p, tm, tf))[0]
    v0 = np.append(1.0, p['go']) if guess is None else guess
    v, info, ier, msg = fsolve(f, v0, full_output=True, xtol=1e-13)
    if ier != 1 or np.max(np.abs(f(v))) > 1e-8:
        raise RuntimeError(msg)
    w, Y = v[0], v[1:]
    pa = prices(w, p, tm, tf)
    lab = quantities(Y, w, p, tm, tf, pa)[1]
    return dict(w=w, Y=Y, L=lab, v=v, pr=pa)


def goods_emp_change(p, s, tau=0.25, use='both'):
    tm, tf = np.zeros(5), np.zeros(5)
    sl = list(GOODS) if s == 'All' else [s]
    if use in ('both', 'm'):
        tm[sl] = tau
    if use in ('both', 'f'):
        tf[sl] = tau
    b = solve(p)
    c = solve(p, tm, tf, guess=b['v'])
    return 100 * (c['L'][:4].sum() / b['L'][:4].sum() - 1)


def simple_benchmark(d, alpha, zeta=ZETA, eps=1e-6):
    """The simple benchmark of the appendix's Proposition: no intermediates, labor share
    theta_j = 1 - alpha_j, capital elastic at a fixed rental with income to the household, Cobb-Douglas
    spending shares, U.S. home sales, exports, and imports (imports scaled so that trade is balanced).
    Returns the Proposition's terms for a small tariff on each goods sector."""
    from scipy.optimize import brentq
    th = 1 - np.asarray(alpha, float)[:5]
    x = d['x'].copy()
    m = d['imp'] * x.sum() / d['imp'].sum()
    dom = d['go'][:5] - x
    h = dom / (dom + m)
    beta = (dom + m) / (dom + m).sum()
    L = (th * d['go'][:5]).sum()

    def eqm(s, tau):
        t = np.zeros(5)
        if s is not None:
            t[s] = tau

        def f(w):
            pr = w ** th
            hs = h * pr ** (1 - zeta) / (h * pr ** (1 - zeta) + (1 - h) * (1 + t) ** (1 - zeta))
            xv = x * pr ** (1 - zeta)
            k = ((1 - th) * hs * beta).sum() + (t * (1 - hs) * beta / (1 + t)).sum()
            E = (w * L + ((1 - th) * xv).sum()) / (1 - k)
            return xv.sum() - ((1 - hs) * beta * E / (1 + t)).sum(), hs, E, xv
        w = brentq(lambda w: f(w)[0], 0.3, 3.0)
        _, hs, E, xv = f(w)
        lab = th * (hs * beta * E + xv) / w
        return lab[:4].sum(), E / w
    kt = th * (x + (zeta - 1) * th * (x + h * m))
    kap_g = kt[:4].sum() / kt.sum()
    c = th * h * beta
    C, c_g = c.sum(), c[:4].sum() / c.sum()
    g0, e0 = eqm(None, 0.0)
    out = {}
    for s, n in enumerate(NAMES):
        g1, e1 = eqm(s, eps)
        dE = (e1 - e0) / eps
        imp_sub = (1 - kap_g) * th[s] * (zeta[s] - 1) * h[s] * m[s]
        inc = -C * (kap_g - c_g) * dE
        out[n] = {'import substitution': imp_sub, 'income effect': inc, 'dE': dE,
                  'formula': imp_sub + inc, 'numerical': (g1 - g0) / eps,
                  '25% tariff': 100 * (eqm(s, 0.25)[0] / g0 - 1)}
    return pd.DataFrame(out).T, kap_g, c_g


if __name__ == '__main__':
    d = load_data()
    p0 = calibrate(d)
    b = solve(p0)
    print('baseline check: w = %.10f, max |Y - GO| = %.2e' % (b['w'], np.max(np.abs(b['Y'] - d['go']))))
    goods_alpha = (p0['alpha'][:4] * d['va'][:4]).sum() / d['va'][:4].sum()
    alpha_u = p0['alpha'].copy()
    alpha_u[:4] = goods_alpha
    zeta_u = ZETA.copy()
    zeta_u[:4] = 3.25
    # uniform goods input structure: each goods sector gets the goods-sector average VA share and
    # input coefficients (output-weighted)
    wgo = d['go'][:4] / d['go'][:4].sum()
    lam_u, a_u = p0['lam_v'].copy(), p0['a'].copy()
    lam_u[:4] = (p0['lam_v'][:4] * wgo).sum()
    a_u[:, :4] = (p0['a'][:, :4] @ wgo)[:, None]
    both = calibrate(d, zeta=zeta_u, alpha=alpha_u)
    specs = [('Static model', calibrate(d)),
             ('Uniform elasticities', calibrate(d, zeta=zeta_u)),
             ('Uniform capital shares', calibrate(d, alpha=alpha_u)),
             ('Uniform elasticities and capital shares', both),
             ('Uniform input structure', calibrate(d, io=(lam_u, a_u)))]
    pols = list(zip(list(GOODS) + ['All'], NAMES + ['All goods']))
    rows = {}
    # full model: steady state with net foreign assets fixed at their initial level (Appendix H)
    ps = pd.read_csv(os.path.join(OUT, 'tab_psib.csv'))
    ps = ps[ps['model'] == 'Steady state'].set_index('policy')['emp_lr']
    rows['Full model'] = {n: ps[n] for _, n in pols}
    for lab, p in specs:
        rows[lab] = {n: goods_emp_change(p, s) for s, n in pols}
    tab = pd.DataFrame(rows).T
    print(tab.round(2))
    for use, lab in [('m', 'intermediates only'), ('f', 'final goods only')]:
        print('Static model, %s:' % lab, {n: round(goods_emp_change(p0, s, use=use), 2) for s, n in pols})
    print('capital shares', dict(zip(NAMES + ['Svcs', 'Cons'], p0['alpha'].round(3))),
          'goods avg %.3f' % goods_alpha)
    sb, kap_g, c_g = simple_benchmark(d, p0['alpha'])
    print('Simple benchmark: kappa_G %.3f  c_G %.3f' % (kap_g, c_g))
    print(sb.round(3))
    sb.to_csv(os.path.join(OUT, 'tab_simple_benchmark.csv'))
    print('Simple benchmark, labor only:')
    print(simple_benchmark(d, np.zeros(NS))[0].round(3))
    tab.to_csv(os.path.join(OUT, 'tab_static.csv'))
    with open(os.path.join(OUT, 'tab_static.tex'), 'w') as f:
        f.write('\\begin{tabular}{l' + 'r' * len(pols) + '}\n\\toprule\n')
        f.write(' & ' + ' & '.join("``%s''" % n if n in NAMES else n for _, n in pols) + ' \\\\\n\\midrule\n')
        for i, (lab, r) in enumerate(tab.iterrows()):
            if i == 1:
                f.write('\\midrule\n')
            f.write(lab + ' & ' + ' & '.join('%.2f' % v for v in r.values) + ' \\\\\n')
        f.write('\\bottomrule\n\\end{tabular}\n')
