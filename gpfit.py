"""Richardson fits of the convergence ladders and the requirements they define.

Every routine here takes arrays or a per-sweep frame and returns numbers; none reads a
file, holds state or draws anything. GP_CALIBRATION.ipynb supplies the data and the plots.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar, lsq_linear

# ── Calibrate s from L vs n_m fits  (Richardson, order k chosen by AIC) ──────
# Theory:  n_m^req = C_m · D_y^(s-q) · n_x·n_y·n_z           (eq. nm_req)
#
# Fit form (leading order p, truncation order k):
#   L(n_m) = L_inf − C1/(n_m/N₀)^p − … − Ck/(n_m/N₀)^(p·k)
# where N₀ = 1e6 is a scale used purely to keep the design-matrix columns
# O(1) — without it the matrix at 5–8 points with n_m^3 entries spans 1e10
# in column norms and the 4-parameter LS solve becomes ill-conditioned for
# the already-well-converged emittances (8 / 12 / 20 nm).
#
# p AND k. `p` is the LEADING convergence order — the physical rate at which
# the discretisation/sampling error dies. `k` is how many terms of the
# sub-leading series p, 2p, 3p… are kept; it does not change the rate, only how
# far down the pre-asymptotic tail the fit can reach. The two are partly
# degenerate — (p=1, k=2) and (p=2, k=1) both contain a 1/n² term — which is why
# k cannot be selected without first fixing p, and why p must be measured at
# k=1 (its defining limit) rather than at the selected k. See `profile_p`.
#
# p IS DERIVED FROM THE SCAN DESIGN, not asserted per sweep. Two mechanisms
# carry the error, with exponents fixed by physics:
#   P_SHOT = 1 — macroparticle shot noise. The controlling quantity is
#          macroparticles per CELL, m = n_m/(n_x n_y n_z) — which is exactly why
#          the theory reads n_m^req ∝ n_x n_y n_z — and the bias goes as 1/m.
#   P_CIC  = 2 — CIC deposition and bilinear field interpolation are both
#          second-order in the cell size h ∝ 1/n, so the term is O(h²). NGP
#          deposition would give 1 here.
# Which of the two LEADS depends on how m moves along the particular sweep, and
# that is a property of the scan design alone — no luminosity involved. Measure
#          alpha = dln(m)/dln(n)
# over the sampled points; the shot-noise term then decays as n^(-alpha*P_SHOT)
# and the CIC term as n^(-P_CIC), so
#          p = min(alpha*P_SHOT, P_CIC),   and p = P_CIC when alpha <= 0.
# alpha <= 0 means the shot-noise term does not decay with n: at alpha = 0 it is
# constant and is absorbed into L_inf, and at alpha < 0 it GROWS with n and
# cannot be a Richardson term at all. Either way only the CIC term converges.
# Measured on this data set (`derive_p` prints it):
#   n_m sweeps  alpha = +1.000 exactly (m ∝ n_m at fixed grid)        -> p = 1
#   n_y locus   alpha ~ 0 (-0.001/-0.006/-0.044 at eps_y = 1/2/4 nm)  -> p = 2
#               the locus is built to hold m fixed, and it does
#   n_x, n_z    alpha = -1 (m ∝ 1/n at fixed n_m)                     -> p = 2
# So n_m and the spatial sweeps genuinely need DIFFERENT p, and the reason is
# structural rather than a matter of taste.
# P_POLICY decides whether those numbers are ASSUMED or MEASURED. Either way
# `profile_p` runs and reports the measured p against the prior, so the prior is
# always audited. DEFAULT IS 'physics', because the measurement fails on this
# data set: see `profile_p` — n_m's chi2 is monotone in p (no interior minimum
# at all) and the k=1 model is mis-specified for n_x and n_y. Only n_z yields a
# clean measurement, p = 2.4 +- 0.1 against the CIC prior of 2. Under
# 'profiled' an unidentifiable profile falls back to the prior anyway.
#
# p IS THE DOMINANT SYSTEMATIC and it is not data-determined: at k=1 the
# n_m^req of the eps_y = 1 nm sweep moves 1.9e7 -> 2.2e6 as p goes 0.5 -> 2.
# It is therefore varied explicitly at the end of this cell rather than being
# buried inside a fitted error bar.
#
# n_m^req: last crossing of |L(n) − L_inf| = L_inf·tol  (brentq).
# Uncertainty: σ(log n_m^req) by MC marginalised over k, p and β.
# s regression: weighted log-log linear fit using 1/σ_log² as weights.


TOLERANCE_PCT = 5.0              # +-5 % convergence band
TOL_NM = TOLERANCE_PCT / 100.0   # 0.05
K_MAX  = 3                        # largest Richardson order offered to AIC
N_SCALE = 1e6                     # rescale n_m → n_m/N₀ for numerical stability

# Leading-order policy — see the header. 'profiled' measures one shared p per
# sweep family; 'physics' pins p to the prior below and reports the profiled
# value alongside as a check.
P_POLICY = 'derived'      # 'derived' = from the scan design (see header)
                          # 'profiled' = adopt the k=1 chi2 profile where it is
                          #              identifiable (it is not, on this data)
P_SHOT, P_CIC = 1.0, 2.0  # the only two physics inputs
ALPHA_TOL     = 1e-2      # |alpha| below this counts as "m held fixed"
P_BASE        = P_SHOT    # starting point for the n_m profile
P_DERIVED     = {}        # filled per sweep family by derive_p()

L_FRAC_NM = 0.75



def fit_Linf_fixed_p(g, n_col, p):
    """Weighted LS fit of L = L_inf - C / n^p with p fixed.
    Linear in (L_inf, C); closed-form solution weighted by 1/SEM^2.
    Returns L_inf, sigma_L_inf, C, sigma_C, chi2_red, n_points.
    """
    n   = g[n_col].values.astype(float)
    L   = g['mean_L'].values.astype(float)
    sem = g['sem_L'].values.astype(float)
    X = np.column_stack([np.ones_like(n), -1.0 / n**p])
    W = 1.0 / sem**2
    XtWX = (X * W[:, None]).T @ X
    XtWy = (X * W[:, None]).T @ L
    beta = np.linalg.solve(XtWX, XtWy)
    cov  = np.linalg.inv(XtWX)
    L_inf, C = beta[0], beta[1]
    resid = (L - (L_inf - C / n**p)) / sem
    dof   = len(L) - 2
    chi2_red = (resid**2).sum() / dof if dof > 0 else np.nan
    return dict(L_inf=L_inf, sigma_L_inf=np.sqrt(cov[0,0]),
                C=C, sigma_C=np.sqrt(cov[1,1]),
                chi2_red=chi2_red, n_points=len(L), p=p)


def scan_alpha(g, n_col, anchor):
    """alpha = dln(m)/dln(n) along a sweep, with m = n_m/(n_x n_y n_z).

    A property of the SCAN DESIGN only: it uses the grid and n_m actually run at
    each point and never touches the luminosity. `anchor` supplies the held
    values; any column present in `g` overrides it (the n_y locus carries its own
    per-point n_m)."""
    cols = {}
    for k in ('n_x', 'n_y', 'n_z', 'n_m'):
        cols[k] = (g[k].values.astype(float) if k in g.columns
                   else np.full(len(g), float(anchor[k])))
    cols[n_col] = g[n_col].values.astype(float)
    n = cols[n_col]
    m = cols['n_m'] / (cols['n_x'] * cols['n_y'] * cols['n_z'])
    if len(np.unique(n)) < 2 or not np.all(np.isfinite(m)) or np.any(m <= 0):
        return float('nan')
    return float(np.polyfit(np.log(n), np.log(m), 1)[0])


def effective_p(alpha):
    """Leading order implied by alpha. Derivation in the cell header."""
    if not np.isfinite(alpha) or alpha <= ALPHA_TOL:
        return P_CIC
    return float(min(alpha * P_SHOT, P_CIC))


def derive_p(datasets, n_col, anchors, label=''):
    """Derive the shared leading order for one sweep family from its design.

    datasets {key: g}, anchors {key: dict with n_x/n_y/n_z/n_m}. If the per-eps_y
    answers disagree the SMALLEST p is adopted (the slowest-decaying mechanism
    leads) and it is flagged, because a disagreement means the sweeps in the
    family were not all built the same way."""
    alphas = {e: scan_alpha(datasets[e], n_col, anchors[e]) for e in datasets}
    p_each = {e: effective_p(a) for e, a in alphas.items()}
    vals   = [v for v in p_each.values() if np.isfinite(v)]
    p      = float(min(vals)) if vals else P_CIC
    consistent = len(set(round(v, 6) for v in vals)) <= 1
    P_DERIVED[n_col] = p
    print(f'  p({label or n_col}) = {p:.2f}   from the scan design:  alpha = '
          + ', '.join(f'{e:g}:{alphas[e]:+.3f}' for e in alphas))
    if not consistent:
        print('      NOTE: per-eps_y p differ ('
              + ', '.join(f'{e:g}:{p_each[e]:.2f}' for e in p_each)
              + ') - the slowest-decaying is adopted')
    off = [e for e, a in alphas.items()
           if np.isfinite(a) and abs(a) > 0.1 and n_col == 'n_y']
    if off:
        print('      NOTE: alpha far from 0 at eps_y = '
              + ', '.join(f'{e:g}' for e in off)
              + ' - those sweeps are off the n_m locus (NM_FLOOR / wide'
                ' acceptance window), so m is not held fixed there')
    return dict(p=p, alphas=alphas, p_each=p_each, consistent=consistent)


def nm_asymptotic_points(g):
    """The n_m points inside the asymptotic range of their own sweep.

    One rule, one constant, no per-eps_y special cases. Returns `g` unchanged if
    it is empty."""
    if not len(g):
        return g
    return g[g['mean_L'] >= L_FRAC_NM * g['mean_L'].max()].reset_index(drop=True)


def _fit_kterm_core(g, n_col, p, k, n_scale, bounds=None):
    """One pass of weighted LS at k=k. bounds=None → unconstrained,
    bounds=tuple → constrained via lsq_linear."""
    n = g[n_col].values.astype(float) / n_scale
    L = g['mean_L'].values.astype(float)
    e = g['sem_L'].values.astype(float)
    N = len(L)
    if N < k + 1:
        return None
    X = np.column_stack([np.ones_like(n)] + [-1.0 / n**(p * i) for i in range(1, k + 1)])
    sw  = 1.0 / e
    Xw  = X * sw[:, None]
    Lw  = L * sw
    if bounds is None:
        try:
            beta = np.linalg.solve(Xw.T @ Xw, Xw.T @ Lw)
        except np.linalg.LinAlgError:
            return None
    else:
        try:
            beta = lsq_linear(Xw, Lw, bounds=bounds, method='trf', verbose=0).x
        except (ValueError, np.linalg.LinAlgError):
            return None
    resid = (L - X @ beta) / e
    dof   = N - (k + 1)
    chi2r = (resid**2).sum() / dof if dof > 0 else float('nan')
    try:
        cov = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        cov = np.eye(k + 1)
    return dict(beta=beta, cov=cov, chi2_r=chi2r, dof=dof, N=N, n_scale=n_scale)


def fit_kterm_local(g, n_col, p, k, n_scale=N_SCALE):
    """Hybrid k-term weighted LS:
       1. Try UNCONSTRAINED. If the resulting cubic has a positive crossing
          AND no Cᵢ is dramatically negative, keep it (preserves full k=3
          shape information when the data supports it).
       2. Otherwise fall back to CONSTRAINED with Cᵢ ≥ 0 (which always
          gives a defined crossing but may peg some Cᵢ to 0).
       Returns the chosen fit dict with an extra key 'fit_mode'."""
    fit_unc = _fit_kterm_core(g, n_col, p, k, n_scale, bounds=None)
    if fit_unc is not None:
        nm_unc = n_m_req_kterm(fit_unc['beta'], TOL_NM, p, k, n_scale=n_scale)
        # accept unconstrained if cubic has a positive root within the bracket
        if np.isfinite(nm_unc) and nm_unc > 0:
            fit_unc['fit_mode'] = 'unconstrained'
            return fit_unc

    lb = np.array([-np.inf] + [0.0] * k)
    ub = np.array([ np.inf] * (k + 1))
    fit_con = _fit_kterm_core(g, n_col, p, k, n_scale, bounds=(lb, ub))
    if fit_con is not None:
        fit_con['fit_mode'] = 'constrained (Cᵢ≥0)'
    return fit_con


def n_m_req_kterm(beta, tol, p, k, n_scale=N_SCALE):
    """Smallest positive n_m (original units) where Σᵢ Cᵢ/(n/n_scale)^(p·i) = L_inf·tol."""
    L_inf = beta[0]
    Cs    = beta[1:k + 1]
    target = L_inf * tol
    def f(n_norm):
        return sum(Cs[i - 1] / n_norm**(p * i) for i in range(1, k + 1)) - target
    try:
        n_norm_root = brentq(f, 1e-5, 1e5)
        return n_norm_root * n_scale
    except (ValueError, RuntimeError):
        return float('nan')


def n_req_kterm_signed(beta, tol, p, k, n_scale):
    """`n_m_req_kterm` rewritten on |L(n) − L_inf|.

    `n_m_req_kterm` solves Σᵢ Cᵢ/(n/n₀)^(p·i) = L_inf·tol, which assumes L rises
    to L_inf from below (all Cᵢ > 0) — true for every n_m sweep, so this returns
    identical answers there. Some n_x / n_z sweeps instead converge from ABOVE,
    giving Cᵢ < 0 and a target the original equation can never reach. Solving on
    the absolute deviation keeps the identical criterion — "the n beyond which
    the fit stays inside the ±tol band" — for both approach directions.

    Takes the LAST crossing so a non-monotone series cannot report convergence
    at an n where it later leaves the band again.
    """
    L_inf  = beta[0]
    Cs     = beta[1:k + 1]
    target = abs(L_inf) * tol

    def excess(n_norm):
        dev = sum(Cs[i - 1] / n_norm**(p * i) for i in range(1, k + 1))
        return abs(dev) - target

    # Vectorised bracket scan (this is called thousands of times by the MC),
    # then brentq refines the root exactly, so the grid only has to be fine
    # enough to isolate the crossing.
    grid = np.geomspace(1e-5, 1e5, 1201)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        vals = excess(grid)
    outside = np.where(np.isfinite(vals) & (vals > 0))[0]
    if len(outside) == 0:
        return float('nan')          # fit never leaves the band
    j = int(outside[-1])
    if j >= len(grid) - 1:
        return float('nan')          # still outside at the top of the bracket
    try:
        return float(brentq(excess, grid[j], grid[j + 1])) * n_scale
    except (ValueError, RuntimeError):
        return float('nan')


def fit_kterm_signed(g, n_col, p, k, n_scale):
    """`fit_kterm_local` with the same unconstrained→(Cᵢ≥0) hybrid, but judging
    the unconstrained fit with `n_req_kterm_signed`, so a converging-from-above
    fit is not rejected as "no positive crossing" and forced into the Cᵢ≥0
    branch (which can only bend upward, and degenerates to a flat line)."""
    fit_unc = _fit_kterm_core(g, n_col, p, k, n_scale, bounds=None)
    if fit_unc is not None:
        n_unc = n_req_kterm_signed(fit_unc['beta'], TOL_NM, p, k, n_scale)
        if np.isfinite(n_unc) and n_unc > 0:
            fit_unc['fit_mode'] = 'unconstrained'
            return fit_unc
    lb = np.array([-np.inf] + [0.0] * k)
    ub = np.array([np.inf] * (k + 1))
    fit_con = _fit_kterm_core(g, n_col, p, k, n_scale, bounds=(lb, ub))
    if fit_con is not None:
        fit_con['fit_mode'] = 'constrained (Cᵢ≥0)'
    return fit_con


def richardson_terms(beta, p, k, n, n_scale):
    """Contribution of each expansion term to L at abscissa `n`.

    T_i = C_i·(n/n₀)^(−p·i). Invariant under the choice of n₀ (rescaling n₀
    rescales C_i by exactly the compensating factor), so this is a property of
    the fitted curve, not of the parametrisation."""
    return np.array([beta[i] * (n / n_scale)**(-p * i) for i in range(1, k + 1)])


def is_asymptotic(beta, p, k, n_min, n_scale, n_max=None):
    """Is the fitted series actually an ASYMPTOTIC expansion over the range fitted?

    A Richardson expansion L_inf − ΣC_i·n^(−p·i) is an asymptotic series: it is
    only meaningful where its terms DECREASE, |T_1| > |T_2| > … > |T_k|. Checked
    at the smallest n in the fit, the worst case.

    WHY THIS MATTERS — the "snag". The i-th basis function varies by
    (n_max/n_min)^(p·i) across the fitted range. For the n_y sweeps (p=2, range
    32→512) the k=3 term varies by 1.7e7: it is numerically zero at every point
    except the lowest, where it can be tuned freely. AIC happily buys that χ²
    improvement, and the fitted curve hooks sharply at the lowest point while
    passing smoothly through the rest. Deleting that point does not help — the
    same term simply re-localises on the new lowest point.

    The hook is the visible symptom of terms that have stopped decreasing.
    Measured before this test was added, at ε_y = 1 nm, n_y = 32:
    T = (+6.39, −16.31, +11.39) against L_inf = 3.81 — three terms each larger
    than the luminosity itself, cancelling to pass through one point. That is a
    numerical artefact with no physical content, and any n^req extrapolated
    through it inherits the artefact.

    Enforcing this is not extra conservatism, it is the validity condition of the
    model already in use.

    NOT ALSO ENFORCED: strict monotonicity of |L(n) − L_inf| across the window.
    For k=2 that is exactly |T_2/T_1| <= 1/2 at n_min, and it was tested — it is
    too strong here. Three sweeps sit marginally above (0.51-0.58), and rejecting
    k=2 for them leaves k=1 with chi2 = 158/3, 220/3, 143/3, i.e. a plainly
    mis-specified fit, and the calibration degrades (q -> 0.00 +- 0.03,
    chi2 32.7/5). Decreasing terms is the defensible line; it bounds the residual
    curvature instead of forbidding it, putting the turnover at or just below the
    lowest fitted point rather than 2.5x above it.
    """
    _ = n_max          # accepted for signature stability; see note above
    if k < 2:
        return True
    T = np.abs(richardson_terms(beta, p, k, n_min, n_scale))
    return bool(np.all(np.isfinite(T)) and np.all(T[1:] < T[:-1]))


def select_k_aic(g, n_col, p, n_scale, k_max=K_MAX, require_asymptotic=True):
    """Choose the Richardson order k by AIC. Returns (best, table).

    Orders whose series is not asymptotic over the fitted range (`is_asymptotic`)
    are INADMISSIBLE — they are recorded in the table with `admissible=False` and
    excluded from both the selection and the Akaike weights used by the MC.

    The per-point errors are known (sample STD over the 15 canonical seeds), so
    the Gaussian log-likelihood gives, for K = k+1 free parameters and N points,

        AIC  = chi2 + 2K
        AICc = AIC + 2K(K+1)/(N − K − 1)        (reported, NOT used to select)

    Candidates are k = 1 … min(k_max, N−2): every order that still leaves at
    least one degree of freedom.

    WHY PLAIN AIC AND NOT AICc. AICc is the usual small-N recommendation and was
    tried here first, but at these N its correction is not a mild refinement —
    it is the dominant term. At N=5, K=3 it adds 24 to the score, so no
    attainable chi2 improvement can ever justify k=2. It therefore pinned k=1 on
    scans where k=1 is plainly wrong: the n_x sweep at eps_y = 20 nm gives
    chi2 = 18.6/3 at k=1 against 0.15/2 at k=2, and the n_z sweep at 0.5 nm gives
    37.1/2 against 0.04/1. AICc is derived for a correctly-specified model whose
    variance is estimated from the same data; here the variances are measured
    independently from 15 seeds and the k=1 model is often genuinely
    mis-specified, which is exactly where that penalty misbehaves. Both scores
    are printed so the choice stays auditable.

    Returns best=None when no order is admissible (N < 3).
    """
    N     = len(g)
    n_min = float(np.min(g[n_col].values.astype(float)))
    n_max = float(np.max(g[n_col].values.astype(float)))
    rows  = []
    for k in range(1, min(k_max, N - 2) + 1):
        fit = fit_kterm_signed(g, n_col, p, k, n_scale)
        if fit is None or not np.isfinite(fit['chi2_r']):
            continue
        K    = k + 1
        chi2 = float(fit['chi2_r'] * fit['dof'])
        aic  = chi2 + 2 * K
        adm  = (not require_asymptotic) or is_asymptotic(fit['beta'], p, k,
                                                         n_min, n_scale, n_max)
        rows.append(dict(k=k, K=K, N=N, dof=int(fit['dof']), chi2=chi2, aic=aic,
                         aicc=(aic + 2 * K * (K + 1) / (N - K - 1)
                               if N - K - 1 > 0 else float('inf')),
                         fit_mode=fit['fit_mode'], fit=fit, admissible=adm))
    live = [r for r in rows if r['admissible']]
    if not live:
        return None, rows
    best = min(live, key=lambda r: r['aic'])
    for r in rows:
        r['delta_aic'] = (r['aic'] - best['aic'] if r['admissible']
                          else float('inf'))
    return best, live


def fit_asymptotic_window(g, n_col, p, n_scale_fn=None, k_max=K_MAX, min_pts=3):
    """Fit on the LARGEST top-window over which the expansion stays asymptotic.

    The Richardson series is valid as n→∞. If no admissible order exists on the
    full sweep, the lowest point is outside the expansion's domain — the correct
    response is to drop it, NOT to add a term that spikes there. So: try all N
    points; if `select_k_aic` finds no admissible k, drop the lowest n and retry,
    down to `min_pts`.

    This is the standard asymptotic-range determination of a grid-convergence
    study, and it removes the snag by construction: with the terms required to
    decrease, the curve cannot hook at its own lower endpoint.

    Returns dict(window, n_scale, best, table, n_dropped, dropped) or None.
    """
    gs = g.sort_values(n_col).reset_index(drop=True)
    N  = len(gs)
    if n_scale_fn is None:
        n_scale_fn = n_scale_geom
    for m in range(N, min_pts - 1, -1):
        w  = gs.iloc[N - m:].reset_index(drop=True)
        ns = float(n_scale_fn(w[n_col].values))
        best, table = select_k_aic(w, n_col, p, ns, k_max=k_max)
        if best is not None:
            return dict(window=w, n_scale=ns, best=best, table=table,
                        n_dropped=N - m,
                        dropped=gs.iloc[:N - m][n_col].tolist())
    return None


def n_scale_geom(n):
    """Per-sweep rescale n → n/n₀ with n₀ = geometric mean of the sampled n.
    Conditioning only: no effect on L_inf, on n^req, or on any χ²."""
    return float(np.exp(np.mean(np.log(np.asarray(n, dtype=float)))))


def profile_p(datasets, n_col, scales=None, p_bounds=(0.3, 4.0), label=''):
    """Measure the shared leading order p for a family of sweeps (METHODS §3).

    `datasets` is {key: g} — one aggregated sweep per ε_y, same swept variable,
    same geometry. Deposition scheme and field solver do not change with
    emittance, so the leading order should not either; sharing p across ε_y
    costs one parameter for ~20-40 points.

    PROFILED AT k = 1, NOT AT THE AIC k. This is not a convenience — p and k are
    degenerate and the joint profile is UNBOUNDED. `p` is defined as the rate in
    L_inf − L(n) ~ C·n^(−p) as n→∞, which is the k=1 statement; the k≥2 terms
    sit at 2p, 3p and can absorb the leading behaviour themselves, so lowering p
    while raising k always buys χ². Measured here: at k=2 the n_x and n_y total
    χ² minima slide from p≈1.5 down to p≈0.75, and n_z goes flat (0.6-1.0 across
    p = 0.75-3.5). Profiling at the AIC k therefore runs p to the lower bound and
    produces n^req extrapolations that are pure model artefacts (n_m^req ~ 6e9).
    So: profile p at k=1, then select k by AIC at that p.

    IDENTIFIABILITY IS CHECKED, NOT ASSUMED. Two ways the measurement can fail,
    both reported:
      * the minimum runs to a bound → the data do not determine p at all. This
        is the case for n_m: its total χ² decreases monotonically all the way
        down (60.8 at p=0.3, 179.6 at p=1.0, 536.0 at p=2.0), i.e. the sweeps
        prefer an ever-slower tail because they are not deep enough into the
        asymptotic regime to measure a decay rate.
      * an interior minimum exists but the k=1 model is mis-specified there
        (χ²/ndf ≫ 1) → the "measured" p is describing curvature the model cannot
        represent, not a convergence rate.
    Either way `identifiable` comes back False and the physics prior stands.

    σ_p from Δχ²_total = 1 about the minimum — a 1σ interval only because χ² is
    built from `sem_L`, the error on the 15-seed MEAN. On an STD-scaled χ² the
    same Δχ²=1 would span ~√n_seeds ≈ 3.9 σ.
    """
    keys = list(datasets)
    if scales is None:
        scales = {e: n_scale_geom(datasets[e][n_col].values) for e in keys}

    def total_chi2(p):
        tot = 0.0
        for e in keys:
            f = fit_kterm_signed(datasets[e], n_col, p, 1, scales[e])
            if f is None or not np.isfinite(f['chi2_r']):
                return np.inf
            tot += f['chi2_r'] * f['dof']
        return tot

    r = minimize_scalar(total_chi2, bounds=p_bounds, method='bounded')
    p_hat, chi2_min = float(r.x), float(r.fun)

    def excess(p):
        return total_chi2(p) - chi2_min - 1.0
    try:
        p_lo = brentq(excess, p_bounds[0], p_hat)
    except (ValueError, RuntimeError):
        p_lo = p_bounds[0]
    try:
        p_hi = brentq(excess, p_hat, p_bounds[1])
    except (ValueError, RuntimeError):
        p_hi = p_bounds[1]
    sigma_p = max((p_hi - p_lo) / 2.0, 1e-3)

    n_pts   = sum(len(datasets[e]) for e in keys)
    ndf     = n_pts - 2 * len(keys) - 1          # k=1 → 2 params each, −1 for p
    chi2_r  = chi2_min / ndf if ndf > 0 else float('nan')
    at_edge = (min(abs(p_hat - p_bounds[0]), abs(p_hat - p_bounds[1])) < 0.02)
    well_specified = np.isfinite(chi2_r) and chi2_r < 3.0
    identifiable   = bool((not at_edge) and well_specified)

    prior = P_DERIVED.get(n_col, np.nan)
    pull  = (p_hat - prior) / sigma_p if np.isfinite(prior) else np.nan
    why   = ('' if identifiable else
             ('  <- NOT identifiable: profile at lower bound' if at_edge else
              f'  <- NOT identifiable: k=1 mis-specified (χ²_r={chi2_r:.1f})'))
    print(f'  p({label or n_col}) = {p_hat:.2f} ± {sigma_p:.2f}   '
          f'[prior {prior:g}, pull {pull:+.1f}σ]   '
          f'k=1 total χ²/ndf = {chi2_min:.1f}/{ndf}{why}')

    return dict(p=p_hat, sigma_p=sigma_p, chi2=chi2_min, ndf=ndf, chi2_r=chi2_r,
                n_pts=n_pts, scales=scales, prior=prior,
                identifiable=identifiable, at_edge=at_edge)


def p_in_use(prof, n_col):
    """The leading order actually adopted.

    'physics'  → the prior always (profile still reported as a check).
    'profiled' → the measured value, but ONLY where the profile is
                 identifiable; otherwise the prior, because an unidentified p is
                 not a measurement and letting it through silently replaces a
                 physical model with a numerical artefact.
    """
    if P_POLICY == 'profiled' and prof.get('identifiable'):
        return float(prof['p'])
    return float(prof['p_derived'])


def sigma_p_in_use(prof, n_col):
    """Width for the p marginalisation (METHODS §6). Zero when p is an input
    rather than a measurement — under 'physics', or where the profile failed."""
    if P_POLICY == 'profiled' and prof.get('identifiable'):
        return float(prof['sigma_p'])
    return 0.0        # a derived p is an input, not a measurement


def mc_sigma_log_n_req(g, n_col, p_hat, sigma_p, tol, n_scale,
                       n_samples=2000, rng_seed=42, k_max=K_MAX):
    """MC σ(log n^req), marginalised over k, p AND β  (METHODS §6).

    Per draw:  p ~ N(p̂, σ_p)  →  k ~ Akaike weights at p̂  →  refit β̂(p,k) by
    weighted LS (linear at fixed (p,k), so one small solve)  →  β ~ MVN(β̂,
    cov·max(1,χ²_r))  →  evaluate n^req.

    WHY EACH IS MARGINALISED
    k: not determined by the data. The Akaike weights run k=1: 0.5-0.7,
       k=2: 0.25-0.36, k=3: 0.09-0.13, so the runner-up carries a quarter to a
       half of the evidence; k=1 has 2 free parameters against k=3's 4, so
       pinning it both halves σ_log and moves the central value.
    p: fitted from the data (`profile_p`), not known. Holding it at p̂ inside
       the MC understates σ_log, and would leave σ(log n_y^req) omitting a
       source that σ(log n_m^req) also has — the two calibrations' χ² would not
       be comparable. Under P_POLICY='physics' p is an input, σ_p → 0, and this
       term correctly vanishes.

    Returns (sigma_log, median). Returns (nan, nan) on <50 surviving draws.
    """
    rng = np.random.default_rng(rng_seed)

    best, table = select_k_aic(g, n_col, p_hat, n_scale, k_max=k_max)
    if best is None:
        return float('nan'), float('nan')
    ks = np.array([r['k'] for r in table])
    w  = np.array([np.exp(-0.5 * r['delta_aic']) for r in table], float)
    w /= w.sum()

    sig_p = float(sigma_p) if (np.isfinite(sigma_p) and P_POLICY == 'profiled') else 0.0

    vals = []
    for _ in range(n_samples):
        p_d = rng.normal(p_hat, sig_p) if sig_p > 0 else p_hat
        if not np.isfinite(p_d) or p_d < 0.05:
            continue
        k_d = int(rng.choice(ks, p=w))
        f   = fit_kterm_signed(g, n_col, p_d, k_d, n_scale)
        if f is None:
            continue
        sc = max(1.0, f['chi2_r']) if np.isfinite(f['chi2_r']) else 1.0
        try:
            b = rng.multivariate_normal(f['beta'], f['cov'] * sc)
        except (ValueError, np.linalg.LinAlgError):
            continue
        v = n_req_kterm_signed(b, tol, p_d, k_d, n_scale)
        if np.isfinite(v) and v > 0:
            vals.append(v)
    if len(vals) < 50:
        return float('nan'), float('nan')
    a = np.array(vals)
    return float(np.log(a).std()), float(np.median(a))


def sigma_window(g, n_col, p, tol, n_scale_fn, k_max=K_MAX, min_pts=3):
    """σ_win — the SELECTION term (METHODS §6.2).

    `fit_asymptotic_window` returns one answer: the largest window over which the
    expansion stays asymptotic. But where the asymptotic range starts is a
    judgement the data only weakly constrain, and n^req is extrapolated FROM that
    window. So refit on every top-window down to `min_pts` — drop the lowest n,
    refit; drop the lowest two, refit; … — and take the spread of ln n^req across
    them. That spread is how much the answer depends on where the range was cut.

    Not hypothetical: with no lower cut at all the ε_y = 0.5 nm n_m sweep returns
    an n_m^req an order of magnitude below the adopted one.

    Returns 0.0 if only one window is usable, nan if none is."""
    gs = g.sort_values(n_col).reset_index(drop=True)
    N   = len(gs)
    vals = []
    for m in range(N, min_pts - 1, -1):
        w  = gs.iloc[N - m:].reset_index(drop=True)
        ns = float(n_scale_fn(w[n_col].values))
        best, _tab = select_k_aic(w, n_col, p, ns, k_max=k_max)
        if best is None:
            continue
        v = n_req_kterm_signed(best['fit']['beta'], tol, p, best['k'], ns)
        if np.isfinite(v) and v > 0:
            vals.append(v)
    if not vals:
        return float('nan')
    if len(vals) == 1:
        return 0.0
    return float(np.log(vals).std(ddof=1))


def sigma_jackknife(g, n_col, p, tol, n_scale_fn, k_max=K_MAX, min_pts=3):
    """σ_jack — the DATA-LIMITATION term (METHODS §6.3).

    Leave one abscissa out and redo the WHOLE procedure on the remainder —
    window search, AIC choice of k, crossing — then form the usual jackknife
    spread of ln n^req:

        σ_jack² = (N−1)/N · Σ_i ( ln n^req_(−i) − mean )²

    Answers the question the parametric MC cannot: how much would this move if
    the sweep had been sampled slightly differently? For a correctly specified
    model σ_jack ≈ σ_MC; when it is much larger, individual points are steering
    the extrapolation, which is exactly the regime these sweeps are in.

    Returns nan on fewer than 3 usable leave-one-out fits."""
    gs = g.sort_values(n_col).reset_index(drop=True)
    N  = len(gs)
    vals = []
    for i in range(N):
        sub = gs.drop(index=i).reset_index(drop=True)
        if len(sub) < min_pts:
            continue
        aw = fit_asymptotic_window(sub, n_col, p, n_scale_fn=n_scale_fn,
                                   k_max=k_max, min_pts=min_pts)
        if aw is None:
            continue
        v = n_req_kterm_signed(aw['best']['fit']['beta'], tol, p,
                               aw['best']['k'], aw['n_scale'])
        if np.isfinite(v) and v > 0:
            vals.append(v)
    if len(vals) < 3:
        return float('nan')
    a = np.log(vals)
    n = len(a)
    return float(np.sqrt((n - 1) / n * ((a - a.mean())**2).sum()))


def sigma_log_total(g, n_col, p_hat, sigma_p, tol, n_scale_fn,
                    k_max=K_MAX, rng_seed=42, min_pts=3):
    """The adopted error on ln n^req: three terms in quadrature (METHODS §6).

        σ_total = sqrt( σ_MC² + σ_win² + σ_jack² )

    σ_MC   (§6.1) parameter noise inside one model — k, β, and p when it is
                  fitted rather than derived. Computed on the ADOPTED window.
    σ_win  (§6.2) where the asymptotic range was cut.
    σ_jack (§6.3) which simulations happen to exist.

    The last two are computed on the FULL point set `g`, since both are about
    varying the point set. They are not strictly independent — both respond to
    leverage from the low-n end — so the quadrature sum is mildly conservative,
    and that is the intended direction.

    NOT included: the leading order p. One shared p moves every ε_y coherently,
    so it rotates the calibration line instead of scattering points about it, and
    it is quoted as an explicit systematic (§7.1) rather than buried here.

    Returns dict(mc, win, jack, total)."""
    aw = fit_asymptotic_window(g, n_col, p_hat, n_scale_fn=n_scale_fn,
                               k_max=k_max, min_pts=min_pts)
    if aw is None:
        return dict(mc=float('nan'), win=float('nan'), jack=float('nan'),
                    total=float('nan'))
    s_mc = mc_sigma_log_n_req(aw['window'], n_col, p_hat, sigma_p, tol,
                              aw['n_scale'], rng_seed=rng_seed, k_max=k_max)[0]
    s_win  = sigma_window(g, n_col, p_hat, tol, n_scale_fn, k_max, min_pts)
    s_jack = sigma_jackknife(g, n_col, p_hat, tol, n_scale_fn, k_max, min_pts)
    parts  = [x for x in (s_mc, s_win, s_jack) if np.isfinite(x)]
    total  = float(np.sqrt(sum(x**2 for x in parts))) if parts else float('nan')
    return dict(mc=s_mc, win=s_win, jack=s_jack, total=total)
