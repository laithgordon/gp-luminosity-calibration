"""Conservative-locus constants -- FROZEN.

These numbers define the conservative production runs (the `conservative` block
of data/gp_luminosity_for_wx.csv) and tab:ny_recommendations. They were fixed
when the conservative decks were written and the runs submitted (decks
2026-09-11, runs 2026-09-12) and must NOT be recomputed from the current fit.
The fitted constants in data/gp_constants_export.csv (C_y_fit = 26.2,
q_n_fit = 0.402 +- 0.152, C_m_fit, s_fit) are a separate, descriptive result
that moves if the ladder data change; this locus does not.

    n_y_cons = 2^ceil(log2(C_Y_CONS * D_y^Q_N_CONS))                  rounded up
    n_m_cons = C_M_CONS * D_y^NM_EXPONENT_CONS * n_x * n_y_cons * n_z

NM_EXPONENT_CONS is the effective exponent s_cons - q_p = 3.300. It, not s_cons,
defines the runs; with Q_P below, s_cons = NM_EXPONENT_CONS + Q_P = 3.5153.
D_y per eps_y is data/D_y_table.json.
"""
import math

C_Y_CONS = 50.0            # n_y normalisation: a chosen margin, ~1.9x the fitted C_y
Q_N_CONS = 0.402           # n_y exponent
C_M_CONS = 2.305e-7        # n_m normalisation
NM_EXPONENT_CONS = 3.300   # effective n_m exponent, s_cons - q_p
N_X_CONS = 512             # transverse cells held fixed on the locus
N_Z_CONS = 64              # longitudinal slices held fixed on the locus

# q_p is the envelope-integration (first-waist) pinch exponent, shared with the WarpX
# repository; the predicted n_y exponent q_n^pred = q_p + 1/4 is derived from it.
Q_P = 0.2153
Q_N_PRED = Q_P + 0.25


def n_y_cons(d_y: float) -> int:
    """Vertical cells on the conservative locus (power of two, rounded up)."""
    return 2 ** math.ceil(math.log2(C_Y_CONS * d_y ** Q_N_CONS))


def n_m_cons(d_y: float, n_y: int, n_x: int = N_X_CONS, n_z: int = N_Z_CONS) -> float:
    """Macroparticles per beam on the conservative locus (not rounded)."""
    return C_M_CONS * d_y ** NM_EXPONENT_CONS * n_x * n_y * n_z
