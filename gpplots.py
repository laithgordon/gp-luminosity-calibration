"""Figure output for GP_CALIBRATION.ipynb: where figures are written, and the one
sweep plot that is drawn the same way several times.
"""
from pathlib import Path

import matplotlib.pyplot as plt

PLOT_DIR = Path('plots')         # every figure is written here as PDF + PNG
PLOT_DIR.mkdir(exist_ok=True)
L_UNIT = 1e34                    # report L in units of 10^34 cm^-2 s^-1


def save_fig(fig, name, dpi=300):
    """Save `fig` to plots/<name>.pdf and plots/<name>.png (tight bbox)."""
    for ext in ('pdf', 'png'):
        fig.savefig(PLOT_DIR / f'{name}.{ext}', bbox_inches='tight', dpi=dpi)
    print(f'Saved: plots/{name}.pdf, plots/{name}.png')


def plot_sweep(g, sweep_col, nominal_val, xlabel, title):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = g[sweep_col].values
    y = g['mean_L'].values / L_UNIT
    yerr = g['std_L'].values / L_UNIT     # DISPLAY error = seed STD
    ax.errorbar(x, y, yerr=yerr, fmt='o-', capsize=3, color='C0',
                ecolor='C0', markersize=6, lw=1.4, label='mean ± STD')
    # annotate n_seeds for any point that isn't exactly 15
    for xi, yi, ns in zip(x, y, g['n_seeds'].values):
        if ns != 15:
            ax.annotate(f'n={ns}', (xi, yi), textcoords='offset points',
                        xytext=(0, 8), ha='center', fontsize=8, color='gray')
    ax.set_xscale('log', base=2)
    ax.xaxis.set_major_formatter(plt.matplotlib.ticker.FuncFormatter(
        lambda v, _: f'{int(v)}'))
    ax.axvline(nominal_val, color='k', ls=':', lw=1, alpha=0.6,
               label=f'nominal {sweep_col}={nominal_val}')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$L$  [$10^{34}\,$cm$^{-2}\,$s$^{-1}$]')
    ax.set_title(title)
    ax.grid(True, which='both', alpha=0.25)
    ax.legend(loc='best', fontsize=9)
    plt.tight_layout()
    return fig, ax
