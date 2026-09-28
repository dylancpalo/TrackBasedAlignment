import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

if len(sys.argv) != 2:
    raise SystemExit("Usage: python plot_DOCA_XYV_folder_input.py CSV_<name>_<name2>")

csv_folder = Path(sys.argv[1])
if not csv_folder.is_dir():
    raise SystemExit(f"CSV folder does not exist: {csv_folder}")

def find_one(prefix):
    matches = sorted(csv_folder.glob(f"{prefix}*.csv"))
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one {prefix}*.csv in {csv_folder}; found {len(matches)}")
    return matches[0]

doca_file = find_one("DOCA_Histograms_")
x_file = find_one("XErrorPanel_")
y_file = find_one("YErrorPanel_")
v_file = find_one("VErrorPanel_")

print(f"CSV folder: {csv_folder}")
print(f"DOCA input: {doca_file.name}")
print(f"X input: {x_file.name}")
print(f"Y input: {y_file.name}")
print(f"V/dV input: {v_file.name}")

N_STATIONS = 18
PANELS_PER_STATION = 12
A2_OVER_A1 = 0.2
POINT_SIZE = 10

DOCA_PANELS = [0, 1, 2, 3, 4, 5]
X_PANELS = [1, 4, 7, 10]
Y_PANELS = [0, 5, 6, 11]
V_PANELS = list(range(PANELS_PER_STATION))

DOCA_LIMITS = {
    "sigma1": (0, 500),
    "sigma2": (00, 1000),
    "mu1": (-200, 200),
    "mu2": (-200, 200),
}

PARAMETERS = [
    ("sigma1", r"$\sigma_1$", "Gaussian width [um]"),
    ("sigma2", r"$\sigma_2$", "Gaussian width [um]"),
    ("mu1", r"$\mu_1$", "Gaussian mean [um]"),
    ("mu2", r"$\mu_2$", "Gaussian mean [um]"),
]

def load_csv(filename):
    raw = pd.read_csv(filename, header=None, skipinitialspace=True)
    if raw.shape[1] < 6:
        raise ValueError(f"{filename}: expected at least 6 columns")
    df = raw.iloc[:, :6].copy()
    df.columns = ["x", "A1", "sigma1", "sigma2", "mu1", "mu2"]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[np.isfinite(df["x"])].copy()
    idx = df["x"].astype(int)
    df["station"] = idx // 12
    df["panel_in_station"] = idx % 12
    df["panel_in_plane"] = idx % 6
    return df

def good_doca(df):
    fitted_entries = (
        np.sqrt(2*np.pi) * df["A1"]
        * (df["sigma1"] + 0.2*df["sigma2"]) / 40.0
    )
    low = fitted_entries < 1000.0
    repeat = (
        df["mu1"].eq(df["mu1"].shift(1))
        & df["sigma1"].eq(df["sigma1"].shift(1))
    )
    finite = np.isfinite(
        df[["A1","sigma1","sigma2","mu1","mu2"]]
    ).all(axis=1)
    return finite & ~low & ~repeat

def good_xy(df):
    A2 = 0.2 * df["A1"]
    integrated = (
        np.sqrt(2*np.pi) * df["A1"] * df["sigma1"]
        + np.sqrt(2*np.pi) * A2 * df["sigma2"]
    )
    fitcols = ["A1","sigma1","sigma2","mu1","mu2"]
    repeat = pd.Series(False, index=df.index)
    if len(df) > 1:
        a = df[fitcols].iloc[1:].reset_index(drop=True)
        b = df[fitcols].iloc[:-1].reset_index(drop=True)
        repeat.iloc[1:] = a.eq(b).all(axis=1).to_numpy()
    finite = np.isfinite(df[fitcols]).all(axis=1)
    return finite & (integrated >= 500.0) & (df["sigma1"] >= 152.0) & ~repeat

def colors(n):
    return plt.rcParams["axes.prop_cycle"].by_key()["color"][:n]

def scatter_by_panel(ax, df, good, col, panel_values, panel_col, cs):
    for c, p in zip(cs, panel_values):
        m = good & (df[panel_col] == p)
        ax.scatter(df.loc[m,"x"], df.loc[m,col], s=POINT_SIZE, color=c)

def make_v1(df, good, label, panel_values, panel_col, limits=None):
    cs = colors(len(panel_values))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (col, title, ylabel) in zip(axes.ravel(), PARAMETERS):
        scatter_by_panel(ax, df, good, col, panel_values, panel_col, cs)
        for station in range(N_STATIONS):
            m = good & (df["station"] == station) & df[panel_col].isin(panel_values)
            vals = df.loc[m, col].dropna()
            if len(vals):
                ax.hlines(vals.mean(), station*12, station*12+11,
                          color="black", linewidth=1.6)
        ax.set_xlim(0, 215)
        if limits and col in limits:
            ax.set_ylim(*limits[col])
        ax.set_xlabel("Index")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25)
    handles = [
        Line2D([0],[0], marker="o", linestyle="None", color=c,
               label=f"Panel {p}", markersize=5)
        for c,p in zip(cs,panel_values)
    ] + [Line2D([0],[0], color="black", label="Station average")]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1,0.5), fontsize=9)
    fig.suptitle(f"{label} Double-Gaussian Fit Parameters — Station Averages", fontsize=14)
    fig.tight_layout(rect=(0,0,0.89,1))
    return fig

def make_v2(df, good, label, panel_values, panel_col, limits=None):
    cs = colors(len(panel_values))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (col, title, ylabel) in zip(axes.ravel(), PARAMETERS):
        scatter_by_panel(ax, df, good, col, panel_values, panel_col, cs)
        for c, p in zip(cs, panel_values):
            m = good & (df[panel_col] == p)
            vals = df.loc[m, col].dropna()
            if len(vals):
                ax.axhline(vals.mean(), color=c, linewidth=1.5)
        ax.set_xlim(0, 215)
        if limits and col in limits:
            ax.set_ylim(*limits[col])
        ax.set_xlabel("Index")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25)
    handles = [
        Line2D([0],[0], marker="o", linestyle="-", color=c,
               label=f"Panel {p} average", markersize=5)
        for c,p in zip(cs,panel_values)
    ]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1,0.5), fontsize=9)
    fig.suptitle(f"{label} Double-Gaussian Fit Parameters — Panel Averages", fontsize=14)
    fig.tight_layout(rect=(0,0,0.89,1))
    return fig

def run_one(filename, label, panel_values, panel_col, good_func, limits=None):
    df = load_csv(filename)
    good = good_func(df) & df[panel_col].isin(panel_values)

    fig1 = make_v1(df, good, label, panel_values, panel_col, limits)
    fig2 = make_v2(df, good, label, panel_values, panel_col, limits)

    stem = filename.stem
    out1 = filename.with_name(stem + "_V1_station_averages.pdf")
    out2 = filename.with_name(stem + "_V2_panel_averages.pdf")
    fig1.savefig(out1, bbox_inches="tight")
    fig2.savefig(out2, bbox_inches="tight")
    plt.close(fig1)
    plt.close(fig2)

    print(f"{label}:")
    print(f"  accepted = {good.sum()} / {len(df)}")
    print(f"  {out1}")
    print(f"  {out2}")

# Produces 8 PDFs total: station-average and panel-average summaries
# for DOCA, X, Y, and V/dV. Each PDF contains sigma1, sigma2, mu1, mu2.

run_one(doca_file, "DOCA", DOCA_PANELS, "panel_in_plane", good_doca, DOCA_LIMITS)
run_one(x_file, "X Error", X_PANELS, "panel_in_station", good_xy, None)
run_one(y_file, "Y Error", Y_PANELS, "panel_in_station", good_xy, None)
run_one(v_file, "V/dV Error", V_PANELS, "panel_in_station", good_xy, None)
