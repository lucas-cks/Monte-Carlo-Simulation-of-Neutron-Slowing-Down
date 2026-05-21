import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("default")


# Data helpers

def load_case(filename: str, y_max: float = 9.5, flux_cut: float = 1e-8) -> pd.DataFrame:
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    df = pd.read_csv(filename)
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    required = {"y", "flux"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{filename} is missing columns: {sorted(missing)}")

    if "flux_err" not in df.columns:
        df["flux_err"] = 0.0

    df = df[(df["y"] > 0.0) & (df["y"] < y_max) & (df["flux"] > flux_cut)].copy()
    df = df.sort_values("y").reset_index(drop=True)

    df["X"] = 16 * df["y"]
    return df


def load_convergence(filename: str = "convergence_histories.csv") -> pd.DataFrame:
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    df = pd.read_csv(filename)
    required = {"histories", "peak_y", "peak_flux", "tail_slope"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{filename} is missing columns: {sorted(missing)}")

    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    df = df.sort_values("histories").reset_index(drop=True)
    return df


def mb_ref(y: np.ndarray) -> np.ndarray:
    """Normalized Maxwell-Boltzmann speed-shape reference: y^2 exp(-y^2)."""
    y = np.asarray(y, dtype=float)
    ref = y**2 * np.exp(-y**2)
    m = np.nanmax(ref)
    return ref / m if m > 0 else ref


def one_over_v_ref(y: np.ndarray, scale: float = 0.30) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    return scale / y


def calculate_ak_product(label: str) -> float:
    if "Hydrogen" in label and "K=0.18" in label:
        return 1.0 * 0.18
    if "High Abs" in label and "K=0.36" in label:
        return 1.0 * 0.36
    if "Carbon" in label and "K=0.03125" in label:
        return 12.0 * 0.03125
    if "Helium" in label and "K=0.03125" in label:
        return 2.0 * 0.03125
    if "Silicon" in label and "K=0.03125" in label:
        return 16.0 * 0.03125
    if "Manganese" in label and "K=0.03125" in label:
        return 25.0 * 0.03125
    return np.nan


def safe_peak(df: pd.DataFrame, y_lo: float = 0.5, y_hi: float = 2.5):
    region = df[(df["y"] >= y_lo) & (df["y"] <= y_hi)]
    if len(region) == 0:
        return np.nan, np.nan
    idx = region["flux"].idxmax()
    return float(region.loc[idx, "y"]), float(region.loc[idx, "flux"])


def safe_tail_slope(df: pd.DataFrame, y_lo: float = 3.0, y_hi: float = 7.0) -> float:
    tail = df[(df["y"] >= y_lo) & (df["y"] <= y_hi)]
    if len(tail) < 6:
        return np.nan

    x = np.log(tail["y"].to_numpy())
    y = np.log(tail["flux"].to_numpy())
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return np.nan

    denom = len(x) * np.sum(x * x) - np.sum(x) ** 2
    if np.abs(denom) < 1e-30:
        return np.nan

    slope = (len(x) * np.sum(x * y) - np.sum(x) * np.sum(y)) / denom
    return float(slope)


def sparse_for_errorbars(df: pd.DataFrame, max_points: int = 30) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    step = max(1, len(df) // max_points)
    return df.iloc[::step].copy()


# Main figure

def analyze_flux_results():
    files = {
        "Hydrogen (A=1, K=0.18)": "neutron_flux_hydrogen.csv",
        "High Abs (A=1, K=0.36)": "neutron_flux_highK.csv",
        "Carbon (A=12, K=0.03125)": "neutron_flux_carbon.csv",
        "Helium (A=2, K=0.03125)": "neutron_flux_helium.csv",
        "Silicon (A=16, K=0.03125)": "neutron_flux_silicon.csv",
        "Manganese (A=25, K=0.03125)": "neutron_flux_manganese.csv",
    }

    colors = {
        "Hydrogen (A=1, K=0.18)": "green",
        "High Abs (A=1, K=0.36)": "red",
        "Carbon (A=12, K=0.03125)": "blue",
        "Helium (A=2, K=0.03125)": "orange",
        "Silicon (A=16, K=0.03125)": "purple",
        "Manganese (A=25, K=0.03125)": "brown"
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    y_ref = np.linspace(0.01, 10.0, 800)
    mb = mb_ref(y_ref)

    y_tail = np.linspace(2.0, 8.0, 200)
    inv_v = one_over_v_ref(y_tail)

    print("Neutron Flux Simulation Analysis")

    peak_positions = {}
    tail_slopes = {}
    hardening_factors = {}

    for label, filename in files.items():
        try:
            df = load_case(filename)
        except FileNotFoundError:
            print(f"Missing file: {filename}")
            continue

        if len(df) == 0:
            print(f"No valid data for {label}")
            continue

        # Panel 1: log-log spectrum
        axes[0, 0].loglog(
            df["X"], df["flux"],
            color=colors[label], linewidth=2, label=label
        )

        err_df = sparse_for_errorbars(df, max_points=25)
        axes[0, 0].errorbar(
            err_df["X"], err_df["flux"],
            yerr=err_df["flux_err"],
            fmt="none",
            ecolor=colors[label],
            alpha=0.25,
            capsize=2,
            linewidth=1
        )

        # Panel 2: thermal region
        thermal_df = df[df["y"] <= 3.0].copy()
        axes[0, 1].errorbar(
            thermal_df["X"], thermal_df["flux"],
            yerr=thermal_df["flux_err"],
            color=colors[label],
            linewidth=2,
            capsize=2,
            label=label
        )

        # Panel 3: deviation from MB
        interp_flux = np.interp(y_ref, df["y"].to_numpy(), df["flux"].to_numpy(), left=np.nan, right=np.nan)
        deviation = interp_flux - mb
        axes[1, 0].plot(16 * y_ref, deviation, color=colors[label], linewidth=2, label=label)

        # Metrics
        peak_y, peak_flux = safe_peak(df, 0.5, 2.5)
        slope = safe_tail_slope(df, 3.0, 7.0)

        peak_positions[label] = peak_y
        tail_slopes[label] = slope

        print(
            f"{label:<28}"
            f"Peak y = {peak_y:.3f} | "
            f"Tail slope = {slope:.3f}"
        )

    # Panel 1 formatting
    X_tail = 16 * y_tail
    axes[0, 0].loglog(X_tail, inv_v, "k--", linewidth=2, alpha=0.7, label="Theoretical 1/v")
    axes[0, 0].set_title("Neutron Flux Spectrum (Log-Log Scale)", fontsize=14, fontweight="bold")
    axes[0, 0].set_xlabel(r"Variable $X=16(v/v_T)$", fontsize=12)
    axes[0, 0].set_ylabel(r"Normalized Flux $\phi(X)$", fontsize=12)
    axes[0, 0].grid(True, which="both", alpha=0.3)
    axes[0, 0].legend(fontsize=10)
    axes[0, 0].set_xlim(1.6, 160)
    axes[0, 0].set_ylim(1e-3, 10)

    # Panel 2 formatting
    axes[0, 1].plot(16* y_ref, mb, "k:", linewidth=3, alpha=0.7, label="Maxwell-Boltzmann")
    axes[0, 1].set_title("Thermal Peak Region (Linear Scale)", fontsize=14, fontweight="bold")
    axes[0, 1].set_xlabel(r"Variable $X=16(v/v_T)$", fontsize=12)
    axes[0, 1].set_ylabel(r"Normalized Flux $\phi(X)$", fontsize=12)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=10)
    axes[0, 1].set_xlim(0.0, 40.0)
    axes[0, 1].set_ylim(0.0, 1.2)

    # Panel 3 formatting
    axes[1, 0].axhline(0.0, color="k", linestyle=":", alpha=0.5)
    axes[1, 0].set_title("Deviation from Maxwell-Boltzmann", fontsize=14, fontweight="bold")
    axes[1, 0].set_xlabel(r"Variable $X=16(v/v_T)$", fontsize=12)
    axes[1, 0].set_ylabel(r"$\phi_{\mathrm{sim}}-\phi_{\mathrm{MB}}$", fontsize=12)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=10)
    axes[1, 0].set_xlim(0.0, 48.0)

    # Panel 4: summary table
    axes[1, 1].axis("off")
    table_data = []
    for label in files.keys():
        if label in peak_positions and np.isfinite(peak_positions[label]):
            ak = calculate_ak_product(label)
            table_data.append([
                label,
                f"{peak_positions[label]:.3f}",
                f"{tail_slopes[label]:.3f}",
                f"{ak:.3f}",
            ])

    if table_data:
        # FIX: colLabels now has 4 entries, matching the 4 columns in table_data
        table = axes[1, 1].table(
            cellText=table_data,
            colLabels=["Case", "Peak y", "Tail Slope", "A×K"],
            cellLoc="center",
            loc="center",
            bbox=[0.06, 0.25, 0.88, 0.62],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)

    axes[1, 1].set_title("Spectral Hardening Analysis", fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig("comprehensive_flux_analysis.png", dpi=300, bbox_inches="tight")
    plt.show()

    print("\nConsistency Check")


    for label in files.keys():
        if label not in peak_positions or not np.isfinite(peak_positions[label]):
            continue

        ak = calculate_ak_product(label)
        expected = np.sqrt(1.0 + 1.11 * ak)

        print(f"\n{label}")
        print(f"  A×K                         = {ak:.3f}")
        print(f"  Expected y_peak             = {expected:.3f}")
        print(f"  Measured y_peak             = {peak_positions[label]:.3f}")
        print(f"  Deviation percentage        = {(peak_positions[label] - expected) / expected * 100:.1f}%")
        print(f"  Tail slope                  = {tail_slopes[label]:.3f}")

        if np.isfinite(tail_slopes[label]) and abs(tail_slopes[label] + 1.0) < 0.15:
            print("  1/v tail behavior          : GOOD")
        else:
            print("  1/v tail behavior          : WEAK")



# Convergence study

def plot_convergence_study(filename: str = "convergence_histories.csv"):
    try:
        df = load_convergence(filename)
    except FileNotFoundError:
        print(f"Convergence file not found: {filename}")
        return
    except ValueError as e:
        print(f"Convergence file invalid: {e}")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    histories = df["histories"].to_numpy()

    axes[0].plot(histories, df["peak_y"], marker="o", linewidth=2)
    axes[0].set_xscale("log")
    axes[0].set_title("Peak Position vs Histories")
    axes[0].set_xlabel("Histories")
    axes[0].set_ylabel("Peak y")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(histories, df["peak_flux"], marker="o", linewidth=2)
    axes[1].set_xscale("log")
    axes[1].set_title("Peak Flux vs Histories")
    axes[1].set_xlabel("Histories")
    axes[1].set_ylabel("Peak Flux")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(histories, df["tail_slope"], marker="o", linewidth=2, label="Measured")
    axes[2].axhline(-1.0, color="k", linestyle="--", alpha=0.7, label="1/v slope")
    axes[2].set_xscale("log")
    axes[2].set_title("Tail Slope vs Histories")
    axes[2].set_xlabel("Histories")
    axes[2].set_ylabel("Tail Slope")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()

    plt.tight_layout()
    plt.savefig("convergence_study.png", dpi=300, bbox_inches="tight")
    plt.show()


# Individual comparison figures

def plot_individual_comparisons():
    files = {
        "Hydrogen (A=1, K=0.18)": "neutron_flux_hydrogen.csv",
        "High Abs (A=1, K=0.36)": "neutron_flux_highK.csv",
        "Carbon (A=12, K=0.03125)": "neutron_flux_carbon.csv",
        "Helium (A=2, K=0.03125)": "neutron_flux_helium.csv",
        "Silicon (A=16, K=0.03125)": "neutron_flux_silicon.csv",
        "Manganese (A=25, K=0.03125)": "neutron_flux_manganese.csv",
    }

    colors = {
        "Hydrogen (A=1, K=0.18)": "green",
        "High Abs (A=1, K=0.36)": "red",
        "Carbon (A=12, K=0.03125)": "blue",
        "Helium (A=2, K=0.03125)": "orange",
        "Silicon (A=16, K=0.03125)": "purple",
        "Manganese (A=25, K=0.03125)": "brown"
    }

    n_files = len(files)
    n_cols = min(3, n_files)
    n_rows = (n_files + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    y_ref = np.linspace(0.01, 3.0, 400)
    mb = mb_ref(y_ref)

    for idx, (label, filename) in enumerate(files.items()):
        ax = axes[idx]
        try:
            df = load_case(filename)
        except FileNotFoundError:
            ax.text(
                0.5, 0.5, f"File not found:\n{filename}",
                ha="center", va="center", transform=ax.transAxes
            )
            ax.set_title(label)
            continue

        if len(df) == 0:
            ax.text(
                0.5, 0.5, f"No valid data\nfor {label}",
                ha="center", va="center", transform=ax.transAxes
            )
            ax.set_title(label)
            continue

        ax.plot(df["X"], df["flux"], color=colors[label], linewidth=2, label="Simulation")

        err_df = sparse_for_errorbars(df, max_points=20)
        ax.errorbar(
            err_df["X"], err_df["flux"],
            yerr=err_df["flux_err"],
            fmt="none",
            ecolor=colors[label],
            alpha=0.25,
            capsize=2
        )

        ax.plot(16* y_ref, mb, "k--", linewidth=2, alpha=0.7, label="Maxwell-Boltzmann")

        peak_y, peak_flux = safe_peak(df, 0.5, 2.5)
        
        peak_X = 16 * peak_y
        if np.isfinite(peak_y):
            ax.axvline(x=peak_X, color=colors[label], linestyle=":", alpha=0.7)
            ax.plot(peak_X, peak_flux, "o", color=colors[label], markersize=8)

        ax.set_xlim(0.0, 32.0)
        ax.set_ylim(0.0, 1.1)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel(r"Variable $X=16(v/v_T)$", fontsize=11)
        ax.set_ylabel("Normalized Flux", fontsize=11)
        ax.set_title(f"{label}\nPeak at X = {peak_X:.2f}", fontsize=12)
        ax.legend()

    # Hide any unused subplots
    for idx in range(len(files), len(axes)):
        axes[idx].axis("off")

    plt.tight_layout()
    plt.savefig("individual_flux_comparisons.png", dpi=300, bbox_inches="tight")
    plt.show()


# Simple comparison

def plot_simple_comparison():
    files = {
        "Hydrogen (A=1, K=0.18)": "neutron_flux_hydrogen.csv",
        "High Abs (A=1, K=0.36)": "neutron_flux_highK.csv",
        "Carbon (A=12, K=0.03125)": "neutron_flux_carbon.csv",
        "Helium (A=2, K=0.03125)": "neutron_flux_helium.csv",
        "Silicon (A=16, K=0.03125)": "neutron_flux_silicon.csv",
        "Manganese (A=25, K=0.03125)": "neutron_flux_manganese.csv",
    }

    colors = {
        "Hydrogen (A=1, K=0.18)": "green",
        "High Abs (A=1, K=0.36)": "red",
        "Carbon (A=12, K=0.03125)": "blue",
        "Helium (A=2, K=0.03125)": "orange",
        "Silicon (A=16, K=0.03125)": "purple",
        "Manganese (A=25, K=0.03125)": "brown"
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    y_ref = np.linspace(0.01, 10.0, 600)
    mb = mb_ref(y_ref)

    print("Simple Flux Analysis")

    for label, filename in files.items():
        try:
            df = load_case(filename)
        except FileNotFoundError:
            print(f"Missing file: {filename}")
            continue

        if len(df) == 0:
            print(f"No valid data: {label}")
            continue

        ax1.loglog(df["X"], df["flux"], color=colors[label], linewidth=2, label=label)

        err_df = sparse_for_errorbars(df, max_points=20)
        ax1.errorbar(
            err_df["X"], err_df["flux"],
            yerr=err_df["flux_err"],
            fmt="none",
            ecolor=colors[label],
            alpha=0.2,
            capsize=2
        )

        thermal_df = df[df["X"] <= 40.0]
        ax2.errorbar(
            thermal_df["X"], thermal_df["flux"],
            yerr=thermal_df["flux_err"],
            color=colors[label],
            linewidth=2,
            capsize=2,
            label=label
        )

        peak_y, _ = safe_peak(thermal_df, 0.5, 2.5)
        print(f"{label:<28} Peak y = {peak_y:.3f}")

    y_tail = np.linspace(2.0, 8.0, 200)
    ax1.loglog(y_tail, 0.30 / y_tail, "k--", alpha=0.7, label="1/v theory")
    ax1.loglog(y_ref, mb, "k:", alpha=0.7, label="Maxwell-Boltzmann")
    ax2.plot(y_ref, mb, "k:", linewidth=2, alpha=0.7, label="Maxwell-Boltzmann")

    ax1.set_title("Neutron Flux Spectrum (Log-Log)")
    ax1.set_xlabel(r"Variable $X=16(v/v_T)$")
    ax1.set_ylabel(r"Flux $\phi(X)$")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    ax1.set_xlim(1.6, 160)
    ax1.set_ylim(1e-3, 2)

    ax2.set_title("Thermal Peak Region")
    ax2.set_xlabel(r"Variable $X=16(v/v_T)$")
    ax2.set_ylabel(r"Flux $\phi(X)$")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0.0, 40.0)
    ax2.set_ylim(0.0, 1.2)

    plt.tight_layout()
    plt.savefig("simple_flux_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()


# Main

if __name__ == "__main__":
    print("Neutron Flux Simulation Analysis")
    print("Coveyou, Bate & Osborn (1956)")
    print()

    try:
        analyze_flux_results()
    except Exception as e:
        print(f"Comprehensive analysis failed: {e}")
        print("Falling back to simple comparison")
        plot_simple_comparison()

    try:
        plot_individual_comparisons()
    except Exception as e:
        print(f"Individual comparisons failed: {e}")

    try:
        plot_convergence_study()
    except Exception as e:
        print(f"Convergence study failed: {e}")

    print("\nAnalysis complete.")