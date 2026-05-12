import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def analyze_flux_results():
    """Comprehensive analysis of neutron flux simulation results"""
    
    # File names and labels
    files = {
        'Hydrogen (A=1, K=0.18)': 'neutron_flux_hydrogen.csv',
        'Carbon (A=12, K=0.1)': 'neutron_flux_carbon.csv', 
        'High Abs (A=1, K=0.36)': 'neutron_flux_highK.csv'
    }
    
    colors = {
        'Hydrogen (A=1, K=0.18)': 'green',
        'Carbon (A=12, K=0.1)': 'blue',
        'High Abs (A=1, K=0.36)': 'red'
    }
    
    # Initialize plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Maxwell-Boltzmann reference
    y_ref = np.linspace(0.01, 10, 500)
    mb_flux = y_ref**2 * np.exp(-y_ref**2)
    mb_flux_normalized = mb_flux / np.max(mb_flux)
    
    # Theoretical 1/v reference
    y_tail = np.linspace(2, 8, 100)
    one_over_v = 0.3 / y_tail  # Arbitrary normalization
    
    print("=" * 70)
    print("NEUTRON FLUX ANALYSIS - Comparison with 1956 Paper")
    print("=" * 70)
    
    # Analysis for each case
    peak_positions = {}
    tail_slopes = {}
    hardening_factors = {}
    
    for label, filename in files.items():
        try:
            df = pd.read_csv(filename)
            
            # Clean data - remove source spike and zero fluxes
            clean_df = df[(df['y'] < 9.5) & (df['flux'] > 1e-4)].copy()
            
            if len(clean_df) == 0:
                print(f"Warning: No valid data for {label}")
                continue
                
            # 1. LOG-LOG PLOT (Tail behavior)
            axes[0,0].loglog(clean_df['y'], clean_df['flux'], 
                           label=label, color=colors[label], linewidth=2, marker='o', markersize=3)
            
            # 2. LINEAR PLOT (Thermal region)
            thermal_df = clean_df[clean_df['y'] <= 3.0]
            axes[0,1].plot(thermal_df['y'], thermal_df['flux'], 
                         label=label, color=colors[label], linewidth=2)
            
            # 3. DEVIATION FROM MAXWELL-BOLTZMANN
            # Interpolate simulation data to reference grid
            sim_flux_interp = np.interp(y_ref, clean_df['y'], clean_df['flux'], 
                                      left=np.nan, right=np.nan)
            deviation = sim_flux_interp - mb_flux_normalized
            
            axes[1,0].plot(y_ref, deviation, label=label, color=colors[label], linewidth=2)
            
            # Calculate statistics
            peak_idx = thermal_df['flux'].idxmax()
            peak_y = thermal_df.loc[peak_idx, 'y']
            peak_flux = thermal_df.loc[peak_idx, 'flux']
            
            # Tail slope analysis (y > 3)
            tail_data = clean_df[(clean_df['y'] >= 3.0) & (clean_df['y'] <= 7.0)]
            if len(tail_data) > 5:
                log_y = np.log(tail_data['y'])
                log_flux = np.log(tail_data['flux'])
                slope, intercept = np.polyfit(log_y, log_flux, 1)
                tail_slope = slope
            else:
                tail_slope = np.nan
            
            # Hardening factor (peak position relative to MB peak at y=1)
            hardening = peak_y - 1.0
            
            peak_positions[label] = peak_y
            tail_slopes[label] = tail_slope
            hardening_factors[label] = hardening
            
            print(f"{label:<25} | Peak: {peak_y:.3f} | Hardening: {hardening:+.3f} | Tail slope: {tail_slope:.3f}")
            
        except FileNotFoundError:
            print(f"Error: Could not find {filename}")
            continue
    
    # Plot formatting and reference curves
    
    # 1. Log-log plot
    axes[0,0].loglog(y_tail, one_over_v, 'k--', label='Theoretical 1/v', alpha=0.7, linewidth=2)
    axes[0,0].set_title('Neutron Flux Spectrum (Log-Log Scale)', fontsize=14, fontweight='bold')
    axes[0,0].set_xlabel('Speed Ratio $y = v/v_T$', fontsize=12)
    axes[0,0].set_ylabel('Normalized Flux $\\phi(y)$', fontsize=12)
    axes[0,0].grid(True, which="both", alpha=0.3)
    axes[0,0].legend(fontsize=10)
    axes[0,0].set_xlim(0.1, 10)
    axes[0,0].set_ylim(1e-3, 2)
    
    # 2. Linear thermal region plot
    axes[0,1].plot(y_ref, mb_flux_normalized, 'k:', label='Maxwell-Boltzmann', linewidth=3, alpha=0.7)
    axes[0,1].set_title('Thermal Peak Region (Linear Scale)', fontsize=14, fontweight='bold')
    axes[0,1].set_xlabel('Speed Ratio $y = v/v_T$', fontsize=12)
    axes[0,1].set_ylabel('Normalized Flux $\\phi(y)$', fontsize=12)
    axes[0,1].set_xlim(0, 2.5)
    axes[0,1].set_ylim(0, 1.2)
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].legend(fontsize=10)
    
    # 3. Deviation from Maxwell-Boltzmann
    axes[1,0].axhline(y=0, color='k', linestyle=':', alpha=0.5)
    axes[1,0].set_title('Deviation from Maxwell-Boltzmann', fontsize=14, fontweight='bold')
    axes[1,0].set_xlabel('Speed Ratio $y = v/v_T$', fontsize=12)
    axes[1,0].set_ylabel('$\\phi_{sim} - \\phi_{MB}$', fontsize=12)
    axes[1,0].set_xlim(0, 3)
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].legend(fontsize=10)
    
    # 4. Hardening analysis table
    axes[1,1].axis('off')
    
    # Create summary table
    table_data = []
    for label in files.keys():
        if label in peak_positions:
            ak_product = calculate_ak_product(label)
            table_data.append([
                label,
                f"{peak_positions[label]:.3f}",
                f"{hardening_factors[label]:+.3f}",
                f"{tail_slopes[label]:.3f}",
                f"{ak_product:.3f}"
            ])
    
    if table_data:
        table = axes[1,1].table(
            cellText=table_data,
            colLabels=['Case', 'Peak y', 'Hardening', 'Tail Slope', 'A×K'],
            cellLoc='center',
            loc='center',
            bbox=[0.1, 0.3, 0.8, 0.6]
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        axes[1,1].set_title('Spectral Hardening Analysis', fontsize=14, fontweight='bold', y=0.9)
    
    plt.tight_layout()
    plt.savefig('comprehensive_flux_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print detailed analysis
    print("\n" + "=" * 70)
    print("PHYSICS CONSISTENCY CHECK")
    print("=" * 70)
    
    # Check consistency with paper's predictions
    for label in files.keys():
        if label in peak_positions:
            ak = calculate_ak_product(label)
            expected_hardening = 1.11 * ak  # From paper's formula
            
            print(f"\n{label}:")
            print(f"  A×K = {ak:.3f}")
            print(f"  Expected hardening (1.11×A×K) = {expected_hardening:.3f}")
            print(f"  Actual peak shift = {hardening_factors[label]:.3f}")
            print(f"  Tail slope (should be ≈ -1.0) = {tail_slopes[label]:.3f}")
            
            if abs(tail_slopes[label] + 1.0) < 0.1:
                print("  ✓ 1/v tail behavior: GOOD")
            else:
                print("  ✗ 1/v tail behavior: POOR")
                
            # Check if hardening trend is physically reasonable
            if ak > 0.3 and hardening_factors[label] > 0.1:
                print("  ✓ Hardening trend: REASONABLE")
            elif ak < 0.2 and abs(hardening_factors[label]) < 0.2:
                print("  ✓ Hardening trend: REASONABLE") 
            else:
                print("  ✗ Hardening trend: QUESTIONABLE")

def calculate_ak_product(label):
    """Calculate A×K product from case label"""
    if 'Hydrogen' in label and 'K=0.18' in label:
        return 1.0 * 0.18
    elif 'Carbon' in label and 'K=0.1' in label:
        return 12.0 * 0.1
    elif 'High Abs' in label and 'K=0.36' in label:
        return 1.0 * 0.36
    else:
        return 0.0

def plot_individual_comparisons():
    """Plot individual cases with Maxwell-Boltzmann comparison"""
    
    files = {
        'Hydrogen (A=1, K=0.18)': 'neutron_flux_hydrogen.csv',
        'Carbon (A=12, K=0.1)': 'neutron_flux_carbon.csv',
        'High Abs (A=1, K=0.36)': 'neutron_flux_highK.csv'
    }
    
    colors = {'Hydrogen (A=1, K=0.18)': 'green', 
              'Carbon (A=12, K=0.1)': 'blue', 
              'High Abs (A=1, K=0.36)': 'red'}
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Maxwell-Boltzmann reference
    y_ref = np.linspace(0.01, 3, 300)
    mb_flux = y_ref**2 * np.exp(-y_ref**2)
    mb_flux_normalized = mb_flux / np.max(mb_flux)
    
    for idx, (label, filename) in enumerate(files.items()):
        try:
            df = pd.read_csv(filename)
            clean_df = df[(df['y'] < 9) & (df['flux'] > 1e-4)]
            
            axes[idx].plot(clean_df['y'], clean_df['flux'], 
                         color=colors[label], linewidth=2, label='Simulation')
            axes[idx].plot(y_ref, mb_flux_normalized, 'k--', 
                         linewidth=2, alpha=0.7, label='Maxwell-Boltzmann')
            
            # Mark peak position
            peak_idx = clean_df[clean_df['y'] <= 2]['flux'].idxmax()
            peak_y = clean_df.loc[peak_idx, 'y']
            peak_flux = clean_df.loc[peak_idx, 'flux']
            
            axes[idx].axvline(x=peak_y, color=colors[label], linestyle=':', alpha=0.7)
            axes[idx].plot(peak_y, peak_flux, 'o', color=colors[label], markersize=8)
            
            axes[idx].set_xlim(0, 2)
            axes[idx].set_ylim(0, 1.1)
            axes[idx].set_xlabel('Speed Ratio $y = v/v_T$', fontsize=11)
            axes[idx].set_ylabel('Normalized Flux', fontsize=11)
            axes[idx].set_title(f'{label}\nPeak at y = {peak_y:.2f}', fontsize=12)
            axes[idx].grid(True, alpha=0.3)
            axes[idx].legend()
            
        except FileNotFoundError:
            axes[idx].text(0.5, 0.5, f'File not found:\n{filename}', 
                         ha='center', va='center', transform=axes[idx].transAxes)
            axes[idx].set_title(label)
    
    plt.tight_layout()
    plt.savefig('individual_flux_comparisons.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_simple_comparison():
    """Simple comparison plot without complex analysis"""
    
    files = {
        'Hydrogen (A=1, K=0.18)': 'neutron_flux_hydrogen.csv',
        'Carbon (A=12, K=0.1)': 'neutron_flux_carbon.csv',
        'High Abs (A=1, K=0.36)': 'neutron_flux_highK.csv'
    }
    
    colors = {'Hydrogen (A=1, K=0.18)': 'green', 
              'Carbon (A=12, K=0.1)': 'blue', 
              'High Abs (A=1, K=0.36)': 'red'}
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Maxwell-Boltzmann reference
    y_ref = np.linspace(0.01, 10, 500)
    mb_flux = y_ref**2 * np.exp(-y_ref**2)
    mb_flux_normalized = mb_flux / np.max(mb_flux)
    
    print("Simple Flux Analysis:")
    print("-" * 50)
    
    for label, filename in files.items():
        try:
            df = pd.read_csv(filename)
            clean_df = df[(df['y'] < 9) & (df['flux'] > 1e-4)]
            
            # Log-log plot
            ax1.loglog(clean_df['y'], clean_df['flux'], label=label, 
                     color=colors[label], linewidth=2)
            
            # Linear plot (thermal region)
            thermal_df = clean_df[clean_df['y'] <= 3.0]
            ax2.plot(thermal_df['y'], thermal_df['flux'], label=label,
                   color=colors[label], linewidth=2)
            
            # Calculate peak
            peak_idx = thermal_df['flux'].idxmax()
            peak_y = thermal_df.loc[peak_idx, 'y']
            
            print(f"{label:<25} | Peak at y = {peak_y:.3f}")
            
        except FileNotFoundError:
            print(f"File not found: {filename}")
            continue
    
    # Add references
    y_tail = np.linspace(2, 8, 100)
    one_over_v = 0.3 / y_tail
    ax1.loglog(y_tail, one_over_v, 'k--', label='1/v theory', alpha=0.7)
    ax1.loglog(y_ref, mb_flux_normalized, 'k:', label='Maxwell-Boltzmann', alpha=0.7)
    
    ax2.plot(y_ref, mb_flux_normalized, 'k:', label='Maxwell-Boltzmann', linewidth=2, alpha=0.7)
    
    # Formatting
    ax1.set_title('Neutron Flux Spectrum (Log-Log)')
    ax1.set_xlabel('Speed Ratio $y = v/v_T$')
    ax1.set_ylabel('Flux $\\phi(y)$')
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    ax1.set_xlim(0.1, 10)
    ax1.set_ylim(1e-3, 2)
    
    ax2.set_title('Thermal Peak Region')
    ax2.set_xlabel('Speed Ratio $y = v/v_T$')
    ax2.set_ylabel('Flux $\\phi(y)$')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 2.5)
    ax2.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig('simple_flux_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("Neutron Flux Simulation Analysis")
    print("Based on: Coveyou, Bate & Osborn (1956) - Effect of Moderator Temperature")
    print()
    
    try:
        # Try comprehensive analysis first
        analyze_flux_results()
    except Exception as e:
        print(f"Comprehensive analysis failed: {e}")
        print("Falling back to simple comparison...")
        plot_simple_comparison()
    
    # Always try individual comparisons
    try:
        plot_individual_comparisons()
    except Exception as e:
        print(f"Individual comparisons failed: {e}")
    
    print("\n" + "=" * 70)
    print("Analysis complete. Check generated plots.")
    print("=" * 70)