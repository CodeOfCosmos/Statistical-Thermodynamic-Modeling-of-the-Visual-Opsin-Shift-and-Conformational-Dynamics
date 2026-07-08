import numpy as np
import matplotlib.pyplot as plt

def simulate_thermodynamic_shift():
    
    # 1. DEFINE SYSTEM CONFIGURATIONS & REACTION COORDINATE
  
    # xi represents the reaction coordinate (0 = Inactive Rhodopsin, 1 = Active Metarhodopsin II)
    xi = np.linspace(-0.2, 1.2, 500)

    # Constants for thermal energy scale (arbitrary energy units scaled to kT)
    # At room temperature, kT is roughly 0.6 kcal/mol
    kT = 1.0

    
    # 2. DEFINE THE DARK FREE ENERGY SURFACE (A_dark)
  
    # In the dark, the inactive state (xi=0) is a deep global minimum.
    # The active state (xi=1) is an unstable, high-energy state.
    A_dark = 15.0 * (xi**2) + 5.0 * xi

    
    # 3. DEFINE THE RETINAL MECHANICAL PUSH WORK (W_push)
  
    # When light hits, the straight all-trans retinal introduces severe
    # steric strain at xi=0, but fits perfectly at xi=1.
    # We model this as a decreasing work function as the protein opens.
    W_push = 22.0 * (1.0 - xi)
    # Ensure work doesn't artificially go negative past xi=1 for simulation stability
    W_push[xi > 1.0] = 0.0

    # The light-perturbed free energy surface
    A_light = A_dark + W_push

    -
    # 4. COMPUTE BOLTZMANN PROBABILITIES & PARTITION FUNCTIONS
   
    # Unnormalized probabilities (Boltzmann factors: e^(-A/kT))
    boltz_factors_dark = np.exp(-A_dark / kT)
    boltz_factors_light = np.exp(-A_light / kT)

    # Calculate Partition Functions (Z) via numerical integration (trapezoidal rule)
    d_xi = xi[1] - xi[0]
    Z_dark = np.trapz(boltz_factors_dark, dx=d_xi)
    Z_light = np.trapz(boltz_factors_light, dx=d_xi)

    # Normalized Probability Distributions (P(xi))
    P_dark = boltz_factors_dark / Z_dark
    P_light = boltz_factors_light / Z_light

    
    # 5. CALCULATE EQUILIBRIUM CONSTANTS (K_eq)
    
    # Extract values near the macrostates to find the population ratio
    idx_inactive = np.abs(xi - 0.0).argmin()
    idx_active = np.abs(xi - 1.0).argmin()

    K_eq_dark = P_dark[idx_active] / P_dark[idx_inactive]
    K_eq_light = P_light[idx_active] / P_light[idx_inactive]

   
    # 6. PLOT THE RESULTS
   
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Free Energy Surfaces
    ax1.plot(xi, A_dark, 'b-', label='Dark State ($A_{dark}$)', linewidth=2)
    ax1.plot(xi, A_light, 'r-', label='Light Perturbed State ($A_{light}$)', linewidth=2)
    ax1.set_title('Free Energy Landscape ($A(\\xi)$)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Reaction Coordinate ($\\xi$)\n[0 = Inactive, 1 = Active]', fontsize=10)
    ax1.set_ylabel('Free Energy ($A$ / $k_BT$)', fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend()

    # Plot 2: Boltzmann Probability Distributions
    ax2.plot(xi, P_dark, 'b--', label='Dark Population ($P_{dark}$)', linewidth=2)
    ax2.plot(xi, P_light, 'r--', label='Light Population ($P_{light}$)', linewidth=2)
    ax2.set_title('Boltzmann Population Distribution ($P(\\xi)$)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Reaction Coordinate ($\\xi$)\n[0 = Inactive, 1 = Active]', fontsize=10)
    ax2.set_ylabel('Probability Density', fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()

    plt.tight_layout()
    plt.show()

   
    # 7. PRINT THERMODYNAMIC SUMMARY
    
    print("=== STATISTICAL MECHANICS SIMULATION RESULTS ===")
    print(f"Dark Partition Function (Z_dark):   {Z_dark:.4f}")
    print(f"Light Partition Function (Z_light): {Z_light:.4f}")
    print(f"Dark Equilibrium Constant (K_eq):   {K_eq_dark:.4e} (Pop. heavily favors Inactive)")
    print(f"Light Equilibrium Constant (K_eq):  {K_eq_light:.4e} (Pop. heavily favors Active)")
    print("=================================================")

if __name__ == "__main__":
    simulate_thermodynamic_shift()
