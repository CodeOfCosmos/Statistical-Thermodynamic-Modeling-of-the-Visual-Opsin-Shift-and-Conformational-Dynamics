import numpy as np

def calculate_state_free_energies():
    # 1. Setup fine grid for numerical integration
    xi = np.linspace(-0.2, 1.2, 2000)
    d_xi = xi[1] - xi[0]
    kT = 1.0  # Energy units scaled to kT
    
    # 2. Define the landscapes (calibrated constants from previous section)
    A_dark = 15.0 * (xi**2) + 5.0 * xi
    W_push = 22.0 * (1.0 - xi)
    W_push[xi > 1.0] = 0.0
    A_light = A_dark + W_push
    
    # 3. Create boolean masks to separate Cis and Trans zones
    cis_mask = xi <= 0.5
    trans_mask = xi > 0.5
    
    # 4. DARK STATE EVALUATION
    # Boltzmann factors in the dark
    bf_dark_cis = np.exp(-A_dark[cis_mask] / kT)
    bf_dark_trans = np.exp(-A_dark[trans_mask] / kT)
    
    # Integrate to find Z_cis and Z_trans
    Z_cis_dark = np.trapz(bf_dark_cis, dx=d_xi)
    Z_trans_dark = np.trapz(bf_dark_trans, dx=d_xi)
    
    # Calculate Helmholtz Free Energy for each pool
    F_cis_dark = -kT * np.log(Z_cis_dark)
    F_trans_dark = -kT * np.log(Z_trans_dark)
    delta_F_dark = F_trans_dark - F_cis_dark
    
    # 5. LIGHT STATE EVALUATION
    # Boltzmann factors in the light
    bf_light_cis = np.exp(-A_light[cis_mask] / kT)
    bf_light_trans = np.exp(-A_light[trans_mask] / kT)
    
    # Integrate to find Z_cis and Z_trans
    Z_cis_light = np.trapz(bf_light_cis, dx=d_xi)
    Z_trans_light = np.trapz(bf_light_trans, dx=d_xi)
    
    # Calculate Helmholtz Free Energy for each pool
    F_cis_light = -kT * np.log(Z_cis_light)
    F_trans_light = -kT * np.log(Z_trans_light)
    delta_F_light = F_trans_light - F_cis_light
    
    # 6. Display Results
    print("=== DISCRETE TWO-STATE THERMODYNAMIC MODELLING ===")
    print("\n[1] DARK STATE (Inactive Rhodopsin Baseline):")
    print(f"  Z_cis   = {Z_cis_dark:.6f}  ->  F_cis   = {F_cis_dark:.4f} kT")
    print(f"  Z_trans = {Z_trans_dark:.6f}  ->  F_trans = {F_trans_dark:.4f} kT")
    print(f"  Net Driving Force (Del F) = {delta_F_dark:+.4f} kT")
    print("  Interpretation: Del F is highly POSITIVE. System cannot spontaneously activate.")
    
    print("\n[2] LIGHT STATE (Post-Isomerization Strain):")
    print(f"  Z_cis   = {Z_cis_light:.6f}  ->  F_cis   = {F_cis_light:.4f} kT")
    print(f"  Z_trans = {Z_trans_light:.6f}  ->  F_trans = {F_trans_light:.4f} kT")
    print(f"  Net Driving Force (Del F) = {delta_F_light:+.4f} kT")
    print("  Interpretation: Del F is now NEGATIVE! System spontaneously flows into the Trans state.")
    print("==================================================")

if __name__ == "__main__":
    calculate_state_free_energies()
