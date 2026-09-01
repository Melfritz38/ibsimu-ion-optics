import subprocess
import numpy as np
from scipy.optimize import minimize

iteration_count = 0

def zielfunktion(v):
    global iteration_count
    iteration_count += 1
    
    v_l1, v_l2, v_l3, v_l4 = v
    
    # Ausführen des C++ Simulators
    cmd = [
        "./simulation_opt",
        f"{v_l1:.2f}",
        f"{v_l2:.2f}",
        f"{v_l3:.2f}",
        f"{v_l4:.2f}"
    ]
    
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    try:
        rms = float(res.stdout.strip())
    except ValueError:
        rms = 999.0

    print(f"Iter {iteration_count:03d} | L1: {v_l1:7.1f}V | L2: {v_l2:7.1f}V | L3: {v_l3:7.1f}V | L4: {v_l4:7.1f}V ==> RMS: {rms:.4f} mm")
    return rms

# Startwerte basierend auf deinen bisherigen Messergebnissen
start_voltages = [-500.0, -1200.0, -10.0, -500.0]

print("Starte 4D Nelder-Mead Optimierung...")
res = minimize(
    zielfunktion, 
    start_voltages, 
    method='Nelder-Mead',
    options={'xatol': 1.0, 'fatol': 0.01, 'maxiter': 200}
)

opt_l1, opt_l2, opt_l3, opt_l4 = res.x

print("\n==========================================")
print("OPTIMIERUNG ABGESCHLOSSEN!")
print(f"Beste Spannungen:")
print(f"  V_L1 = {opt_l1:.1f} V")
print(f"  V_L2 = {opt_l2:.1f} V")
print(f"  V_L3 = {opt_l3:.1f} V")
print(f"  V_L4 = {opt_l4:.1f} V")
print(f"Minimaler RMS Strahlradius: {res.fun:.4f} mm")
print("==========================================\n")
