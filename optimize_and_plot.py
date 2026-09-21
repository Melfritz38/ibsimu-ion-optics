import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import minimize

iteration_count = 0

def zielfunktion(v):
    global iteration_count
    iteration_count += 1
    v_l1, v_l2, v_l3, v_l4 = v
    
    cmd = ["./simulation_opt", f"{v_l1:.2f}", f"{v_l2:.2f}", f"{v_l3:.2f}", f"{v_l4:.2f}"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    try:
        rms = float(res.stdout.strip())
    except ValueError:
        rms = 999.0

    print(f"Iter {iteration_count:03d} | L1: {v_l1:7.1f}V | L2: {v_l2:7.1f}V | L3: {v_l3:7.1f}V | L4: {v_l4:7.1f}V ==> RMS: {rms:.4f} mm")
    return rms

# 1. OPTIMIERUNG
start_voltages = [-500.0, -1200.0, -10.0, -500.0]
print("Starte 4D Nelder-Mead Optimierung...")
res = minimize(zielfunktion, start_voltages, method='Nelder-Mead', options={'xatol': 1.0, 'fatol': 0.01, 'maxiter': 200})

opt_l1, opt_l2, opt_l3, opt_l4 = res.x
traj_file = "trajectories_best.dat"

# Optimum erneut simulieren (speichert trajectories_best.dat und epot.dat)
subprocess.run(["./simulation_opt", f"{opt_l1:.2f}", f"{opt_l2:.2f}", f"{opt_l3:.2f}", f"{opt_l4:.2f}", traj_file])

# 2. PLOT ERSTELLEN
fig, ax = plt.subplots(figsize=(13, 6))

# A. Potenzialfeld & Äquipotenziallinien zeichnen
try:
    pot_data = np.loadtxt('epot.dat')
    x_vals = np.unique(pot_data[:, 0])
    y_vals = np.unique(pot_data[:, 1])
    z_vals = pot_data[:, 2].reshape(len(x_vals), len(y_vals)).T

    X, Y = np.meshgrid(x_vals, y_vals)
    contours = ax.contour(X, Y, z_vals, levels=30, cmap='coolwarm', alpha=0.45, linewidths=0.7)
    ax.clabel(contours, inline=True, fontsize=6, fmt='%1.0f V')
except Exception as e:
    print("Hinweis: Potenzialfeld konnte nicht geladen werden:", e)

# B. Linsenelektroden einzeichnen
lenses = [
    (0.013, 0.018, 0.0050, 0.035, 'Repeller (30V)'),
    (0.018, 0.020, 0.0075, 0.035, f'L1 ({opt_l1:.0f}V)'),
    (0.030, 0.094, 0.0125, 0.035, f'L2 ({opt_l2:.0f}V)'),
    (0.096, 0.099, 0.0075, 0.035, f'L3 ({opt_l3:.0f}V)'),
    (0.101, 0.341, 0.0130, 0.035, f'L4 ({opt_l4:.0f}V)'),
    (0.361, 0.363, 0.0150, 0.035, 'L5 (0V)')
]

for x1, x2, r_in, r_out, label in lenses:
    width = x2 - x1
    height = r_out - r_in
    rect = patches.Rectangle((x1, r_in), width, height, linewidth=1, edgecolor='black', facecolor='darkgray', alpha=0.8, zorder=3)
    ax.add_patch(rect)
    ax.text(x1 + width/2, r_out + 0.001, label, rotation=35, ha='left', va='bottom', fontsize=8, zorder=4)

# C. Probenoberfläche (x = 0.380 m)
x_probe = 0.380
ax.axvline(x=x_probe, color='darkblue', linestyle='--', linewidth=2, label='Probenoberfläche', zorder=3)
probe_rect = patches.Rectangle((x_probe, 0.0), 0.005, 0.035, facecolor='blue', alpha=0.2, zorder=2)
ax.add_patch(probe_rect)

# D. Trajektorien plotten
with open(traj_file) as f:
    lines = f.readlines()

x_curr, y_curr = [], []
for line in lines:
    if line.strip() == '':
        if x_curr:
            ax.plot(x_curr, y_curr, color='red', alpha=0.5, linewidth=0.9, zorder=4)
            x_curr, y_curr = [], []
    else:
        parts = line.split()
        x_curr.append(float(parts[0]))
        y_curr.append(float(parts[1]))

ax.axhline(0, color='black', linestyle=':', linewidth=1, label='Optische Achse')
ax.set_xlim(0.0, 0.400)
ax.set_ylim(0.0, 0.045)
ax.set_xlabel('Z-Position / Flugrichtung x (m)')
ax.set_ylabel('Strahlradius r (m)')
ax.set_title(f'Optimierter Ionenstrahlverlauf mit Feldlinien (RMS: {res.fun:.3f} mm)')
ax.grid(True, linestyle='--', alpha=0.3)
ax.legend(loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig('ion_optics_fields.png', dpi=300)
plt.show()

print("Fertiges Diagramm als 'ion_optics_fields.png' gespeichert.")
