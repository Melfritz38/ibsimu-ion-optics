import matplotlib.pyplot as plt
import numpy as np

# 1. Sweep-Kurve plotten
sweep_data = np.loadtxt('sweep_data.dat')
plt.figure(figsize=(8, 4))
plt.plot(sweep_data[:, 0], sweep_data[:, 1], 'b-o', markersize=3, label='Treffer (r <= 2mm)')
plt.xlabel('Spannung V_L3 (V)')
plt.ylabel('Anzahl fokussierter Ionen')
plt.title('Fokussierungskurve in Abhängigkeit von V_L3')
plt.grid(True)
plt.savefig('sweep_plot.png', dpi=300)
plt.close()

# 2. Beste Trajektorien plotten
plt.figure(figsize=(10, 4))
with open('trajectories_best.dat') as f:
    lines = f.readlines()

x, y = [], []
for line in lines:
    if line.strip() == '':
        if x:
            plt.plot(x, y, 'r-', alpha=0.3, linewidth=0.8)
            x, y = [], []
    else:
        parts = line.split()
        x.append(float(parts[0]))
        y.append(float(parts[1]))

plt.xlabel('x (m)')
plt.ylabel('Radius r (m)')
plt.title('Ionenstrahl-Trajektorien beim Optimalwert')
plt.grid(True)
plt.savefig('trajectories_plot.png', dpi=300)
plt.close()

print("Diagramme erfolgreich als 'sweep_plot.png' und 'trajectories_plot.png' gespeichert!")
