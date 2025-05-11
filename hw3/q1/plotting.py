import numpy as np
import matplotlib.pyplot as plt
import os

# Ensure output directory exists
os.makedirs("./plots", exist_ok=True)

# Plot computed grids
for i in range(3):
    grid = np.loadtxt(f"./grids/grid{i}.txt")
    fig, ax = plt.subplots()
    heatmap = ax.imshow(grid, cmap='coolwarm')
    plt.colorbar(heatmap, label='')
    ax.set_title(f"Heatmap for grid{i}")
    plt.savefig(f"./plots/heatmap_plot{i}.png", dpi=300)
    plt.close(fig)

# Plot analytic solution
analytic_grid = np.loadtxt("./grids/analytic_grid.txt")
fig, ax = plt.subplots()
heatmap = ax.imshow(analytic_grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')
ax.set_title("Analytic Solution Heatmap")
plt.savefig("./plots/analytic_heatmap_plot.png", dpi=300)
plt.close(fig)

