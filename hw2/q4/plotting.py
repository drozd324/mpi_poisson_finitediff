import numpy as np
import matplotlib.pyplot as plt

grid = np.loadtxt("./grids/grid.txt")

fig1, ax1 = plt.subplots()
heatmap = ax1.imshow(grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')
plt.savefig("./plots/heatmap_plot.png", dpi=300)

plt.clf()
# analytic soution
analytic_grid = np.loadtxt("./grids/analytic_grid.txt")
fig2, ax2 = plt.subplots()
heatmap = ax2.imshow(analytic_grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')
plt.savefig("./plots/analytic_heatmap_plot.png", dpi=300)
