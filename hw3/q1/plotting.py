import numpy as np
import matplotlib.pyplot as plt
import os

# plot computed grids
for i in range(3):
	grid = np.loadtxt(f"./grids/grid{i}.txt")
	fig, ax = plt.subplots()
	heatmap = ax.imshow(grid, cmap='coolwarm')
	plt.colorbar(heatmap, label='')
	ax.set_title(rf"Heatmap for grid{i} on $\Omega = [0,1] \times [0,1]$")
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	
	plt.savefig(f"./plots/heatmap_plot{i}.png", dpi=300)
	plt.close(fig)

# plot analytic solution
analytic_grid = np.loadtxt("./grids/analytic_grid.txt")
fig, ax = plt.subplots()
heatmap = ax.imshow(analytic_grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')
ax.set_title(r"Analytic Solution Heatmap on $\Omega = [0,1] \times [0,1]$")
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("x")
ax.set_ylabel("y")
	

plt.savefig("./plots/analytic_heatmap_plot.png", dpi=300)
plt.close(fig)

