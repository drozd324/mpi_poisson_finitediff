import numpy as np
import matplotlib.pyplot as plt

grid = np.loadtxt("./grids/grid.txt")

#x_labels = ['A', 'B', 'C', 'D', 'E']
#y_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May']

# Plot
fig, ax = plt.subplots()
heatmap = ax.imshow(grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')

# Customize
#ax.set_xticks(np.arange(len(x_labels)))
#ax.set_yticks(np.arange(len(y_labels)))
#ax.set_xticklabels(x_labels)
#ax.set_yticklabels(y_labels)
#plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
#ax.set_title("Annotated Heatmap")
#plt.tight_layout()
plt.savefig("./plots/heatmap_plot.png", dpi=300)


# analytic soution
analytic_grid = np.loadtxt("./grids/analytic_grid.txt")

plt.clf()

fig, ax = plt.subplots()
heatmap = ax.imshow(analytic_grid, cmap='coolwarm')
plt.colorbar(heatmap, label='')
plt.savefig("./plots/analytic_heatmap_plot.png", dpi=300)
 
