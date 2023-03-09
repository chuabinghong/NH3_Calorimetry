from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,10))
ax = Axes3D(fig)
colors = 'b', 'r'
labels = 'Group1', 'Group2'

for i, c, label in zip(range(len(labels)), colors, labels):
    ax.scatter(tsne_fit[data['Group']==i, 0], tsne_fit[data['Group']==i, 1], tsne_fit[data['Group']==i, 2], s=30, c=c, label=label, alpha=0.5)
fig.legend()