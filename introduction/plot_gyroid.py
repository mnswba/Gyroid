"""Render the gyroid nodal surface and a sheet-gyroid solid to PNG.

Run:  python introduction/plot_gyroid.py
Out:  introduction/figures/gyroid.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from skimage.measure import marching_cubes

N, CELLS = 48, 2                       # grid points, unit cells per axis
g = np.linspace(0, 2 * np.pi * CELLS, N)
x, y, z = np.meshgrid(g, g, g, indexing="ij")
F = np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)
step = g[1] - g[0]
ls = LightSource(azdeg=315, altdeg=45)


def draw(ax, level, color, title):
    v, f, _, _ = marching_cubes(F, level=level, spacing=(step,) * 3)
    ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2], triangles=f, color=color,
                    lightsource=ls, linewidth=0, antialiased=False, shade=True)
    ax.set_title(title, fontsize=11)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=35)
    ax.set_axis_off()


fig = plt.figure(figsize=(15, 5.2), dpi=130)
ax1 = fig.add_subplot(131, projection="3d")
draw(ax1, 0.0, "#4C8BF5", "Gyroid surface  F = 0\n(2 x 2 x 2 unit cells)")

ax2 = fig.add_subplot(132, projection="3d")
draw(ax2, 0.3, "#E4572E", "")
draw(ax2, -0.3, "#2E86AB",
     "Sheet gyroid  |F| <= 0.3\n(red face: hot side, blue face: cold side)")

ax3 = fig.add_subplot(133)
s = np.linspace(0, 2 * np.pi * CELLS, 400)
X, Y = np.meshgrid(s, s)
sec = np.sin(X) * np.cos(Y) + np.sin(Y)          # F at z = 0
ax3.contourf(X, Y, sec, levels=[-3, -0.3], colors=["#2E86AB"])
ax3.contourf(X, Y, sec, levels=[0.3, 3], colors=["#E4572E"])
ax3.contourf(X, Y, sec, levels=[-0.3, 0.3], colors=["#444444"])
ax3.set_aspect("equal")
ax3.set_title("Slice at z = 0: two channels (red / blue)\n"
              "separated by the wall (grey)", fontsize=11)
ax3.set_xlabel("x (rad)")
ax3.set_ylabel("y (rad)")

plt.tight_layout(rect=(0, 0.02, 1, 0.93))
plt.savefig("introduction/figures/gyroid.png")
print("saved")
