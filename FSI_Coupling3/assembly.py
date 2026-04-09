import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# ---- True function ----
def f(x):
    return np.sin(x) + 0.2 * np.cos(3*x)

# Training points (progressively added)
X_all = np.atleast_2d([1, 3, 5, 7, 9]).T
y_all = f(X_all).ravel()

# Test points
X = np.atleast_2d(np.linspace(0, 10, 200)).T

# Kernel for GP
kernel = C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0)
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, alpha=0.05)

# ---- Figure setup ----
fig, ax = plt.subplots(figsize=(7,4))
line_true, = ax.plot(X, f(X), 'r--', lw=1.5, label="True function $f(x)$")
line_mean, = ax.plot([], [], 'b-', lw=2, label="GP mean")
band = ax.fill_between(X.ravel(), 0, 0, color="blue", alpha=0.2)
points = ax.scatter([], [], c="k", s=50, zorder=10, label="Observations")
ax.set_ylim(-3, 3)
ax.set_xlim(0, 10)
ax.set_xlabel("$x$")
ax.set_ylabel("$f(x)$")
ax.legend(loc="upper right")

# ---- Animation functions ----
def init():
    line_mean.set_data([], [])
    points.set_offsets(np.empty((0, 2)))  # Fix: empty 2D array for scatter offsets
    return line_mean, points, band

def animate(i):
    Xi, yi = X_all[:i+1], y_all[:i+1]
    gp.fit(Xi, yi)
    y_pred, sigma = gp.predict(X, return_std=True)

    # Update GP mean
    line_mean.set_data(X, y_pred)

    # Update uncertainty band
    ax.collections.clear()  # clear old bands
    ax.fill_between(X.ravel(),
                    y_pred - 1.96*sigma,
                    y_pred + 1.96*sigma,
                    color="blue", alpha=0.2)

    # Update points
    points.set_offsets(np.c_[Xi.ravel(), yi])

    return line_mean, points

ani = animation.FuncAnimation(fig, animate, frames=len(X_all),
                              init_func=init, blit=False, repeat=False)

# ---- Save video ----
ani.save("gp_explainer.mp4", writer="ffmpeg", fps=1)
plt.close()
