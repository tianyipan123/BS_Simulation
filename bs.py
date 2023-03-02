# %%
import yfinance as yf
from datetime import date, timedelta
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# %% Gather Stock Information
ticker = "MSFT"
end_date = date(2023, 2, 22)
start_date = end_date - timedelta(days=7)
df = yf.download(ticker, start_date, end_date, interval="1m")["Adj Close"]
print("average = " + str(df.sum() / df.size))
volatility = (df.pct_change().std() * np.sqrt(252 * 390)) ** 2
print("std = " + str(volatility))

# plt.figure()
# plt.plot(df)
# plt.show()

# %% Explicit Method
# setting parameters
T = 1
S_max = 300
S_min = 0
m = 100
n = 33
sigma2 = volatility
r = 0.04
K = 250

# initialize targets
t = np.linspace(0, T, m + 1)
S = np.linspace(S_min, S_max, n + 1)
dt = T / m
dS = S_max / n
V = np.zeros((t.size, S.size))

# set boundary conditions
# value of option when S=S_min is already 0, so no need to set again
V[-1, :] = np.maximum(S.T - K, 0)
V[:, -1] = S_max - K * np.exp(-r * (T - t))

# apply formula
for i in range(m, 0, -1):
    for j in range(n - 1, 0, -1):
        d2V_dS2 = (V[i, j + 1] + V[i, j - 1] - 2 * V[i, j]) / (dS ** 2)
        dV_dS = (V[i, j + 1] - V[i, j - 1]) / (2 * dS)
        V[i - 1, j] = V[i, j] + dt * (0.5 * sigma2 * (S[j] ** 2) * d2V_dS2
                                      + r * S[j] * dV_dS - r * V[i, j])

# plot approximation
t_mesh, S_mesh = np.meshgrid(t, S)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_approx = ax.plot_surface(t_mesh, S_mesh, V.T, linewidth=0,
                              antialiased=False, label="approximation")
surf_approx._edgecolors2d = surf_approx._edgecolor3d
surf_approx._facecolors2d = surf_approx._facecolor3d
ax.set_xlabel("t")
ax.set_ylabel("S")
ax.set_zlabel("V")

# plot true result
t_true, S_true = np.meshgrid(t[:-1], S[1:])  # slicing to avoid division by zero
d1 = (np.log(S_true / K) + (r + 0.5 * sigma2) * (T - t_true)) / (
    np.sqrt(sigma2 * (T - t_true)))
d2 = d1 - np.sqrt(sigma2 * (T - t_true))
V_true = S_true * norm.cdf(d1) - K * np.exp(-r * (T - t_true)) * norm.cdf(d2)
surf_true = ax.plot_surface(t_true, S_true, V_true, linewidth=0,
                            antialiased=False, label="true")
surf_true._edgecolors2d = surf_true._edgecolor3d
surf_true._facecolors2d = surf_true._facecolor3d

plt.title("Black-Scholes Equation, Explicit Method")
ax.legend()
plt.show()

# plot error histogram
plt.figure()
plt.hist((V_true - V.T[:-1, 1:]).reshape(-1), bins=30)
plt.title("Explicit Method, Error Estimate")
plt.show()

# %% Implicit Method
# setting parameters
T = 1
S_max = 300
S_min = 0
m = 100
n = 100
sigma2 = volatility
r = 0.04
K = 250

# initialize targets
t = np.linspace(0, T, m + 1)
S = np.linspace(S_min, S_max, n + 1)
dt = T / m
dS = S_max / n
V = np.zeros((t.size, S.size))

# set boundary conditions
# value of option when S=S_min is already 0, so no need to set again
V[-1, :] = np.maximum(S.T - K, 0)
V[:, -1] = S_max - K * np.exp(-r * (T - t))

# apply formula
S_eye = np.diag(S[1:-1])
d2V_mat = np.eye(n - 1, n - 1, k=-1) + np.eye(n - 1, n - 1, k=1) \
          - 2 * np.eye(n - 1)
dV_mat = np.eye(n - 1, n - 1, k=1) - np.eye(n - 1, n - 1, k=-1)
A = r * np.eye(n - 1) - 0.5 * sigma2 * (S_eye * S_eye / dS / dS) @ d2V_mat \
    - r * S_eye @ dV_mat / (2 * dS)

inv = np.linalg.inv(np.eye(n - 1) + A * dt)
print("cond(A) = " + str(np.linalg.cond(inv)))

for i in range(m, 0, -1):
    # compute f
    f = np.zeros(n - 1)
    # f[0] = 0 since S(t = 0) = 0 and V(S = 0) = 0
    f[-1] = -(0.5 * sigma2 * (S[-1] ** 2) / (dS ** 2) + r * S[-1] / (2 * dS)) \
        * V[i - 1, -1]
    f = f.reshape(-1, 1)

    # compute V
    value = (inv @ (V[i, 1:-1].reshape(-1, 1) - f * dt)).flatten()
    # value = np.linalg.solve(np.eye(n - 1) + A * dt,
    #                         V[i, 1:-1].reshape(-1, 1) - f * dt).flatten()
    V[i - 1, 1:-1] = value

# plot approximation
t_mesh, S_mesh = np.meshgrid(t, S)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_approx = ax.plot_surface(t_mesh, S_mesh, V.T, linewidth=0,
                              antialiased=False, label="approximation")
surf_approx._edgecolors2d = surf_approx._edgecolor3d
surf_approx._facecolors2d = surf_approx._facecolor3d
ax.set_xlabel("t")
ax.set_ylabel("S")
ax.set_zlabel("V")

# plot true result
t_true, S_true = np.meshgrid(t[:-1], S[1:]) # slicing to avoid division by zero
d1 = (np.log(S_true / K) + (r + 0.5 * sigma2) * (T - t_true)) / (np.sqrt(sigma2 * (T - t_true)))
d2 = d1 - np.sqrt(sigma2 * (T - t_true))
V_true = S_true * norm.cdf(d1) - K * np.exp(-r * (T - t_true)) * norm.cdf(d2)
surf_true = ax.plot_surface(t_true, S_true, V_true, linewidth=0, antialiased=False, label="true")
surf_true._edgecolors2d = surf_true._edgecolor3d
surf_true._facecolors2d = surf_true._facecolor3d

plt.title("Black-Scholes Equation, Implicit Method")
ax.legend()
plt.show()

# plot error histogram
plt.figure()
plt.hist((V_true - V.T[:-1, 1:]).reshape(-1), bins=30)
plt.title(f"Implicit Method, Error Estimate, m = {m}, n = {n}")
plt.show()
