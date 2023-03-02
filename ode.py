#%%
import numpy as np
import matplotlib.pyplot as plt

#%% Logistic Equation: Euler's Method
hs = [0.2, 0.1, 0.05]
a = 0.1
plt.figure()
# Euler's method
for h in hs:
    value = [a]
    steps = np.arange(0, 1, h)
    for _ in steps:
        x = value[-1]
        value.append(x + x * (1 - x) * h)
    plt.plot(np.append(steps, 1), value)
# True solution
steps = np.linspace(0, 1, 100)
x = a * np.exp(steps) / (1 - a + a * np.exp(steps))
plt.plot(steps, x)
legends = ["h = " + str(h) for h in hs]
legends.append("true solution")
plt.legend(legends)
plt.title("Euler's Method")
plt.xlabel("t")
plt.ylabel("x")
plt.show()

#%% Stiff Problem: Euler's method
hs = [0.2, 0.1, 0.05]
a = 1
lam = -10
plt.figure()
# Euler's method
for h in hs:
    value = [a]
    steps = np.arange(0, 1, h)
    for _ in steps:
        x = value[-1]
        value.append(x + lam * x * h)
    plt.plot(np.append(steps, 1), value)
# True solution
steps = np.linspace(0, 1, 100)
x = a * np.exp(lam * steps)
plt.plot(steps, x)
legends = ["h = " + str(h) for h in hs]
legends.append("true solution")
plt.legend(legends)
plt.title("Euler's Method, Overshoot, $\lambda=-10$")
plt.xlabel("t")
plt.ylabel("x")
plt.show()

#%% Stiff Problem: Implicit Method
hs = [0.2, 0.1, 0.05]
a = 1
lam = -10
plt.figure()
# Euler's method
for h in hs:
    steps = np.arange(0, 1, h)
    mult = 1 / (1 - h * lam)
    value = a * np.power(mult, np.arange(steps.size + 1))
    plt.plot(np.append(steps, 1), value)
# True solution
steps = np.linspace(0, 1, 100)
x = a * np.exp(lam * steps)
plt.plot(steps, x)
legends = ["h = " + str(h) for h in hs]
legends.append("true solution")
plt.legend(legends)
plt.title("Implicit Method, $\lambda=-10$")
plt.xlabel("t")
plt.ylabel("x")
plt.show()
