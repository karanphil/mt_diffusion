import numpy as np
import matplotlib.pyplot as plt
from modules.io import plot_init

plot_init()

theta = np.arange(0,360,1)
theta_rad = np.deg2rad(theta)

theta_a = np.arange(0,90,1)
theta_a_rad = np.deg2rad(theta_a)

theta_n = np.zeros((90,360))
for i in range(90):
    theta_n[i,:] = np.arccos(np.sin(theta_a_rad[i])*np.cos(theta_rad))

theta_n_deg = np.rad2deg(theta_n)

top_line = np.max(theta_n_deg, axis=1)
bot_line = np.min(theta_n_deg, axis=1)

plt.figure()
plt.fill_between(theta_a,bot_line, top_line, color=(0,167/255.,89/255.), alpha=0.75)
plt.show()