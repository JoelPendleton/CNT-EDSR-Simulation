import numpy as np
import matplotlib.pyplot as plt

f = np.load("./Simulation-13-09-21/simulation.out/B_eff000000.npy")
# print(f.shape)


# Effective magnetic field vectors across the wire
B_eff_vector = np.square(f[0,14:19,:,49]) + np.square(f[1,14:19,:,49]) + np.square(f[2,14,:,49]) # The 3 values represent the vector directions


# Effective magnetic field magnitude across the wire
B_eff_mag= np.sqrt(B_eff_vector)

# Effective magnetic field mean magnitude at each y-coord
B_eff_mag_mean = np.mean(B_eff_mag, axis=0)

# Effective magnetic field magnitude s.d. at each y-coord
B_eff_mag_sd = np.std(B_eff_mag, axis=0)

x_coords = list(range(1, 783))

plt.figure()
plt.plot(x_coords, B_eff_mag_mean, color="black")
plt.ylabel("$B_{eff}$ $(mT)$")
plt.xlabel("$y$ $(10^{-8}m)$")

plt.savefig("../figures/B-eff-nanotube.svg")
plt.errorbar(x_coords, B_eff_mag_mean, yerr=B_eff_mag_sd, label='both limits (default)', color="black")
plt.savefig("../figures/B-eff-nanotube-with-errors.svg")

print("Plots saved!")

