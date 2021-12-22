import numpy as np
import matplotlib.pyplot as plt

f = np.load("./Simulation-13-09-21/simulation.out/B_eff000000.npy")
# print(f.shape)


# Effective magnetic field vectors across the wire
B_x = f[0,14:19,:,49]
B_y = f[1,14:19,:,49]
B_z = f[2,14:19,:,49]
B_eff_vector = np.square(B_x) + np.square(B_y) + np.square(B_z) # The 3 values represent the vector directions


# Effective magnetic field magnitude across the wire
B_eff_mag= np.sqrt(B_eff_vector)

# Effective magnetic field mean magnitude at each y-coord
B_eff_mag_mean = np.mean(B_eff_mag, axis=0)
B_x_mean = np.mean(B_x, axis=0)
B_y_mean = np.mean(B_y, axis=0)
B_z_mean = np.mean(B_z, axis=0)


# Effective magnetic field magnitude s.d. at each y-coord
B_eff_mag_sd = np.std(B_eff_mag, axis=0)
B_x_sd = np.std(B_x, axis=0)
B_y_sd = np.std(B_y, axis=0)
B_z_sd = np.std(B_z, axis=0)

x_coords = list(range(1, 783))

plt.figure()
plt.plot(x_coords, B_eff_mag_mean, color="black")
plt.ylabel("$B_{eff}$ $(T)$")
plt.xlabel("$y$ $(10{-8}m)$")

plt.savefig("../figures/B-eff-magnitude-nanotube.svg")
plt.errorbar(x_coords, B_eff_mag_mean, yerr=B_eff_mag_sd, label='both limits (default)', color="black")
plt.savefig("../figures/B-eff-magnitude-nanotube-with-errors.svg")

plt.figure()
plt.axvline(x=429, linestyle='--', label="bottom magnet", color="red")
plt.axvline(x=353, linestyle='--', label="top magnets", color="blue") # I need to test whether this is the case, perform simulation
                                                        # one fewer magnet on top. If the double feature disappears
                                                        # then the first line is the top magnets

plt.plot(x_coords, B_x_mean, label="$B_x$")
plt.plot(x_coords, B_y_mean, label="$B_y$")
plt.plot(x_coords, B_z_mean, label="$B_z$")
plt.ylabel("$B_{eff}$ $(T)$")
plt.xlabel("$y$ $(10{-8}m)$")
plt.legend()
plt.savefig("../figures/B-eff-nanotube.svg")


print("Plots saved!")

