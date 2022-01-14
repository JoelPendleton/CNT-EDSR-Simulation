import numpy as np
import matplotlib.pyplot as plt

x_coords = list(range(1, 783))

##### B_0 Field Parallel to Magnets ####

f = np.load("./Simulation-09-01-22/simulation.out/B_eff000000.npy")

# Effective magnetic field vectors across the wire
B_x = f[0,14:19,:,49:54]
B_y = f[1,14:19,:,49:54]
B_z = f[2,14:19,:,49:54]

# Effective magnetic field mean magnitude at each y-coord
B_x_mean = np.mean(B_x, axis=(0,2))
B_y_mean = np.mean(B_y, axis=(0,2))
B_z_mean = np.mean(B_z, axis=(0,2))

# Effective magnetic field magnitude s.d. at each y-coord
B_x_sd = np.std(B_x, axis=(0,2))
B_y_sd = np.std(B_y, axis=(0,2))
B_z_sd = np.std(B_z, axis=(0,2))

fig_parallel = plt.figure()
plt.axvline(x=429, linestyle='--', label="bottom magnet", color="red")
plt.axvline(x=353, linestyle='--', label="top magnets", color="blue") # I need to test whether this is the case, perform simulation
                                                        # one fewer magnet on top. If the double feature disappears
                                                        # then the first line is the top magnets

plt.plot(x_coords, B_x_mean, label="$B_x$")
# plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
plt.plot(x_coords, B_y_mean, label="$B_y$")
# plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
plt.plot(x_coords, B_z_mean, label="$B_z$")
# plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
# plt.title("Effective Magnetic Field across Nanotube for Field Parallel to x-axis")
plt.ylabel("$B_{eff}$ $(mT)$")
plt.xlabel("$y$ $(10^{-8}m)$")
plt.legend()
plt.savefig("../figures/B-eff-nanotube-parallel-to-mag.eps")
plt.close(fig_parallel)



##### B_0 Field Perpendicular to Magnets ####

f = np.load("./Simulation-13-12-21/simulation.out/B_eff000000.npy")


# Effective magnetic field vectors across the wire
B_x = f[0,14:19,:,49:54]
B_y = f[1,14:19,:,49:54]
B_z = f[2,14:19,:,49:54]

# Effective magnetic field mean magnitude at each y-coord
B_x_mean = np.mean(B_x, axis=(0,2))
B_y_mean = np.mean(B_y, axis=(0,2))
B_z_mean = np.mean(B_z, axis=(0,2))

# Effective magnetic field magnitude s.d. at each y-coord
B_x_sd = np.std(B_x, axis=(0,2))
B_y_sd = np.std(B_y, axis=(0,2))
B_z_sd = np.std(B_z, axis=(0,2))

fig_parallel = plt.figure()
plt.axvline(x=429, linestyle='--', label="bottom magnet", color="red")
plt.axvline(x=353, linestyle='--', label="top magnets", color="blue") # I need to test whether this is the case, perform simulation
                                                        # one fewer magnet on top. If the double feature disappears
                                                        # then the first line is the top magnets
plt.plot(x_coords, B_x_mean, label="$B_x$")
# plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
plt.plot(x_coords, B_y_mean, label="$B_y$")
# plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
plt.plot(x_coords, B_z_mean, label="$B_z$")
# plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
# plt.title("Effective Magnetic Field across Nanotube for Field Perpendicular to x-axis")
plt.ylabel("$B_{eff}$ $(mT)$")
plt.xlabel("$y$ $(10^{-8}m)$")
plt.legend()
plt.savefig("../figures/B-eff-nanotube-perpendicular-to-mag.eps")
plt.close(fig_parallel)



print("Plots saved!")

