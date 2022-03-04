import numpy as np
import matplotlib.pyplot as plt

x_coords = list(range(0, 182))



file_name = "simulation-0-y"
# angle = 5
# def return_indices(y):
#     y_index = y * np.cos(np.radians(angle))
#     x_index = y * np.sin(np.radians(angle))
#
#     return  x_index, y_index
# naming convention = simulation-[degrees CNT rotated clockwise]-[direction of external field]

f = np.load("./Simulation-16-02-22/{}.out/B_eff000000.npy".format(file_name))

y_indices = np.arange(8,190,1)
# print(return_indices(y_indices))
# Effective magnetic field vectors across the wire

# For each y need to add on some x-shift


B_x = f[0,14:19,8:190,49:54]
B_y = f[1,14:19,8:190,49:54]
B_z = f[2,14:19,8:190,49:54]

# Effective magnetic field mean magnitude at each y-coord
B_x_mean = np.mean(B_x, axis=(0,2))
B_y_mean = np.mean(B_y, axis=(0,2))
B_z_mean = np.mean(B_z, axis=(0,2))

# Effective magnetic field magnitude s.d. at each y-coord
B_x_sd = np.std(B_x, axis=(0,2))
B_y_sd = np.std(B_y, axis=(0,2))
B_z_sd = np.std(B_z, axis=(0,2))

fig = plt.figure()
plt.axvline(x=52, linestyle='--', label="middle of bottom magnet")
plt.axvline(x=128, linestyle='--', label="middle of top magnets")

plt.axvline(x=66, linestyle='--', label="End of Dot 1", color="grey")
plt.axvline(x=116, linestyle='--', label="Start of Dot 2", color="grey")

plt.plot(x_coords, B_x_mean, label="$B_x$")
# plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
plt.plot(x_coords, B_y_mean, label="$B_y$")
# plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
plt.plot(x_coords, B_z_mean, label="$B_z$")
# plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
plt.ylabel("$B_{eff}$ $(T)$")
plt.xlabel("$y$ $(10^{-8}m)$")
plt.legend()
plt.savefig("../figures/{}.svg".format(file_name))
plt.close(fig)


print("Plot saved!")

