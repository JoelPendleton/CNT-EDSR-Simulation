def plot_spectrum(syst, Bfields, length):
#
#     def potential(x): # Infinite square well
#         if 0 < x < length:
#             return 0
#         else:
#             return 99999999999
#
#     energies = []
#     for b in Bfields:
#         params = dict(A=0, B=0, C=C_constant, V=potential)
#
#         # Obtain the Hamiltonian as a sparse matrix
#         ham_mat = syst.hamiltonian_submatrix(sparse=True, params=params)
#         k = 4
#         evals, evecs = sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=k * 2, sigma=0))
#
#         # print(vectors) # It's giving 2x the number of eigenvectors.
#
#         # The jumping is because of the zeeman splitting and how it only calculates the lowest two levels
#         evals = evals[::2]
#
#         energies.append(evals)
#
#
#     plt.figure()
#     plt.plot(Bfields, energies, '.k')
#
#
#     plt.xlabel("magnetic field ()")
#     plt.ylabel("energy [H]")
#     plt.savefig("spectrum.png") # With A = 0 we expect straight forward zeeman splitting