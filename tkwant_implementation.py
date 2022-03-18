import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
from numpy import vectorize
from matplotlib import animation
from tkwant import onebody
import tkwant
import os.path
import json
import time as timelib

from csv import writer


def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)


class System:
    def __init__(self, hamiltonian, pertubation_type="sin", magnetic_fields="none", number_of_lattices=50,
                 potential_type=0):

        # Define units
        self.potential_type = potential_type
        self.hamiltonian = hamiltonian
        self.number_of_lattices = number_of_lattices

        # Constants in SI units
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_SI = 0.66e-6
        self.B_0_SI = 250e-3  # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.b_sl_SI = 1.16e6
        self.m_SI = 9.11e-31 / 10  # divide by ten to get the effective mass
        self.mu_B_SI = 9.2740100783e-24
        self.lattice_size_SI = self.total_length_SI / self.number_of_lattices
        self.z_SI = np.arange(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattice_size_SI,
                              dtype=np.double)
        self.magnetic_fields = magnetic_fields
        # Constants in a.u.
        self.B_0_au = self.tesla_to_au(self.B_0_SI)
        self.b_sl_au = self.tesla_to_au(self.b_sl_SI) * (1 / self.m_to_au(1))
        self.hbar_au = 1
        self.g = 2
        self.m_au = 1 / 10
        self.e_au = 1
        self.total_length_au = self.m_to_au(self.total_length_SI)  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # Distance in a.u. between lattice points
        self.mu_B_au = .5
        self.a_0 = 1
        self.z_au = np.arange(-self.total_length_au / 2, self.total_length_au / 2, self.lattice_size_au,
                              dtype=np.double)
        self.pulse_frequency_au = (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au)  # in Hz
        self.pertubation_type = pertubation_type

    def cosine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        return ((eV_0_au * np.cos(2 * np.pi * pulse_frequency_au * time)) * z) / total_length_au

    def sine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        return ((eV_0_au * np.sin(2 * np.pi * pulse_frequency_au * time)) * z) / total_length_au

    def tesla_to_au(self, tesla):
        """
        Function to convert the magnetic flux density from SI units to AU units
        :param tesla: the magnetic flux density in teslas.
        :return: the magnetic flux density in AU.
        """
        return tesla / 2.35e5

    def au_to_tesla(self, au):
        """
        Function to convert the magnetic flux density from SI units to AU units
        :param tesla: the magnetic flux density in teslas.
        :return: the magnetic flux density in AU.
        """
        return au * 2.35e5

    def second_to_au(self, time):
        return time * 4.1341373336493e16

    def au_to_second(self, time):
        return time / 4.1341373336493e16

    def hartree_to_ev(self, hartree):
        return hartree * 2.72114e1

    def ev_to_hartree(self, ev):
        return ev / 2.72114e1

    def au_to_m(self, au):
        return self.a_0_SI * au

    def m_to_au(self, m):
        return m / self.a_0_SI

    def hz_to_au(self, hz):

        return hz * 1.51983e-16

    def au_to_hz(self, au):

        return au / 1.51983e-16

    def import_mumax3_simulations(self):

        file_name = "simulation-0-x-shifted-50nm-right"

        f = np.load("./magnetic-field-simulations/Simulation-10-03-22/{}.out/B_eff000000.npy".format(file_name))
        # y_indices = np.arange(8, 74, 1)
        # print(return_indices(y_indices))
        # Effective magnetic field vectors across the wire

        # For each y need to add on some z-shift
        B_x = f[0, 1, 0:100, 149]
        B_y = f[1, 1, 0:100, 149]
        B_z = f[2, 1, 0:100, 149]

        # Effective magnetic field mean magnitude at each y-coord
        # B_x_mean = np.mean(B_x, axis=(0, 2))
        # B_y_mean = np.mean(B_y, axis=(0, 2))
        # B_z_mean = np.mean(B_z, axis=(0, 2))
        # B_y_mean is actually in the z-direction
        # B_x_mean is actually in the x-direction
        # B_z_mean is actually in the y-direction

        fig = plt.figure()

        plt.plot(self.z_SI, B_x, label="$B_x$")
        # plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
        plt.plot(self.z_SI, B_y, label="$B_y$")
        # plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
        plt.plot(self.z_SI, B_z, label="$B_z$")
        # plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
        plt.ylabel("$B_{eff}$ $(T)$")
        plt.xlabel("$z$ $(m)$")
        plt.legend()
        plt.savefig("./figures/{}.svg".format(file_name))
        plt.close(fig)

        return B_x, B_z, B_y

    def infinite_square_well_potential(self):

        self.total_length_SI = 0.66e-6
        self.E_sl_ev = 1e-6
        self.B_0_SI = 50e-3  # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.B_0_au = self.tesla_to_au(self.B_0_SI)

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (
                self.g * self.mu_B_au * self.total_length_au)  # slanted magnetic field value computed from E_sl
        self.eV_0_au = self.ev_to_hartree(1e-6)  # the value in a.u. of eV_0.

    def parabolic_potential(self):
        self.total_length_SI = 0.66e-6
        self.z_SI = np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.number_of_lattices)
        self.B_0_SI = 5e-3  # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.E_sl_ev = 1e-6

        self.total_length_au = self.total_length_SI / self.a_0_SI  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # The distance in atomic units spanned by a lattice point
        self.z_au = np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.number_of_lattices)

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (
                self.g * self.mu_B_au * self.total_length_au)  # slanted magnetic field value computed from E_sl
        self.omega_0 = self.ev_to_hartree(1e-3)
        self.eV_0_au = self.ev_to_hartree(10e-6)  # the value in a.u. of eV_0.

        print("b_sl is", self.tesla_to_au(self.b_sl_au) / self.au_to_m(1), "T.")
        print("hbar * omega_0 is", self.hartree_to_ev(self.omega_0), "eV.")
        print("eV_0 is", self.hartree_to_ev(self.eV_0_au), "eV.")
        print("E_sl is", self.E_sl_ev, "eV.")

    def potential(self, z, time):  # Infinite square well
        """
        Function to define the potential of the lead.
        :param x: the position in the system.
        :return: the potential energy.
        """

        if self.potential_type == 0:  # if hard-wall potential is used
            total_potential = 0  # define zero for the potential inside the scattering region.
            # outside the potential will be infinity
        elif self.potential_type == 1:  # for a parabolic potential
            total_potential = .5 * self.m_au *(
                    (z * self.omega_0) ** 2)  # define a parabolic potential inside the scattering region.
        if self.pertubation_type == "cos":
            self.pertubation = self.cosine_v_ac
        else:
            self.pertubation = self.sine_v_ac
        total_potential += self.pertubation(time, z, self.eV_0_au, self.pulse_frequency_au, self.total_length_au)

        return total_potential

    def kwant_shape(self, site):
        """
        function to define the shape of the scattering region.
        :param site: the current site.
        :return: the a boolean saying whether the scattering site should be drawn
        """
        (z,) = site.pos
        return (-self.total_length_au / 2 <= z < self.total_length_au / 2)

    def make_system(self):
        """
        Function to create the system
        :param length: the length of the nanotube
        :return: the system object
        """

        if self.potential_type == 1:  # if we want a parabolic potential
            self.potential_text = "parabolic"
            self.parabolic_potential()  # define constants using parabolic potential function
        else:  # if infinite square well
            self.potential_text = "infinite-well"
            self.infinite_square_well_potential()  # define constants using hard-wall potential function

        self.template = kwant.continuum.discretize(self.hamiltonian,
                                                   grid=self.lattice_size_au)  # make template for lattice points based on the inputted hamiltonian/

        self.syst = kwant.Builder()

        # Add the nanotube to the system
        self.syst.fill(self.template, self.kwant_shape, (0,))

        # kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        self.A_constant = self.hbar_au ** 2 / (2 * self.m_au)
        if self.magnetic_fields != "none":
            B_x, B_y, B_z = self.import_mumax3_simulations()

        def B_function(z):
            if self.magnetic_fields != "none":
                index = np.around(z, 3) == np.around(self.z_au, 3)
                return -self.g * self.mu_B_au * self.tesla_to_au(B_z[index]) * self.hbar_au / 2
            else:
                return -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2

        def C_function(z):
            if self.magnetic_fields != "none":
                index = np.around(z, 3) == np.around(self.z_au, 3)
                return -self.g * self.mu_B_au * self.tesla_to_au(B_x[index]) * self.hbar_au / 2
            else:
                return -self.b_sl_au * z

        def D_function(z):
            if self.magnetic_fields != "none":
                index = np.around(z, 3) == np.around(self.z_au, 3)
                return -self.g * self.mu_B_au * self.tesla_to_au(B_y[index]) * self.hbar_au / 2
            else:
                return 0

        self.params = dict(A=self.A_constant, V=self.potential, B=B_function, C=C_function, D=D_function)

        self.tparams = self.params.copy()
        self.tparams['time'] = 0  # the initial time
        print("System intialised.")

        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=self.tparams)
        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        # Sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        self.initial_eigenvalues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        self.psi_1_init = eigenVectors[:, 0]
        self.psi_2_init = eigenVectors[:, 1]

        self.spin_up_state = tkwant.onebody.WaveFunction.from_kwant(syst=self.syst,
                                                                    psi_init=self.psi_1_init,
                                                                    energy=eigenValues[0],
                                                                    params=self.params)

        return self.syst

    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        :param syst: the system object.
        :return: the sorted eigenvalues and eigenvectors.
        """

        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=self.tparams)
        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        # Sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        # # perturb = tkwant.onebody.kernels.PerturbationExtractor(self.syst, time_name='time', time_start=0)#, params=self.params)
        # self.spin_up_state = tkwant.onebody.WaveFunction.from_kwant(syst=self.syst,
        #                                                             psi_init=eigenVectors[:, 0],
        #                                                             energy=eigenValues[0],
        #                                                             params=self.params)

        return eigenValues, eigenVectors

    def initial_pdfs(self):

        """
        Procedure to show the probability density function.
        :param syst: the system object.
        :return:
        """
        eigenValues, eigenVectors = self.eigenstates()

        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.
        psi1 = self.psi_1_init
        psi1_up, psi1_down = psi1[::2], psi1[1::2]
        # even indices give the spin up and odd indices give the spin down states
        density_1 = np.abs(psi1_up) ** 2 + np.abs(psi1_down) ** 2
        psi2 = self.psi_2_init
        psi2_up, psi2_down = psi2[::2], psi2[1::2]
        density_2 = np.abs(psi2_up) ** 2 + np.abs(psi2_down) ** 2

        fig = plt.figure()

        plt.plot(self.z_SI, density_1, label=r'$|1, -\frac{1}{2}\rangle$')
        plt.plot(self.z_SI, density_2, label=r'$|1, +\frac{1}{2}\rangle$')

        plt.xlabel("$x$ (a.u.)")
        plt.ylabel("$|\psi(x)|^2$")
        plt.legend(loc="upper right")

        # save file with name according to potential type.
        plt.savefig("/content/drive/My Drive/initial-pdfs.svg")  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of PDFs at t=0 saved.")

        return density_1, density_2

    def initial_energies(self):
        """
        Procedure to display the potential and energy levels of the system
        :param syst: the system object.
        :return:
        """

        fig = plt.figure()

        energies = self.hartree_to_ev(np.real(self.initial_eigenvalues))
        y = energies[0:100:2]





        # E_1_actual = (1)**2 * (np.pi)**2 * self.hbar_au**2 / (2 * self.m_au * self.total_length_au**2)
        # E_2_actual = (2)**2 * (np.pi)**2 * self.hbar_au**2 / (2 * self.m_au * self.total_length_au**2)
        # E_1_actual = (1/2)*self.hbar_au * self.omega_0
        # E_2_actual = (1+1/2)*self.hbar_au * self.omega_0
        #
        # E_1_sim = y[0]
        # E_2_sim = y[2]
        # zeeman_sim = y[1]-y[0]
        # print("Simulated zeeman splitting is",zeeman_sim)
        #
        # print("Actual zeeman splitting is",self.g * self.mu_B_au * self.B_0_au)
        # print("E_1 actual:", E_1_actual)
        # print("E_2 actual:", E_2_actual)
        # print("E_1 simulated:", E_1_sim)
        # print("E_2 simulated:", E_2_sim)

        # m, b = np.polyfit(x[0::2], y[0::2], 1)
        # a, b, c = np.polyfit(x[::2], y[::2], 2)

        # t = self.hbar_au ** 2 / (2 * self.m_au * (self.lattice_size_au)**2)
        # print("Is the model a good approximation?", energies[2] < t)
        # plt.plot(x[0::2], m * x[0::2] + b, '--')
        # plt.plot(x, a * x**2 + b*x + c, '--')

        # Plot the energies of the different levels.
        # plt.plot([0, 1], [y[0], y[0]], label=r"$G_-$")
        # plt.plot([0, 1], [y[1], y[1]], label=r"$G_+$")
        # plt.plot([0, 1], [y[2], y[2]], label=r"$E_-$")
        # plt.plot([0, 1], [y[3], y[3]], label=r"$E_+$")

        if self.potential_type == 0:
            x = np.array(list(range(1, 100 + 1)))
            a, b, c = np.polyfit(x, y, 2)
            y_best_fit = (x ** 2) * a + x * b + c
            data = {
                "a": a,
                "b": b,
                "c": c
            }
        else:
            x = np.array(list(range(0, 49 + 1)))
            a, b = np.polyfit(x, y, 1)
            y_best_fit = x * a + b
            data = {
                "a": a,
                "b": b
            }
        json_string = json.dumps(data)

        with open('./energies-{0}-{1}.json'.format(self.potential_text, self.number_of_lattices), 'w') as outfile:
            outfile.write(json_string)

        plt.plot(x, y, '.', label="Eigenenergies")
        plt.plot(x,y_best_fit, '--', label="Line of Best Fit")

        plt.legend(loc="best")
        plt.xlabel("$n$")
        plt.ylabel("$E$ (eV)")

        plt.savefig("./energies-{0}-{1}.eps".format(self.potential_text, self.number_of_lattices))  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of eigenenergies at t=0 saved.")
        return True

    def evolve(self, time_steps=100):
        timestr = timelib.strftime("%Y%m%d-%H%M%S")

        def x_onsite(site):  # function to compute the position operator matrix.
            return [site.pos[0]] * np.identity(2)

        # spin matrices
        sigma_x = np.array([[0, 1],
                            [1, 0]])
        sigma_y = np.array([[0, 1j],
                            [-1j, 0]])
        sigma_z = np.array([[1, 0],
                            [0, -1]])

        # extract lowest two energies (our qubit states)
        E_1 = np.real(self.initial_eigenvalues[0])
        E_2 = np.real(self.initial_eigenvalues[1])

        # extract the state vectors corresponding to these lowest eigenenergies.
        psi1 = self.psi_1_init
        psi2 = self.psi_2_init

        # compute the difference in these energies
        self.delta_E = np.abs(E_2 - E_1)

        # compute the resonant frequency for the rabi oscillations
        omega_res = self.delta_E / self.hbar_au
        self.pulse_frequency_au = omega_res / (
                2 * np.pi)  # set the frequency of the V_AC(t) to be equal to the res. freq.

        eff_B_field_au = self.b_sl_au * self.z_au  # compute the slanted field at each point along the CNT.

        # define the density operator of x
        rho_x = kwant.operator.Density(self.syst, x_onsite, sum=True)
        # compute this density operator on the ground states <2|x|1>:
        rho_x_2_1 = rho_x(psi2, psi1)

        # compute the energy E_x
        E_x = np.real(2 * self.eV_0_au * rho_x_2_1 / self.total_length_au)
        self.t_pi = 2 * np.pi / E_x  # compute t_pi the time required for the state to go from spin-up to spin-down.

        total_osc_time = self.t_pi

        # compute oscillation times.
        times_au = np.linspace(0, total_osc_time, num=time_steps)
        times_SI = np.linspace(0, self.au_to_second(total_osc_time), num=time_steps)

        B_au = np.empty((time_steps,), dtype=np.double)

        # compute the spin operators at new time
        rho_sz = kwant.operator.Density(self.syst, sigma_z, sum=True)
        rho_sy = kwant.operator.Density(self.syst, sigma_y, sum=True)
        rho_sx = kwant.operator.Density(self.syst, sigma_x, sum=True)

        density_operator = kwant.operator.Density(self.syst)
        psi = self.spin_up_state
        E_1 = np.real(self.initial_eigenvalues[0])
        E_2 = np.real(self.initial_eigenvalues[1])

        data = {
            'B_0': self.B_0_au,
            'lattice_points': self.number_of_lattices,
            'length': self.total_length_au,
            'eV_0': self.eV_0_au,
            'E_sl': self.E_sl_au,
            'E_x': E_x,
            'E_z': omega_res,
            'perturbation': self.pertubation_type,
            'potential_type': self.potential_text,
            'effective_mass': self.m_au,
            'E_1': E_1,
            'E_2': E_2,
            'states': []
        }
        print("Simulation starting with", time_steps, "time steps from 0.0 to", total_osc_time, "a.u.")

        for time in times_au:
        # time = times_au[0]
            print("Evolving state to time", time, "a.u.")
            psi.evolve(time)
            density = np.abs(self.spin_up_state.evaluate(density_operator)) ** 2
            average_eff_B_field_au = np.trapz(density * eff_B_field_au,
                                              x=self.z_au)  # compute the expectation of the slanted field.

            # compute the expectations of the spin operators on the evolved state and store their values.
            spin_z = np.real(self.spin_up_state.evaluate(rho_sz))
            spin_y = np.real(self.spin_up_state.evaluate(rho_sy))
            spin_x = np.real(self.spin_up_state.evaluate(rho_sx))

            data['states'].append({
                'time': time,
                'pdf': density.tolist(),
                'B_x': average_eff_B_field_au,
                'rho_sx': spin_x,
                'rho_sy': spin_y,
                'rho_sz': spin_z,
            })

        json_string = json.dumps(data)

        with open('/content/drive/My Drive/{}.json'.format(timestr), 'w') as outfile:
            outfile.write(json_string)

        return True

    def visualise(self):
        file_name = "20220309-135220"
        with open('./evolve-output/{}.json'.format(file_name)) as json_file:
            data = json.load(json_file)
        folder_path = "./results/{}".format(file_name)
        if os.path.isdir(folder_path) == False:
            os.mkdir(folder_path)

        times = []
        rho_sx_list = []
        rho_sy_list = []
        rho_sz_list = []
        B_x_list = []
        pdf_list = []
        parameters = {
            'B_0': self.au_to_tesla(data['B_0']),
            'lattice_points': data['lattice_points'],
            'length': self.au_to_m(data['length']),
            'eV_0': self.hartree_to_ev(data['eV_0']),
            'E_sl': self.hartree_to_ev(data['E_sl']),
            'E_x': self.hartree_to_ev(data['E_x']),
            'E_z': self.hartree_to_ev(data['E_z']),
            'perturbation': data['perturbation'],
            'potential_type': data['potential_type'],
            'effective_mass': data['effective_mass']
        }

        json_string = json.dumps(parameters)

        with open('{}/parameters.json'.format(folder_path), 'w') as outfile:
            outfile.write(json_string)

        print("Parameters saved.")

        for state in data["states"]:
            times.append(state['time'])
            rho_sx_list.append(state['rho_sx'])
            rho_sy_list.append(state['rho_sy'])
            rho_sz_list.append(state['rho_sz'])
            B_x_list.append(state['B_x'])
            pdf_list.append(state['pdf'])

        times_au = times
        times_SI = self.au_to_second(np.array(times))
        B_x_list = self.au_to_tesla(np.asarray(B_x_list))

        fontsize = 16



        spin_expectation_fig, ax = plt.subplots(figsize=(10, 6))
        plt.rcParams.update({'font.size': fontsize})
        plt.plot((times_SI), rho_sx_list, label='$\\langle X \\rangle$')
        plt.plot((times_SI), rho_sy_list, label='$\\langle Y \\rangle$')
        plt.plot((times_SI), rho_sz_list, label='$\\langle Z \\rangle$')
        plt.legend(fontsize=fontsize, loc="upper right")
        ax.set_xlabel('$t$ (s)', fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.savefig('{}/spin-expectations.eps'.format(folder_path), bbox_inches='tight')
        plt.close(spin_expectation_fig)
        print("Plot of spin expectations vs. time saved.")

        mag_field_fig, ax = plt.subplots(figsize=(10, 6))

        plt.plot(times_SI, B_x_list, label=r'$\langle B_x \rangle$')
        plt.legend(fontsize=fontsize, loc="upper right")
        ax.set_xlabel('$t$ (s)', fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        ax.set_ylabel('Magnetic Field Strength (T)', fontsize=fontsize)

        plt.savefig('{}/magnetic-fields.eps'.format(folder_path), bbox_inches='tight')
        plt.close(mag_field_fig)
        print("Plot of magnetic field vs. time saved.")

        # create a figure with two subplots
        fig_animation, axs = plt.subplots(2, 2)

        # 1920 x 1080
        w_in_inches = 19.2
        h_in_inches = 10.8
        dpi = 100
        fig_animation.set_size_inches(w_in_inches, h_in_inches, True)
        fig_animation.set_dpi(dpi)
        fig_animation.tight_layout(pad=5.0)  # add padding to subplots.

        # set variables for each of the axes.
        ax1 = axs[0, 0]
        ax2 = axs[1, 0]
        ax3 = axs[0, 1]
        ax4 = axs[1, 1]

        # intialize two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2)
        line2, = ax2.plot([], [], lw=2, color='r', label='$H_1(t)$')
        line3, = ax3.plot([], [], lw=2, label=r'$\langle X \rangle$')
        line4, = ax3.plot([], [], lw=2, label=r'$\langle Y \rangle$')
        line5, = ax3.plot([], [], lw=2, label=r'$\langle Z \rangle$')
        line6, = ax4.plot([], [], lw=2, label=r'$\langle B_x \rangle$')
        line = [line1, line2, line3, line4, line5, line6]

        z_max_SI = self.total_length_SI / 2  # in SI
        z_max_au = self.total_length_au / 2  # in SI

        y2_max = self.hartree_to_ev(
            (self.eV_0_au * z_max_au / self.total_length_au))  # in SI (maximum value for V_AC(t)

        y1 = np.abs(pdf_list[0]) ** 2
        y1_max = np.max(y1) * 2

        if data['perturbation'] == "cos":
            y3_func = self.cosine_v_ac
        else:
            y3_func = self.sine_v_ac

        # PDFs
        ax1.set_xlim(-z_max_SI, z_max_SI)
        ax1.set_ylim(0, y1_max * 1.5)
        ax1.set_ylabel("$|\psi(z,t)|^2$")

        # Potential
        ax2.set_xlim(-z_max_SI, z_max_SI)
        ax2.set_ylim(-y2_max, y2_max)

        ax2.set_ylabel("$E$ (eV)")
        ax2.set_xlabel("z ($m$)")

        # Expectations
        ax3.set_xlim(0, times_SI[-1])
        ax3.set_ylim(-1.0, 1.0)

        # Magnetic Field

        ax4.set_xlim(0, times_SI[-1])
        ax4.set_ylim(-1.1 * np.max(B_x_list), 1.1 * np.max(B_x_list))
        ax4.set_ylabel("Magnetic Field Strength (T)")
        ax4.set_xlabel("t (s)")

        # initialization function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[1].set_data([], [])
            line[2].set_data([], [])
            line[3].set_data([], [])
            line[4].set_data([], [])
            line[5].set_data([], [])

            return line

        def animate(i):
            y3 = self.hartree_to_ev(
                y3_func(times_au[i], self.z_au, self.eV_0_au, self.pulse_frequency_au, self.total_length_au))

            # update line objects.
            line[0].set_data(self.z_SI, np.abs(pdf_list[i]) ** 2)
            line[1].set_data(self.z_SI, y3)
            line[2].set_data(times_SI[0:i + 1], rho_sx_list[0:i + 1])
            line[3].set_data(times_SI[0:i + 1], rho_sy_list[0:i + 1])
            line[4].set_data(times_SI[0:i + 1], rho_sz_list[0:i + 1])
            line[5].set_data(times_SI[0:i + 1], B_x_list[0:i + 1])

            # ax1.legend(loc="upper right")
            ax2.legend(loc="upper right")
            ax3.legend(loc="upper right")
            ax4.legend(loc="upper right")

            # plt.savefig("./figures/temp/animation-{}.svg".format('{:g}'.format(float('{:.{p}g}'.format(self.au_to_second(self.time), p=2)))))

            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig_animation, animate, init_func=init,
                                       frames=len(times), interval=30, blit=True)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html

        anim.save('{}/animation.mp4'.format(folder_path), writer='ffmpeg')

        print("Animation saved.")


def main():
    lattices = [100,200,300,400,500]
    for lattice_num in lattices:
        print(lattice_num)
        potential = 1  # infinite square-well potential
        # magnetic_field_file = "simulation-0-x"
        # system = System("((A * k_z**2) + V(z, time)) * identity(2) + B(z) * sigma_z + C(z) * sigma_x + D(z) * sigma_x",
        #                 pertubation_type="cos", number_of_lattices=lattices,
        #                 potential_type=potential)  # , magnetic_fields=magnetic_field_file) # call system objecy

        system = System("((A * k_z**2) + V(z, time)) * identity(2) + B(z) * sigma_z",
                        pertubation_type="cos", number_of_lattices=lattice_num,
                        potential_type=potential)  # , magnetic_fields=magnetic_field_file) # call system objecy
        system.make_system()

        # Run these both before you evolve.
        system.initial_energies()
        # system.initial_pdfs()

    # system.evolve(100)
    # system.import_mumax3_simulations()
    # system.visualise()


if __name__ == '__main__':
    main()

