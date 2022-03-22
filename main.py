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
    def __init__(self, hamiltonian, pertubation_type="sin", magnetic_field_file="none", number_of_lattices=50,
                 potential_type=0):

        # Read arguments
        self.potential_type = potential_type
        self.hamiltonian = hamiltonian
        self.number_of_lattices = number_of_lattices
        self.magnetic_field_file = magnetic_field_file
        self.pertubation_type = pertubation_type
        self.evolve_state = False

        # Constants in S.I.
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_SI = 0.66e-6
        self.m_e_SI = 9.11e-31
        self.m_SI = self.m_e_SI / 18  # kg
        self.mu_B_SI = 9.2740100783e-24
        self.lattice_size_SI = self.total_length_SI / self.number_of_lattices
        self.z_SI = np.arange(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattice_size_SI,
                              dtype=np.double)  # in m

        self.B_0_SI = 5e-3  # in T
        self.g = 2

        # Units in eV
        self.E_sl_eV = 1e-6  # in eV
        self.omega_0_eV = 1e-3  # in eV
        self.eV_0_eV = 1e-3  # in eV
        # Constants in a.u.
        self.hbar_au = 1
        self.m_au = self.m_SI / self.m_e_SI
        self.e_au = 1
        self.total_length_au = self.m_to_au(self.total_length_SI)  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # Distance in a.u. between lattice points
        self.mu_B_au = .5
        self.z_au = np.arange(-self.total_length_au / 2, self.total_length_au / 2, self.lattice_size_au,
                              dtype=np.double)

        self.B_0_au = self.tesla_to_au(self.B_0_SI)
        self.E_sl_au = self.ev_to_hartree(self.E_sl_eV)
        self.b_sl_au = self.E_sl_au / (
                self.g * self.mu_B_au * self.total_length_au)  # slanted magnetic field value computed from E_sl
        self.eV_0_au = self.ev_to_hartree(self.eV_0_eV)  # the value in a.u. of eV_0.
        self.omega_0_au = self.ev_to_hartree(self.omega_0_eV)
        self.pulse_frequency_au = (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au)  # in Hz

    def cosine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        """
        Function giving the potential energy experience by an electron in a cosine alternating voltage
        :param time: current time evolution
        :param z: position along quantum dot
        :param eV_0_au: magnitude used for potential
        :param pulse_frequency_au: frequency of oscillation
        :param total_length_au: total length of quantum wire
        :return: the time-dependent pertubation due to electric field
        """
        return ((eV_0_au * np.cos(2 * np.pi * pulse_frequency_au * time)) * z) / total_length_au

    def sine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        """
           Function giving the potential energy experience by an electron in a sine alternating voltage
           :param time: current time evolution
           :param z: position along quantum dot
           :param eV_0_au: magnitude used for potential
           :param pulse_frequency_au: frequency of oscillation
           :param total_length_au: total length of quantum wire
           :return: the time-dependent pertubation due to electric field
           """
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
        """
        Function to import magnetic fields from Mumax3 simulation file
        :return: magnetic fields in different directions
        """

        file_name = self.magnetic_field_file

        # f = np.load("/content/drive/My Drive/{}".format(file_name))
        f = np.load("{}".format(file_name))

        # Effective magnetic field vectors across the wire
        B_x = f[0, 1, 0:100, 149]
        B_y = f[1, 1, 0:100, 149]
        B_z = f[2, 1, 0:100, 149]

        if self.evolve_state == False:
          fig = plt.figure()

          plt.plot(self.z_SI, B_x, label="$B_x$")
          plt.plot(self.z_SI, B_y, label="$B_z$")
          plt.plot(self.z_SI, B_z, label="$B_y$")
          plt.ylabel("Effective Magnetic Field Strength (T)")
          plt.xlabel("$z$ (m)")
          plt.legend()
          plt.savefig("./magnetic-fields.eps")
          plt.close(fig)

        return B_x, B_z, B_y

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
            total_potential = .5 * (
                    (z * self.omega_0_au) ** 2)  # define a parabolic potential inside the scattering region.
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

    def B_function(self, z):
        """
        Function to get the Hamiltonian term due to the magnetic field in the z-direction
        :param z: the position along the quantum wire
        :return: the Hamiltonian term
        """
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
            index = np.around(z, 3) == np.around(self.z_au, 3)  # generate an array of booleans where the
            # true value coincides to the position of the z-coordinate within the array
            # found 3 d.p. by testing (not sure if this is the best way)
            return -self.g * self.mu_B_au * self.tesla_to_au(
                self.B_z[index]) * self.hbar_au / 2  # use this array of booleans to extraact the relevant magnetic field
        else:
            return -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2

    def C_function(self, z):
        """
        Function to get the Hamiltonian term due to the magnetic field in the x-direction
        :param z: the position along the quantum wire
        :return: the Hamiltonian term
        """
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
            index = np.around(z, 3) == np.around(self.z_au, 3)
            return -self.g * self.mu_B_au * self.tesla_to_au(self.B_x[index]) * self.hbar_au / 2
        else:
            return -self.b_sl_au * z

    def D_function(self, z):
        """
        Function to get the Hamiltonian term due to the magnetic field in the y-direction
        :param z: the position along the quantum wire
        :return: the Hamiltonian term
        """
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
            index = np.around(z, 3) == np.around(self.z_au, 3)
            return -self.g * self.mu_B_au * self.tesla_to_au(self.B_y[index]) * self.hbar_au / 2
        else:
            return 0

    def make_system(self):
        """
        Function to create the system
        :param length: the length of the nanotube
        :return: the system object
        """
        if self.potential_type == 1:  # if we want a parabolic potential
            self.potential_text = "parabolic"

        else:  # if infinite square well
            self.potential_text = "infinite-well"

        self.template = kwant.continuum.discretize(self.hamiltonian,
                                                   grid=self.lattice_size_au)  # make template for lattice points based on the inputted hamiltonian/

        self.syst = kwant.Builder()

        # Add the nanotube to the system
        self.syst.fill(self.template, self.kwant_shape, (0,))

        # kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        self.A_constant = self.hbar_au ** 2 / (2 * self.m_au)  # coefficient for the kinetic energy term
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
            self.B_x, self.B_y, self.B_z = self.import_mumax3_simulations()  # import the magnetic fields



        # import these function and coefficients for use in the full Hamiltonian used to define the system
        self.params = dict(A=self.A_constant, V=self.potential, B=self.B_function, C=self.C_function, D=self.D_function)

        self.tparams = self.params.copy()  # copy the params array
        self.tparams['time'] = 0  # add another parameter, with the initial time = 0
        print("System intialised.")

        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=self.tparams)
        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        # Sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        self.initial_eigenvalues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        # initial wave functions unperturbed by time dependent part of hamiltonian
        self.psi_1_init = eigenVectors[:, 0]
        self.psi_2_init = eigenVectors[:, 1]

        # an object representing the spin-up state (we care about how this state evolves with time)
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

        return eigenValues, eigenVectors

    def initial_pdfs(self):

        """
        Function to show the initial probability density functions of the spin-up and spin-down ground state
        :param syst: the system object.
        :return: PDFs of the spin-up and spin-down ground state
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
        plt.savefig("./initial-pdfs.eps")  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of PDFs at t=0 saved.")

        return density_1, density_2

    def initial_energies(self):
        """
        Function to display energy levels of the system in eV
        :param syst: the system object.
        :return:
        """

        fig = plt.figure()

        y = self.hartree_to_ev(np.real(self.initial_eigenvalues))
        print("E_1 is", y[0], "eV.")
        print("E_2 is", y[1], "eV.")

        # Plot the energies of the different levels.
        plt.plot([0, 1], [y[0], y[0]], label=r"$G_-$")
        plt.plot([0, 1], [y[1], y[1]], label=r"$G_+$")
        plt.plot([0, 1], [y[2], y[2]], label=r"$E_-$")
        plt.plot([0, 1], [y[3], y[3]], label=r"$E_+$")
        plt.legend(loc="best")
        plt.ylabel("$E$ (a.u.)")

        plt.savefig("./energies.eps")  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of eigenenergies at t=0 saved.")
        return y

    def evolve(self, time_steps=1000):
        timestr = timelib.strftime("%Y%m%d-%H%M%S")
        self.evolve_state = True

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

        total_osc_time = 0.01 * self.t_pi

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

        # data dictionary to store parameters used in simulation (units are a.u.)
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
            print("Evolving state to time", time, "a.u.")
            psi.evolve(time)  # evolved the wavefunction according to TDSE to time
            density = np.abs(self.spin_up_state.evaluate(density_operator)) ** 2  # compute PDF
            average_eff_B_field_au = np.trapz(density * eff_B_field_au,
                                              x=self.z_au)  # compute the expectation of the slanted field.
            average_eff_B_field_au = np.trapz(density * eff_B_field_au,
                                              x=self.z_au)  # compute the expectation of the slanted field.

            # compute the expectations of the spin operators on the evolved state and store their values.
            spin_z = np.real(self.spin_up_state.evaluate(rho_sz))
            spin_y = np.real(self.spin_up_state.evaluate(rho_sy))
            spin_x = np.real(self.spin_up_state.evaluate(rho_sx))

            # for each time step add to the relevant quantities at this time to the data['states'] list
            data['states'].append({
                'time': time,
                'pdf': density.tolist(),
                'B_x': average_eff_B_field_au,
                'B_y': average_eff_B_field_au,
                'B_z': average_eff_B_field_au,

                'rho_sx': spin_x,
                'rho_sy': spin_y,
                'rho_sz': spin_z,
            })

        json_string = json.dumps(data)

        # save the simulation data
        # with open('/content/drive/My Drive/{}.json'.format(timestr), 'w') as outfile:
        #     outfile.write(json_string)

        with open('{}.json'.format(timestr), 'w') as outfile:
            outfile.write(json_string)

        self.evolve_state = False

        return True

    def visualise(self, file_name):
        """
        Function to visualise the output of the EDSR simulator
        :param file_name: The name of the file you which to import and visualise
        :return:
        """

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
        parameters = {  # from the imported data extract the relevant params and convert to SI units
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
            outfile.write(json_string)  # save the params

        print("Parameters saved.")

        for state in data["states"]:  # for each of the different states/times
            # extract the quantities and put them in their own lists
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

        # plot the pauli spin matrices expectation values
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

        # plot the magnetic field in the x-direction
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

            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig_animation, animate, init_func=init,
                                       frames=len(times), interval=30, blit=True)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html

        # save animation
        anim.save('{}/animation.mp4'.format(folder_path), writer='ffmpeg')

        print("Animation saved.")


def main():
    lattices = 100  # number of lattice points

    potential = 0  # infinite square-well potential
    magnetic_field_file = "B_eff000000.npy"
    system = System("((A * k_z**2) + V(z, time)) * identity(2) + B(z) * sigma_z + C(z) * sigma_x + D(z) * sigma_y",
                    pertubation_type="sin", number_of_lattices=lattices,
                    potential_type=potential, magnetic_field_file=magnetic_field_file)  # call system objecy
    system.make_system()

    # Run these both before you evolve.
    system.initial_energies()
    system.initial_pdfs()

    system.evolve(100)
    # # system.import_mumax3_simulations()
    # system.visualise("20220321-160930")


if __name__ == '__main__':
    main()

