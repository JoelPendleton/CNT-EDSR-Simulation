import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
from numpy import vectorize
from matplotlib import animation
from tkwant import onebody

from csv import writer
def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

class System:
    def __init__(self, hamiltonian, lattices=50, potential_type = 0):

        # Define units
        self.evolve = False
        self.potential_type = potential_type
        self.hamiltonian = hamiltonian
        self.lattices = lattices

        # Constants in SI units
        self.time = 0
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_SI = 0.66e-6
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.b_sl_SI = 1.16e6
        self.m_SI = 9.11e-31  # divide by ten to get the effective mass
        self.mu_B_SI= 9.2740100783e-24
        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattices)

        # Constants in a.u.
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)
        self.b_sl_au =self.tesla_SI_to_au(self.b_sl_SI) * (1/self.m_to_au(1))
        self.hbar_au = 1
        self.g = 2
        self.m_au = 1
        self.e_au = 1
        self.total_length_au = self.total_length_SI / self.a_0_SI # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.lattices # The distance in atomic units spanned by a lattice point
        self.mu_B_au = .5
        self.a_0 = 1
        self.z_au =  np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.lattices)
        self.pulse_frequency_au =  (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au) # in Hz


    def import_mumax3_simulations(self):

        file_name = "simulation-0-y"

        f = np.load("./magnetic-field-simulations/Simulation-16-02-22/{}.out/B_eff000000.npy".format(file_name))
        y_coords = list(range(0, 66))

        # y_indices = np.arange(8, 74, 1)
        # print(return_indices(y_indices))
        # Effective magnetic field vectors across the wire

        # For each y need to add on some z-shift

        B_x = f[0, 14:19, 8:74, 49:54]
        B_y = f[1, 14:19, 8:74, 49:54]
        B_z = f[2, 14:19, 8:74, 49:54]

        # Effective magnetic field mean magnitude at each y-coord
        B_x_mean = np.mean(B_x, axis=(0, 2))
        B_y_mean = np.mean(B_y, axis=(0, 2))
        B_z_mean = np.mean(B_z, axis=(0, 2))



        fig = plt.figure()


        plt.plot(y_coords, B_x_mean, label="$B_x$")
        # plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
        plt.plot(y_coords, B_y_mean, label="$B_y$")
        # plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
        plt.plot(y_coords, B_z_mean, label="$B_z$")
        # plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
        plt.ylabel("$B_{eff}$ $(T)$")
        plt.xlabel("$y$ $(10^{-8}m)$")
        plt.legend()
        plt.savefig("./figures/{}.svg".format(file_name))
        plt.close(fig)

    def tesla_SI_to_au(self, tesla):
        """
        Function to convert the magnetic flux density from SI units to AU units
        :param tesla: the magnetic flux density in teslas.
        :return: the magnetic flux density in AU.
        """
        return tesla / 2.35e5

    def tesla_au_to_SI(self, tesla):
        """
        Function to convert the magnetic flux density from SI units to AU units
        :param tesla: the magnetic flux density in teslas.
        :return: the magnetic flux density in AU.
        """
        return tesla * 2.35e5

    def time_SI_to_au(self, time):
        return time * 4.1341373336493e16

    def time_au_to_SI(self, time):
        return time / 4.1341373336493e16

    def hartree_to_ev(self,hartree):
        return hartree * 2.72114e1
    def ev_to_hartree(self,ev):
        return ev / 2.72114e1

    def au_to_m(self, au):
        return 5.2917721090380e-11 * au

    def m_to_au(self, m):
        return m / 5.2917721090380e-11

    def hz_to_au(self, hz):

        return hz * 1.51983e-16

    def au_to_hz(self, au):

        return au / 1.51983e-16

    def infinite_square_well_potential(self):

        self.total_length_SI = 0.66e-6 / 4
        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattices)
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)

        self.E_sl_ev = 1e-6

        self.total_length_au = self.total_length_SI / self.a_0_SI  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.lattices  # The distance in atomic units spanned by a lattice point
        self.z_au = np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.lattices)

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.total_length_au) # slanted magnetic field value computed from E_sl
        self.eV_0_au = self.ev_to_hartree(100e-6) # the value in a.u. of eV_0.

        print("b_sl is", self.tesla_au_to_SI(self.b_sl_au) / self.au_to_m(1), "T.")
        print("eV_0 is", self.hartree_to_ev(self.eV_0_au), "eV.")
        print("E_sl is", self.E_sl_ev, "eV.")


    def parabolic_potential(self):
        self.total_length_SI = 0.66e-6
        self.confinement_length_SI = 0.33e-6
        self.confinement_length_au = self.m_to_au(self.confinement_length_SI)

        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattices)
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.E_sl_ev = 1e-6

        self.total_length_au = self.total_length_SI / self.a_0_SI  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.lattices  # The distance in atomic units spanned by a lattice point
        self.z_au = np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.lattices)

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.confinement_length_au) # slanted magnetic field value computed from E_sl
        self.omega_0 = self.hbar_au/(self.m_au*(self.confinement_length_au)**2)
        self.eV_0_au = self.ev_to_hartree(10e-6) # the value in a.u. of eV_0.

        print("b_sl is", self.tesla_au_to_SI(self.b_sl_au) / self.au_to_m(1), "T.")
        print("hbar * omega_0 is", self.hartree_to_ev(self.omega_0), "eV.")
        print("eV_0 is", self.hartree_to_ev(self.eV_0_au), "eV.")
        print("E_sl is", self.E_sl_ev, "eV.")



    def make_system(self):
        """
        Function to create the system
        :param length: the length of the nanotube
        :return: the system object
        """

        if self.potential_type == 1: # if we want a parabolic potential
            self.potential_text = "parabolic"
            self.parabolic_potential() # define constants using parabolic potential function
        else: # if infinite square well
            self.potential_text = "infinite-well"
            self.infinite_square_well_potential() # define constants using hard-wall potential function

        self.template = kwant.continuum.discretize(self.hamiltonian, grid=self.lattice_size_au) # make template for lattice points based on the inputted hamiltonian/

        def shape(site):
            """
            function to define the shape of the scattering region.
            :param site: the current site.
            :return: the a boolean saying whether the scattering site should be drawn
            """
            (z, ) = site.pos
            return (-self.total_length_au/2 <= z  < self.total_length_au/2)



        self.syst = kwant.Builder()

        #Add the nanotube to the system
        self.syst.fill(self.template, shape, (0, ))

        kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        return self.syst

    def V_AC(self, z):
        """
        Function defining the alternating potential pertubation (H_1) to be added to the hamiltonian.
        :param z: the positions along the CNT.
        :return: the value of the perturbartion.
        """
        return ((self.eV_0_au * np.sin(
            2 * np.pi * self.pulse_frequency_au * self.time)) * z) / self.total_length_au

    def potential(self, z):  # Infinite square well
        """
        Function to define the potential of the lead.
        :param x: the position in the system.
        :return: the potential energy.
        """

        if self.potential_type == 0: # if hard-wall potential is used
            total_potential = 0 # define zero for the potential inside the scattering region.
            # outside the potential will be infinity
        elif self.potential_type == 1: # for a parabolic potential
            total_potential = .5 * ((z * self.omega_0) ** 2) # define a parabolic potential inside the scattering region.

        if self.evolve: # if the state is evolving
            total_potential += self.V_AC(z) # add a cosine term to the potential to represent an applied E field.

        return total_potential








    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        :param syst: the system object.
        :return: the sorted eigenvalues and eigenvectors.
        """
        # Constants used to define the Hamiltonian (using a.u.)
        self.A_constant =  -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2
        self.B_constant = -self.g * self.mu_B_au * self.b_sl_au * self.hbar_au / 2
        self.C_constant = self.hbar_au **2 / (2 * self.m_au)

        params = dict(C=self.C_constant, V=self.potential, A=self.A_constant, B=self.B_constant)

        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=params)

        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        # Sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        return eigenValues, eigenVectors

    def show_energies(self):
        """
        Procedure to display the potential and energy levels of the system
        :param syst: the system object.
        :return:
        """

        # compute the eigenvalues
        eigenValues, eigenVectors = self.eigenstates()

        fig = plt.figure()

        energies = np.real(eigenValues)
        energies = energies[0:10]
        y = energies
        print("E_1 is",  self.hartree_to_ev(y[0]), "eV.")

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
        plt.plot([0, 1], [y[0], y[0]],  label=r"$G_-$")
        plt.plot([0, 1], [y[1],y[1]], label=r"$G_+$")
        plt.plot([0, 1], [y[2],y[2]], label=r"$E_-$")
        plt.plot([0, 1], [y[3],y[3]], label=r"$E_+$")
        plt.legend(loc="best")
        plt.ylabel("$E$ (a.u.)")

        plt.savefig("./figures/energies/energies-{0}-potential-{1}.svg".format(self.lattices, self.potential_text))  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of eigenenergies saved.")
        return True

    def show_wave_function(self, animate=False):

        """
        Procedure to show the probability density function.
        :param syst: the system object.
        :return:
        """
        eigenValues, eigenVectors = self.eigenstates()

        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.
        psi1 = eigenVectors[:, 0]
        psi1_up, psi1_down = psi1[::2], psi1[1::2]
        # even indices give the spin up and odd indices give the spin down states
        density_1 = np.abs(psi1_up) ** 2 + np.abs(psi1_down) ** 2
        psi2 = eigenVectors[:,1]
        psi2_up, psi2_down = psi2[::2], psi2[1::2]
        density_2 =  np.abs(psi2_up) ** 2 + np.abs(psi2_down) ** 2

        if animate == False: # if we do not wish to animate our wave function and wish to plot it.
            fig = plt.figure()


            plt.plot(self.z_SI, density_1, label=r'$|1, -\frac{1}{2}\rangle$')
            plt.plot(self.z_SI, density_2, label=r'$|1, +\frac{1}{2}\rangle$')

            plt.xlabel("$x$ (a.u.)")
            plt.ylabel("$|\psi(x)|^2$")
            plt.legend(loc="upper right")

            # save file with name according to potential type.
            plt.savefig("./figures/pdfs/initial-pdfs-{0}-potential-{1}.svg".format(self.lattices, self.potential_text))  # With A = 0 we expect straight forward zeeman splitting
            plt.close(fig)
            print("Plot of wave functions saved at at t={}s".format('{:g}'.format(float('{:.{p}g}'.format(self.time, p=2)))))

        return density_1, density_2

    def rabi_oscillations(self, animate = False, time_steps = 20, animation_osc_time = 1e-8):
        self.time = 0 # initialise system to time = 0

        def x_onsite(site): # function to compute the position operator matrix.
            return [site.pos[0]]*np.identity(2)

        # spin matrices
        sigma_x = np.array([[0, 1],
                            [1, 0]])
        sigma_y = np.array([[0, 1j],
                            [-1j, 0]])
        sigma_z = np.array([[1, 0],
                            [0, -1]])

        # compute eigenstates
        eigenValues, eigenVectors = self.eigenstates()


        # extract lowest two energies (our qubit states)
        E_1 = np.real(eigenValues[0])
        E_2 = np.real(eigenValues[1])

        # extract the state vectors corresponding to these lowest eigenenergies.
        psi1 = eigenVectors[:, 0]
        psi2 = eigenVectors[:, 1]

        # compute the difference in these energies
        self.delta_E = np.abs(E_2 - E_1)



        # print("ratio of E_z/w0 = ", self.delta_E/self.omega_0)
        # compute the resonant frequency for the rabi oscillations
        omega_res = self.delta_E / self.hbar_au
        self.pulse_frequency_au = omega_res / (2 * np.pi) # set the frequency of the V_AC(t) to be equal to the res. freq.

        eff_B_field_au = self.b_sl_au * self.z_au # compute the slanted field at each point along the CNT.


        # define the density operator of x
        rho_x = kwant.operator.Density(self.syst, x_onsite, sum=True)
        # compute this density operator on the ground states <1|x|2>:
        rho_x_1_2 = rho_x(psi1,psi2)

        # compute the energy E_x
        E_x = np.real(2 * self.eV_0_au *  rho_x_1_2 / self.total_length_au)

        self.t_pi = 2*np.pi/E_x # compute t_pi the time required for the state to go from spin-up to spin-down.
        # print("E_x is", E_x)
        # print("E_z is", omega_res)
        # print(self.time_au_to_SI(self.t_pi))

        if animate == True: # if we want to animate
            total_osc_time = animation_osc_time
        else: # if we don't wish to animate

            total_osc_time = 2 * self.t_pi

        # compute oscillation times.
        times_au = np.linspace(0, total_osc_time, num=time_steps)
        times_SI = np.linspace(0, self.time_au_to_SI(total_osc_time), num=time_steps)

        B_au = np.empty((time_steps,), dtype=np.double)

        probabilities_0 = [] # Probability that it stays in initial state
        probabilities_1 = [] # Probability that it enters the excited state

        # Expectation values of spin operators
        rho_sz_list = []
        rho_sy_list = []
        rho_sx_list = []

        self.evolve = True # allow the system to evolve.


        for i in range(time_steps):

            self.time = times_au[i] # update the time

            # compute the spin operators at new time
            rho_sz = kwant.operator.Density(self.syst, sigma_z, sum=True)
            rho_sy = kwant.operator.Density(self.syst, sigma_y, sum=True)
            rho_sx = kwant.operator.Density(self.syst, sigma_x, sum=True)

            # update eigenstates
            eigenValues, eigenVectors = self.eigenstates()

            psi_evolved = eigenVectors[:, 0] # compute the evolved spin-up state

            psi_evolved_up, psi_evolved_down = psi_evolved[::2], psi_evolved[1::2] # isolate the comprising spin components
            density_1 = np.abs(psi_evolved_up) ** 2 + np.abs(psi_evolved_down) ** 2 # compute the probability density

            average_eff_B_field_au = np.trapz(density_1*eff_B_field_au, x=self.z_au) #compute the expectation of the slanted field.
            B_au[i] = np.real(average_eff_B_field_au) # store this value in the numpy array

            inside_term_0 = np.conjugate(psi1).T @ psi_evolved # compute the overlap of the evolved state with the initial spin-up state
            inside_term_1 = np.conjugate(psi2).T @ psi_evolved # compute the overlap of the evolved state with the initial spin-down state
            # add these to their respective lists.
            probabilities_0.append(np.abs(inside_term_0)**2)
            probabilities_1.append(np.abs(inside_term_1)**2)

            # compute the expecrations of the spin operators on the evolved state and store their values.
            spin_z = np.real(rho_sz(psi_evolved))
            rho_sz_list.append(spin_z)
            spin_y = np.real(rho_sy(psi_evolved))
            rho_sy_list.append(spin_y)
            spin_x = np.real(rho_sx(psi_evolved))
            rho_sx_list.append(spin_x)

        if animate == False: # if not animating

            # Plot magnetic field
            fig1 = plt.figure()
            plt.plot(times_SI, self.tesla_au_to_SI(B_au), label="$B_x$")
            plt.ylabel(r'$\bar{B}_z$ (T)')
            plt.xlabel("$t$ (s)")
            # plt.title("Plot of Average Magnetic Field Varying with Time")
            plt.savefig("./figures/magnetic-field/B-average-v-time-{0}-potential-{1}.svg".format(self.lattices, self.potential_text))
            plt.close(fig1)
            np.save('magnetic_field.npy', B_au)
            np.save('times.npy', times_au)
            print("Magnetic field vs time plotted.")


            # Plot Expectations
            fig2 = plt.figure()
            plt.plot(times_SI, rho_sx_list, label=r'$\rho_x$')
            plt.plot(times_SI, rho_sy_list, label=r'$\rho_y$')
            plt.plot(times_SI, rho_sz_list, label=r'$\rho_z$')

            plt.xlabel("$t$ (s)")
            plt.legend(loc="upper right")
            plt.savefig("./figures/rabi-oscillations/spin-expectations-{0}-potential-{1}.svg".format(self.lattices, self.potential_text))
            plt.close(fig2)
            print("expect vs time plotted.")

            # Plot Probabilities
            fig3 = plt.figure()
            plt.plot(times_SI, probabilities_0, label=r'$|c_1|^2$')
            plt.plot(times_SI, probabilities_1, label=r'$|c_2|^2$')
            plt.axvline(x=self.time_au_to_SI(self.t_pi), linestyle='--', label=r"$t_{\pi}$")

            plt.ylabel("Probability")
            plt.xlabel("$t$ (s)")
            # plt.legend(loc="upper right")
            plt.savefig(
                "./figures/rabi-oscillations/prob-v-time-{0}-potential-{1}.svg".format(self.lattices, self.potential_text))
            plt.close(fig3)
            print("prob vs time plotted.")

        return probabilities_1, B_au

    def save_animation(self):

        self.time = 0 # intialise system to time = 0

        # create a figure with two subplots
        fig, axs = plt.subplots(2, 2)

        # 1920 x 1080
        w_in_inches = 19.2
        h_in_inches = 10.8
        dpi = 100
        fig.set_size_inches(w_in_inches, h_in_inches, True)
        fig.set_dpi(dpi)
        fig.tight_layout(pad=5.0) # add padding to subplots.

        # set variables for each of the axes.
        ax1 = axs[0,0]
        ax2 = axs[1,0]
        ax3 = axs[0,1]
        ax4 = axs[1,1]

        y1,y2 = self.show_wave_function(animate=True) # x is in au


        # intialize two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2, label=r'$|1, +\frac{1}{2}\rangle$')
        line2, = ax1.plot(self.z_SI, y2, lw=2, label=r'$|1, -\frac{1}{2}\rangle$')
        line3, = ax2.plot([], [], lw=2, color='r', label='$V(z,t)$')
        line4, = ax3.plot([], [], lw=2, color='b', label='$|c_2|^2$')
        line5, = ax4.plot([], [], lw=2, label=r'$\bar{B}_z$')
        line = [line1, line2, line3, line4,line5]

        z_max_SI = self.total_length_SI / 2 # in SI
        z_max_au = self.total_length_au / 2 # in SI

        y2_max = self.hartree_to_ev(self.eV_0_au  * z_max_au / self.total_length_au) # in SI (maximum value for V_AC(t)
        number_of_frames = 20 # number of frames to use for animation -> animation is 30 fps so runs for number_of_frames/30 seconds.

        total_osc_time_au = self.t_pi  # runs for t_pi
        total_osc_time_SI = self.time_au_to_SI(total_osc_time_au)
        times_SI = np.linspace(0, total_osc_time_SI, num=number_of_frames)

        # get the probabilities for these times.
        probabilties, mag_field_au = self.rabi_oscillations(animate=True, time_steps=number_of_frames, animation_osc_time=total_osc_time_au)
        # fig.suptitle("$ \omega  =${0} rads^{-1}, e$V_0$ ={1} eV".format(self.au_to_hz(self.pulse_frequency_au), self.hartree_to_ev(self.eV_0_au)))


        # Set limits and labels for plots.

        # PDFs
        ax1.set_xlim(-z_max_SI, z_max_SI)
        ax1.set_ylim(0, np.max(y1) * 1.5)
        ax1.set_ylabel("$|\psi(z,t)|^2$")

        # Potential
        ax2.set_xlim(-z_max_SI, z_max_SI)
        ax2.set_ylim(-y2_max, y2_max)

        ax2.set_ylabel("$E$ (eV)")
        ax2.set_xlabel("z ($m$)")

        # Probability
        ax3.set_xlim(0, total_osc_time_SI)
        ax3.set_ylim(0, 1.1)
        ax3.set_ylabel("Probability")

        # Magnetic Field

        ax4.set_xlim(0, total_osc_time_SI)
        # mag_field_limits = self.tesla_au_to_SI(-1.5e-4)
        mag_field_SI = self.tesla_au_to_SI(mag_field_au)
        ax4.set_ylim(-100, 100)
        ax4.set_ylabel("Magnetic Field Strength (T)")
        ax4.set_xlabel("t (s)")

        # initialization function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[2].set_data([], [])
            line[3].set_data([], [])
            line[4].set_data([], [])


            return line

        def animate(i):

            # update times
            self.time =  (i  * total_osc_time_au) / number_of_frames

            # update wavefunctions
            y1,y2 = self.show_wave_function(animate=True)

            # update V_AC
            V_AC = vectorize(self.V_AC, otypes=[np.float64])
            y3 = self.hartree_to_ev(V_AC(self.z_au))

            # update line objects.
            line[0].set_data(self.z_SI, y1)
            line[2].set_data(self.z_SI,  y3)
            line[3].set_data(times_SI[0:i],  probabilties[0:i])
            line[4].set_data(times_SI[0:i],  mag_field_SI[0:i])

            ax1.legend(loc="upper right")
            ax2.legend(loc="upper right")
            ax3.legend(loc="upper right")
            ax4.legend(loc="upper right")

            # plt.savefig("./figures/temp/animation-{}.svg".format('{:g}'.format(float('{:.{p}g}'.format(self.time_au_to_SI(self.time), p=2)))))

            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate,init_func=init,
                                       frames=number_of_frames, interval=30, blit=True)



        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html

        anim.save('./figures/wavefunction-animation-{}.mp4'.format(self.potential_text), writer='ffmpeg')


        plt.close(fig)
        print("Animation of wave functions saved.")

def main():

    lattices = 100 # number of lattice points
    print("Testing simulation accuracy of", lattices, "lattice points.")
    system = System("((C * k_z**2) + V(z)) * identity(2) + A * sigma_z + B  * z * sigma_x", lattices, potential_type=0) # call system objecy
    system.make_system() # make the system
    system.time = 0 # initialise the time to zero.
    system.show_wave_function()
    system.show_energies()
    system.rabi_oscillations()
    # system.save_animation()
    # system.import_mumax3_simulations()


if __name__ == '__main__':
    main()

# Why does the parabolic potential ground state wave functions become the n=1 and n=2 wave functions when the
# B_0 is set above 250e-3? Might be Paschen-Back regime.