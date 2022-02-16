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
    def __init__(self, hamiltonian, lattices):
        # constants in SI
        self.time = 0
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        # self.total_length_SI = 0.66e-6
        self.total_length_SI = 0.66e-6
        self.hamiltonian = hamiltonian
        self.B_0_SI = 50e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.b_sl_SI = 1.16e6

        self.potential_type = 0
        # Constants in a.u.0
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)
        self.b_sl_au =self.tesla_SI_to_au(self.b_sl_SI) * (1/self.m_to_au(1))
        self.hbar_au = 1
        self.g = 2
        self.m_au = 1
        self.m_SI = 9.11e-31
        self.e_au = 1
        self.mu_B_SI= 9.2740100783e-24
        self.total_length_au = self.total_length_SI / self.a_0_SI # Total distance of nanotube in terms of au
        self.lattices = lattices
        self.lattice_size_au = self.total_length_au / self.lattices # The distance in atomic units spanned by a lattice point
        self.mu_B_au = .5
        self.a_0 = 1

        # self.pulse_frequency_SI = 1e10 # in Hz

        self.pulse_frequency_au =  (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au) # in Hz

        self.gaussian_mu = 0
        self.pulse_type = 2

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
        self.eV_0_au = self.ev_to_hartree(2e-6)

    def parabolic_potential(self):
        self.total_length_SI = 0.66e-6 / 2
        self.total_length_au = self.total_length_SI / self.a_0_SI # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.lattices # The distance in atomic units spanned by a lattice point

        self.E_sl_ev = 5e-4  # eV
        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.total_length_au)

        self.omega_0 = self.ev_to_hartree(1e-5)/self.hbar_au
        self.eV_0_au = self.ev_to_hartree(10e-6)


    def make_system(self):
        """
        Function to create the system
        :param length: the length of the nanotube
        :return: the system object
        """
        if self.potential_type == 1:
            self.parabolic_potential()
        else:
            self.infinite_square_well_potential()

        self.template = kwant.continuum.discretize(self.hamiltonian, grid=self.lattice_size_au)

        # We need to have 1d since the hamiltonian is 1d otherwise it can't be applied
        def shape(site):
            """
            function to define the shape of the scattering region.
            :param site: the current site.
            :return: the a boolean saying whether the scattering site should be drawn
            """

            (x, ) = site.pos
            return (-self.total_length_au/2 <= x  < self.total_length_au/2)



        self.syst = kwant.Builder()

        #Add the nanotube to the system
        self.syst.fill(self.template, shape, (0, ))

        kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        return self.syst

    def potential(self, x):  # Infinite square well
        """
        Function to define the potential of the lead.
        :param x: the position in the system.
        :return: the potential energy.
        """

        noise_potential = 0
        if self.potential_type == 0:
            self.infinite_square_well_potential()
            total_potential = 0
        elif self.potential_type == 1:
            self.parabolic_potential()
            total_potential = .5 * ((x * self.omega_0) ** 2)

        if self.time > 0:
            noise_potential = 0#(10 * self.eV_0_au) * np.random.random()

            total_potential += ((self.eV_0_au * np.cos(
                2 * np.pi * self.pulse_frequency_au * self.time) ) * x) / self.total_length_au

        return total_potential + noise_potential







    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        :param syst: the system object.
        :return: the sorted eigenvalues and eigenvectors.
        """
        self.A_constant =  -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2
        self.B_constant = -self.g * self.mu_B_au * self.b_sl_au * self.hbar_au / 2
        self.C_constant = self.hbar_au **2 / (2 * self.m_au)

        params = dict(C=self.C_constant, V=self.potential, A=self.A_constant, B=self.B_constant)

        # Calculate the wave functions in the system.
        hamiltonian = self.syst.hamiltonian_submatrix(params=params)
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

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
        eigenValues, eigenVectors = self.eigenstates()
        x_coordinates = np.linspace(-self.total_length_au/2, self.total_length_au/2, self.lattices)

        fig = plt.figure()

        # print("Simulation first level energy = ", E_1)
        # print("Defined first level energy", self.omega_0 * self.hbar_au / 2)
           # print("Simulation difference between first two levels", self.delta_E)


        energies = np.real(eigenValues)
        energies = energies[0:7]
        x = np.linspace(0,len(energies), len(energies))
        y = energies

        m, b = np.polyfit(x, y, 1)
        t = self.hbar_au ** 2 / (2 * self.m_au * (self.lattice_size_au)**2)
        print("Is the model a good approximation?", energies[2] < t)

        plt.plot(x, y, 'o')
        plt.plot(x, m * x + b)
        if self.potential_type == 0:
            potential_text = "infinite-well"
        else:
            potential_text = "parabolic"



        plt.xlabel("$n$")
        plt.ylabel("$E (eV)$")
        # plt.legend(loc="upper right")
        plt.savefig("./figures/energies/energies-{0}-potential-{1}.svg".format(self.lattices, potential_text))  # With A = 0 we expect straight forward zeeman splitting
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

        x_coords = np.linspace(-self.total_length_SI/2, self.total_length_SI/2, self.lattices)
        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.
        psi1 = eigenVectors[:, 0]
        psi1_up, psi1_down = psi1[::2], psi1[1::2]
        # even indices give the spin up and odd indices give the spin down states
        density_1 = np.abs(psi1_up) ** 2 + np.abs(psi1_down) ** 2
        psi2 = eigenVectors[:,1]
        psi2_up, psi2_down = psi2[::2], psi2[1::2]
        density_2 =  np.abs(psi2_up) ** 2 + np.abs(psi2_down) ** 2
        if animate:
            return x_coords,density_1, density_2

        else:
            fig = plt.figure()


            plt.plot(x_coords, density_1, label="$\psi_{G_{- }}$")
            plt.plot(x_coords, density_2, label="$\psi_{G_{+ }}$")

            plt.xlabel("z ($m$)")
            plt.ylabel("$|\psi(z)|^2$")
            plt.legend(loc="upper right")
            if self.potential_type == 0:
                potential_text = "infinite-well"
            else:
                potential_text = "parabolic"

            plt.savefig("./figures/pdfs/initial-pdfs-{0}-potential-{1}.svg".format(self.lattices, potential_text))  # With A = 0 we expect straight forward zeeman splitting
            plt.close(fig)
            print("Plot of wave functions saved at at t={}s".format('{:g}'.format(float('{:.{p}g}'.format(self.time, p=2)))))
            return True

    def rabi_oscillations(self, animate = False, time_steps = 50):
        self.time = 0

        eigenValues, eigenVectors = self.eigenstates()
        E_1 = np.real(eigenValues[0])
        E_2 = np.real(eigenValues[1])
        self.delta_E = E_2 - E_1
        omega_res = self.delta_E / self.hbar_au
        self.pulse_frequency_au = omega_res / (2 * np.pi)

        x_au = np.linspace(-self.total_length_au/2, self.total_length_au/2, self.lattices)
        eff_B_field_au = self.b_sl_au * x_au

        psi1 = eigenVectors[:, 0]
        psi2 = eigenVectors[:, 1]
        psi3 = eigenVectors[:, 2]
        psi4 = eigenVectors[:, 3]
        psi5 = eigenVectors[:, 4]


        total_osc_time = 4 / self.pulse_frequency_au

        times_au = np.linspace(0, total_osc_time, num=time_steps)
        B_au = np.empty((time_steps,), dtype=np.double)

        probabilities_0 = [] # Probability that it stays in initial state
        probabilities_1 = []
        probabilities_2 = []
        probabilities_3 = []
        probabilities_4 = []

        for i in range(time_steps):
            self.time = times_au[i]
            eigenValues, eigenVectors = self.eigenstates()
            psi_evolved = eigenVectors[:, 0]

            psi_evolved_up, psi_evolved_down = psi_evolved[::2], psi_evolved[1::2]
            density_1 = np.abs(psi_evolved_up) ** 2 + np.abs(psi_evolved_down) ** 2
            average_eff_B_field_au = np.trapz(density_1*eff_B_field_au, x=x_au)
            B_au[i] = np.real(average_eff_B_field_au)
            inside_term_0 = np.conjugate(psi1).T @ psi_evolved
            inside_term_1 = np.conjugate(psi2).T @ psi_evolved
            inside_term_2 = np.conjugate(psi3).T @ psi_evolved
            inside_term_3 = np.conjugate(psi4).T @ psi_evolved
            inside_term_4 = np.conjugate(psi5).T @ psi_evolved
            #
            probabilities_0.append(np.abs(inside_term_0)**2)
            probabilities_1.append(np.abs(inside_term_1)**2)
            probabilities_2.append(np.abs(inside_term_2)**2)
            probabilities_3.append(np.abs(inside_term_3)**2)
            probabilities_4.append(np.abs(inside_term_4)**2)

        if self.potential_type == 0:
            potential_text = "infinite-well"
        else:
            potential_text = "parabolic"

        if animate != True:
            # Plot magnetic field
            fig1 = plt.figure()
            plt.plot(times_au, B_au, label="$B_z$")

            plt.ylabel(r'$\bar{B}_z$ (au)')
            plt.xlabel("$t$ (au)")
            # plt.title("Plot of Average Magnetic Field Varying with Time")
            plt.savefig("./figures/magnetic-field/B-average-v-time-{0}-potential-{1}.svg".format(self.lattices, potential_text))
            plt.close(fig1)
            print("Magnetic field vs time plotted.")


            # Plot Probability
            fig2 = plt.figure()
            plt.plot(times_au, probabilities_0, label=r'$|c_1|^2$')
            plt.plot(times_au, probabilities_1, label=r'$|c_2|^2$')
            plt.plot(times_au, probabilities_2, label=r'$|c_3|^2$')
            plt.plot(times_au, probabilities_3, label=r'$|c_4|^2$')
            plt.plot(times_au, probabilities_4, label=r'$|c_5|^2$')

            plt.ylabel("Probability")
            plt.xlabel("$t$ (au)")
            plt.legend(loc="upper right")

            plt.title("Plot of Average Magnetic Field Varying with Time")
            plt.savefig("./figures/rabi-oscillations/prob-v-time-{0}-potential-{1}.svg".format(self.lattices, potential_text))
            plt.close(fig2)
            print("prob vs time plotted.")

        return probabilities_1, B_au





    def save_animation(self):
        self.time = 0
        # create a figure with two subplots
        fig, axs = plt.subplots(2, 2)
        w_in_inches = 19.2
        h_in_inches = 10.8
        fig.set_size_inches(w_in_inches, h_in_inches, True)
        dpi = 100
        fig.set_dpi(dpi)
        fig.tight_layout(pad=5.0)

        # fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1 = axs[0,0]
        ax2 = axs[1,0]
        ax3 = axs[0,1]
        ax4 = axs[1,1]

        x, y1,y2 = self.show_wave_function(animate=True)
        vpotential = vectorize(self.potential, otypes=[np.float64])
        y3 = vpotential(x)
        ax2.set_ylim(-2e-18, 2e-18)

        # intialize two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2, label='$n=1$')
        line2, = ax1.plot(x, y2, lw=2, linestyle='--',label='$n=2$')
        line3, = ax2.plot([], [], lw=2, color='r', label='$V(x)$')
        line4, = ax3.plot([], [], lw=2, color='b', label='$|c_2|^2$')
        line5, = ax4.plot([], [], lw=2, color='b', label=r'$\bar{B}_z$')


        line = [line1, line2, line3, line4,line5]
        xmax = (np.max(x))
        ax1.set_xlim(-xmax, xmax)
        ax1.set_ylim(0, np.max(y1) * 1.5)

        ax2.set_xlim(-xmax, xmax)

        plt.xlabel("z ($a_0$)")
        ax1.set_ylabel("$|\psi(x,t)|^2$")
        ax2.set_ylabel("$E$ (H)")
        ax2.set_xlabel("x ($a_0$)")

        ax3.set_ylabel("Probability")
        ax4.set_ylabel("Magnetic Field Strength (au)")
        ax4.set_xlabel("t (au)")

        number_of_frames = 120  # Takes [number_of_frames] time steps at frequency provided for pulse to travel along the
                               # the nanotube.



        probabilties, mag_field = self.rabi_oscillations(animate=True, time_steps=number_of_frames)
        total_osc_time = 4 / self.pulse_frequency_au
        times_au = np.linspace(0, total_osc_time, num=number_of_frames)

        ax3.set_xlim(0, total_osc_time)
        ax3.set_ylim(0, 1.1)
        ax4.set_ylim(-1.5e-4, 1.5e-4)
        ax4.set_xlim(0, total_osc_time)

        # initialization function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[2].set_data([], [])
            line[3].set_data([], [])
            line[4].set_data([], [])


            return line

        def animate(i):
            self.time =  (i  * total_osc_time) / number_of_frames

            x, y1,y2 = self.show_wave_function(animate=True)

            vpotential = vectorize(self.potential, otypes=[np.float64])
            y3 = vpotential(x)

            line[0].set_data(x, y1)
            line[2].set_data(x,  y3)
            line[3].set_data(times_au[0:i],  probabilties[0:i])
            line[4].set_data(times_au[0:i],  mag_field[0:i])

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
        if self.potential_type == 0:
            potential_text = "infinite-well"
        else:
            potential_text = "parabolic"

        anim.save('./figures/wavefunction-animation-{}.mp4'.format(potential_text), writer='ffmpeg')


        plt.close(fig)
        print("Animation of wave functions saved.")

def main():



    for i in range(0,1):
        lattices = 100 * i + 100
        print("Testing simulation accuracy of", lattices, "lattice points.")
        system = System("((C * k_x**2) + V(x)) * identity(2) + A * sigma_x + B  * x * sigma_z", lattices)
        system.potential_type = 0
        system.make_system()
        # system.show_wave_function()
        # system.show_energies()
        # system.rabi_oscillations()
        system.save_animation()




if __name__ == '__main__':
    main()

# Why does the parabolic potential ground state wave functions become the n=1 and n=2 wave functions when the
# B_0 is set above 250e-3? Might be Paschen-Back regime.