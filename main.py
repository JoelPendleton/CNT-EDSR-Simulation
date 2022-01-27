import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
from numpy import vectorize
from matplotlib import animation


class System:
    def __init__(self, hamiltonian, lattices):
        # constants in SI
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_SI =7820e-9
        self.B_0_SI = 50e-3
        self.b_sl_SI = 50e-3
        self.time = 0

        # Constants in a.u.
        self.hbar_au = 1
        self.g = 2
        self.m_au = 1
        self.m_SI = 9.11e-31

        self.e_au = 1
        self.mu_B_SI= 9.2740100783e-24
        self.total_length_au = self.total_length_SI / self.a_0_SI # Total distance of nanotube in terms of au
        self.lattices = lattices
        self.lattice_size_au = self.total_length_au / self.lattices # The distance in atomic units spanned by a lattice point
        self.mu_B_au = self.e_au * self.hbar_au / (2 * self.m_au)
        self.a_0 = 1
        self.gaussian_height = 0.0014699721440278707

        self.pulse_frequency_SI = 1e9 # in Hz
        self.osc_time_SI = (1/self.pulse_frequency_SI)
        self.osc_time_au = (1/ (self.osc_time_SI/2.419e-17))
        self.pulse_frequency_au = 1/self.osc_time_au
        self.pulse_velocity = self.pulse_frequency_au * self.total_length_au # in au/s
        self.hamiltonian = hamiltonian


        self.template = kwant.continuum.discretize(hamiltonian, grid=self.lattice_size_au)
        self.gaussian_mu = 0
        self.pulse_type = 1

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
        return hartree * 27.2114
    def ev_to_hartree(self,ev):
        return ev / 27.2114

    def au_to_nm(self, au):
        return 0.0529177249 * au

    def au_to_m(self, au):
        return 5.2917721090380e-11 * au

    def m_to_au(self, m):
        return m / 5.2917721090380e-11

    def make_system(self):
        """
        Function to create the system
        :param length: the length of the nanotube
        :return: the system object
        """
        # We need to have 1d since the hamiltonian is 1d otherwise it can't be applied
        def shape(site):
            """
            function to define the shape of the scattering region.
            :param site: the current site.
            :return: the a boolean saying whether the scattering site should be drawn
            """

            (x, ) = site.pos
            return (0 <= x < self.total_length_au)

        self.syst = kwant.Builder()

        #Add the nanotube to the system
        self.syst.fill(self.template, shape, (0, ))

        kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        return self.syst



    def gaussian(self, x, sig):
        """
        Function to compute a gaussian pulse.
        :param x: the coordinate value
        :param mu: the centre of the gaussian
        :param sig: the standard deviation
        :return: the gaussian function values.
        """

        return np.exp(-np.power(x - self.gaussian_mu - self.pulse_velocity * self.time , 2.) / (2 * np.power(sig, 2.)))

    def potential(self, x):  # Infinite square well
        """
        Function to define the potential of the lead.
        :param x: the position in the system.
        :return: the potential energy.
        """
        if -self.total_length_au/2 <= x <= self.total_length_au/2:

            if self.pulse_type == 1:
                self.gaussian_mu = -self.total_length_au/2
                self.pulse_velocity = self.pulse_frequency_au * self.total_length_au  # in au/s
                return self.gaussian_height * self.gaussian(x, self.total_length_au / 14)
            elif self.pulse_type == 2:
                self.gaussian_mu = -self.total_length_au/4
                self.pulse_velocity = 0
                return self.gaussian_height * (np.cos(2 * np.pi * self.pulse_frequency_au * self.time) **2) * self.gaussian(x, self.total_length_au / 14)
            else:

                return  0

        else:
            return 999999999

    def kwant_potential(self, x):
        """
        Potential that kwant uses. It inputs lattice points coord. so we need to convert to au to use our
        other functions
        :param x:
        :return: potential after conversion to au
        """
        return self.potential(x - self.total_length_au/2) # Kwant uses lattice points so need scaling

    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        :param syst: the system object.
        :return: the sorted eigenvalues and eigenvectors.
        """
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)
        self.b_sl_au = self.tesla_SI_to_au(self.b_sl_SI) * (1/self.m_to_au(1)) # need to scale length unit too
        self.A_constant =  -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2
        self.B_constant = -self.g * self.mu_B_au * self.b_sl_au * self.hbar_au / 2
        self.C_constant = self.hbar_au ** 2 / (2 * self.m_au)
        params = dict(C=self.C_constant, V=self.kwant_potential, A=self.A_constant, B=self.B_constant)
        # Calculate the wave functions in the system.
        hamiltonian = self.syst.hamiltonian_submatrix(params=params)
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        # For each eigenvalue for each spin we get two eigenvectors and for each position there are two possible spins
        # so we have four eigenvectors for each position along the nanotube.


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
        E_1 = self.hartree_to_ev(np.real(eigenValues[0]))
        E_2 = self.hartree_to_ev(np.real(eigenValues[1]))
        E_3 = self.hartree_to_ev(np.real(eigenValues[2]))
        E_4 = self.hartree_to_ev(np.real(eigenValues[3]))
        E_5 = self.hartree_to_ev(np.real(eigenValues[4]))

        # print("The energies of the ground and excited states are {0} and {1}, respectively.".format(E_1, E_2))
        plt.plot([0, np.max(x_coordinates)], [E_1, E_1], label="$E_1$")
        plt.plot([0, np.max(x_coordinates)], [E_2, E_2], label="$E_2$")
        plt.plot([0, np.max(x_coordinates)], [E_3, E_3], label="$E_3$", linestyle='dashed')
        plt.plot([0, np.max(x_coordinates)], [E_4, E_4], label="$E_4$", linestyle='dashed')
        plt.plot([0, np.max(x_coordinates)], [E_5, E_5], label="$E_5$", linestyle='dashed')
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        # plt.xlabel("$x (au$)")
        plt.ylabel("$E (eV)$")
        plt.legend(loc="upper right")
        plt.savefig("./figures/energies.svg")  # With A = 0 we expect straight forward zeeman splitting
        plt.close(fig)
        print("Plot of eigenenergies saved.")

    def show_wave_function(self, animate=False):

        """
        Procedure to show the probability density function.
        :param syst: the system object.
        :return:
        """

        eigenValues, eigenVectors = self.eigenstates()

        x_coords = np.linspace(-self.total_length_au/2, self.total_length_au/2, self.lattices)

        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.
        psi1 = eigenVectors[:, 0]
        psi1_up, psi1_down = psi1[::2], psi1[1::2]
        # even indices give the spin up and odd indices give the spin down states
        density_1 = np.abs(psi1_up) ** 2 + np.abs(psi1_down) ** 2
        psi2 = eigenVectors[:, 1]
        psi2_up, psi2_down = psi2[::2], psi2[1::2]
        density_2 = np.abs(psi2_up) ** 2 + np.abs(psi2_down) ** 2
        if animate:
            return x_coords,density_1, density_2

        else:
            fig = plt.figure()
            plt.plot((x_coords), density_1, label="$\psi_{G_{- }}$")
            plt.plot((x_coords), density_2, label="$\psi_{G_{+ }}$")

            plt.xlabel("z ($a_0$)")
            plt.ylabel("$|\psi(z)|^2$")
            plt.legend(loc="upper right")
            plt.savefig("./figures/initial-pdfs.svg".format(self.time))  # With A = 0 we expect straight forward zeeman splitting
            plt.close(fig)
            print("Plot of wave functions saved at at t={}s".format('{:g}'.format(float('{:.{p}g}'.format(self.time, p=2)))))
            return True

    def rabi_oscillations(self):
        self.pulse_type = 2
        self.sigma_y = np.array([[0, -1j],
                                 [1j, 0]])
        self.sigma_z = np.array([[1, 0],
                                 [0, -1]])
        self.sigma_x = np.array([[0, 1],
                                 [1, 0]])
        eigenValues, eigenVectors = self.eigenstates()
        E_1 = np.real(eigenValues[0])
        E_2 = np.real(eigenValues[1])
        w_12 = (E_2 - E_1)/self.hbar_au
        # Compte w_12 for B_0 eigenvalues
        # self.pulse_frequency_au = (w_12) / (2 * np.pi)

        x_au = np.linspace(-self.total_length_au/2, self.total_length_au/2, self.lattices)

        eff_B_field_au = self.tesla_SI_to_au(self.b_sl_SI)*(1/self.m_to_au(1)) * x_au
        eff_B_field_SI = self.tesla_au_to_SI(eff_B_field_au)

        psi1 = eigenVectors[:, 0]
        psi2 = eigenVectors[:, 1]

        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.

        sigma_y_density = kwant.operator.Density(self.syst, self.sigma_y, sum=True).__call__
        sigma_y_overlap = sigma_y_density(psi2, psi1)

        # sigma_x_density = kwant.operator.Density(self.syst, self.sigma_x, sum=True).__call__
        # sigma_x_overlap = sigma_x_density(psi2, psi1)

        # sigma_z_density = kwant.operator.Density(self.syst, self.sigma_y, sum=True).__call__
        # sigma_z_overlap = sigma_z_density(psi2, psi1)

        time_steps = 25
        total_osc_time = 1/self.pulse_frequency_au

        times_au = np.linspace(0, total_osc_time, num=time_steps)
        B_SI = np.empty((time_steps,), dtype=np.double)
        B_au = np.empty((time_steps,), dtype=np.double)

        # potential_overlap_term_1 = -(1/2)* self.g * self.mu_B_au * self.B_0_au * sigma_z_overlap
        probabilities = []

        for i in range(time_steps):

            self.time = times_au[i]
            eigenValues, eigenVectors = self.eigenstates()
            psi_evolved = eigenVectors[:, 0]

            psi_evolved_up, psi_evolved_down = psi_evolved[::2], psi_evolved[1::2]
            density_1 = np.abs(psi_evolved_up) ** 2 + np.abs(psi_evolved_down) ** 2
            average_eff_B_field_au = np.trapz(density_1*eff_B_field_au, x=x_au)

            B_au[i] = average_eff_B_field_au


            probabilities.append(np.abs(np.conjugate(psi2).T @ psi_evolved)**2)

        # Plot magnetic field
        fig1 = plt.figure()
        plt.plot(self.time_au_to_SI(times_au), self.tesla_au_to_SI(B_au))
        plt.ylabel(r'$\bar{B}(t)$ (T)')
        plt.xlabel("$t$ (s)")
        # plt.title("Plot of Average Magnetic Field Varying with Time")
        plt.savefig("./figures/B-average-v-time.svg")
        plt.close(fig1)
        print("Magnetic field vs time plotted.")


        # # Plot magnetic field
        fig2 = plt.figure()
        plt.plot(self.time_au_to_SI(times_au), probabilities)
        plt.ylabel(r'$|c_2|^2$')
        plt.xlabel("$t$ (s)")
        # plt.title("Plot of Average Magnetic Field Varying with Time")
        plt.savefig("./figures/prob-v-time.svg")
        plt.close(fig2)
        print("prob vs time plotted.")



        return True

    def animate_wave_function(self):
        self.gaussian_mu = 0 # self.total_length_au / 4
        # create a figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1)
        x, y1,y2 = self.show_wave_function(animate=True)
        vpotential = vectorize(self.potential, otypes=[np.float64])
        y3 = vpotential(x)

        # intialize two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2, label='$n=1$')
        line2, = ax1.plot([], [], lw=2,label='$n=2$')
        line3, = ax2.plot([], [], lw=2, color='r', label='$V(x)$')
        line = [line1, line2, line3]
        xmax = (np.max(x))
        ax1.set_xlim(-xmax, xmax)
        ax1.set_ylim(0, np.max(y1) * 1.2)
        ax2.set_ylim(0, self.hartree_to_ev(np.max(y3) * 1.5))
        ax2.set_xlim(-xmax, xmax)

        plt.xlabel("z ($a_0$)")
        ax1.set_ylabel("$|\psi(z)|^2$")
        ax2.set_ylabel("$E$ (eV)")
        fig.suptitle("$t=0s$")
        number_of_frames = 120  # Takes [number_of_frames] time steps at frequency provided for pulse to travel along the
                               # the nanotube.

        # initialization function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[1].set_data([], [])
            line[2].set_data([], [])

            return line

        def animate(i):
            self.time =  (i / (self.pulse_frequency_au)) / number_of_frames
            fig.suptitle("$t={}s$".format('{:g}'.format(float('{:.{p}g}'.format(self.time_au_to_SI(self.time), p=2)))))

            x, y1,y2 = self.show_wave_function(animate=True)

            vpotential = vectorize(self.potential, otypes=[np.float64])

            y3 = vpotential(x)
            # x_nm = self.au_to_nm(x)
            line[0].set_data(x, y1)
            line[1].set_data(x, y2)
            line[2].set_data(x,  self.hartree_to_ev(y3))
            ax1.legend(loc="upper right")
            ax2.legend(loc="upper right")
            plt.savefig("./figures/temp/animation-{}.svg".format('{:g}'.format(float('{:.{p}g}'.format(self.time_au_to_SI(self.time), p=2)))))

            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate,init_func=init,
                                       frames=number_of_frames, interval=30, blit=True)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        if self.pulse_type == 1:
            anim.save('./figures/wavefunction-animation-pulse-1.mp4', writer='ffmpeg')
        elif self.pulse_type ==2:
            anim.save('./figures/wavefunction-animation-pulse-2.mp4', writer='ffmpeg')


        plt.close(fig)
        print("Animation of wave functions saved.")

def main():
    lattices = 100

    system1 = System("(C * k_x**2 + V(x))* identity(2) + A * sigma_x + B  * x * sigma_y", lattices)
    #+
    system1.make_system()
    system1.time = 0
    ### Phase 1 ###
    system1.show_energies()
    system1.show_wave_function()


    ### Phase 2 ###
    # system1.animate_wave_function()

    ### Phase 3 ###
    system1.rabi_oscillations()

if __name__ == '__main__':
    main()

