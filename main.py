import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
import tkwant
import scipy.sparse.linalg as sla
from numpy import vectorize
from matplotlib import animation

class System:
    def __init__(self, hamiltonian, lattices):
        # constants in SI
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_m = 7820e-9
        self.B_0_SI = 100e-3
        self.b_sl_SI = 250e-3
        self.time = 0
        # Constants in a.u.
        self.hbar = 1
        self.g = 2
        self.m = 1
        self.e = 1
        self.total_length_au = self.total_length_m / self.a_0_SI
        self.lattices = lattices
        self.lattice_size = self.total_length_au / self.lattices
        self.mu_B = self.e * self.hbar / (2 * self.m)
        self.a_0 = 1
        self.pulse_velocity = 0.1
        self.hamiltonian = hamiltonian
        self.template = kwant.continuum.discretize(hamiltonian)
        self.gaussian_mu = self.lattices/4
        self.a = 1

    def tesla_SI_to_au(self, tesla):
        """
        Function to convert the magnetic flux density from SI units to AU units
        :param tesla: the magnetic flux density in teslas.
        :return: the magnetic flux density in AU.
        """
        return tesla / (self.hbar_SI/(self.e_SI * (self.a_0_SI ** 2)))







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
            return (0 <= x < self.lattices)


        def lead_shape(site):
            """
            function to define the shape of the leads.
            :param site: the current site.
            :return: the a boolean saying whether the lead site should be drawn
            """
            (x, ) = site.pos
            return (0 <= x < self.lattices)

        self.syst = kwant.Builder()

        #Add the nanotube to the system
        self.syst.fill(self.template, shape, (0, ));

        # #Attach the left gate to the system - for now we assume it's parallel to the nanotube
        # sym_left_lead = kwant.TranslationalSymmetry((-self.a, ))
        # left_lead = kwant.Builder(sym_left_lead)
        # left_lead.fill(self.template, lead_shape, (-self.lattices,))
        # self.syst.attach_lead(left_lead)
        #
        #
        # #Attach the right gate to the system - for now we assume it's parallel to the nanotube
        # sym_right_lead = kwant.TranslationalSymmetry((self.a, ))
        # right_lead = kwant.Builder(sym_right_lead)
        # right_lead.fill(self.template, lead_shape, (self.lattices+self.a,)) # The leads have the same hamiltonian as the scattering regions
        # self.syst.attach_lead(right_lead)

        kwant.plot(self.syst, file='./figures/shape.png');
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

        return np.exp(-np.power(x - self.gaussian_mu - self.pulse_velocity * self.time, 2.) / (2 * np.power(sig, 2.)))


    def potential(self, x):  # Infinite square well
        """
        Function to define the potential of the lead.
        :param x: the position in the system.
        :return: the potential energy.
        """
        if 0 <= x <= self.lattices:
            return 0.002 * self.gaussian(x, self.lattices / 14)
        else:
            return 999999999


    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        :param syst: the system object.
        :return: the sorted eigenvalues and eigenvectors.
        """
        B_0_au = self.tesla_SI_to_au(self.B_0_SI)
        b_sl_au = self.tesla_SI_to_au(self.b_sl_SI)  # in hbar/(e*(a_0)**2)
        self.A_constant =  -self.g * self.mu_B * B_0_au * self.hbar / 2
        self.B_constant = -self.g * self.mu_B * b_sl_au * self.hbar / 2
        self.C_constant = self.hbar ** 2 / (2 * self.m)
        params = dict(A=self.A_constant, B=self.B_constant, C=self.C_constant, V=self.potential)
        # Calculate the wave functions in the system.
        h = self.syst.hamiltonian_submatrix(params=params)

        eigenValues, eigenVectors = np.linalg.eig(h)

        idx = eigenValues.argsort()
        eigenValues = np.real(eigenValues[idx])
        eigenVectors = eigenVectors[:, idx]

        # print(eigenValues) # Eigenvalues are repeated
        # For each eigenvalue for each spin we get two eigenvectors and for each position there are two possible spins
        # so we have four eigenvectors for each position along the nanotube.
        # print(np.abs(eigenVectors[0::2,1])**2 + np.abs(eigenVectors[1::2,1])**2 +
        #       np.abs(eigenVectors[0::2,2])**2+ np.abs(eigenVectors[1::2,2])**2)

        return eigenValues, eigenVectors

    def showEnergies(self):
        """
        Procedure to display the potential and energy levels of the system
        :param syst: the system object.
        :return:
        """
        eigenValues, eigenVectors = self.eigenstates()
        x_coordinates = np.linspace(0, self.lattices, self.lattices)
        vpotential = vectorize(self.potential)

        y_coordinates = vpotential(x_coordinates)
        plt.figure()
        plt.plot(x_coordinates, y_coordinates, label="$V(x)$")
        E_1 = eigenValues[0]
        E_2 = eigenValues[1]
        E_3 = eigenValues[2]

        # print("The energies of the ground and excited states are {0} and {1}, respectively.".format(E_1, E_2))
        plt.plot([0,100], [E_1, E_1], label="$E_1$")
        plt.plot([0,100], [E_2, E_2], label="$E_2$")
        plt.plot([0,100], [E_3, E_3], label="$E_3$", linestyle='dashed')
        plt.xlabel("$x (au$)")
        plt.ylabel("$E (H)$")
        plt.legend(loc="upper right")
        plt.savefig("./figures/Energies/energies-{}.svg".format(self.time))  # With A = 0 we expect straight forward zeeman splitting
        plt.close()
        print("Plot of energies saved.")

    def showWavefunction(self, animate=False):

        """
        Procedure to show the probability density function.
        :param syst: the system object.
        :return:
        """

        eigenValues, eigenVectors = self.eigenstates()


        x_coords = np.linspace(0, self.total_length_au, self.lattices)
        psi1 = np.abs(eigenVectors[0::2,0])**2 + np.abs(eigenVectors[1::2,0])**2 +\
               np.abs(eigenVectors[0::2,1])**2+ np.abs(eigenVectors[1::2, 1])**2

        if animate:
            return x_coords,psi1

        else:
            plt.figure()

            plt.plot(x_coords, psi1, label="n=1")

            psi2 = np.abs(eigenVectors[0::2,2])**2 + np.abs(eigenVectors[1::2,2])**2 +\
                   np.abs(eigenVectors[0::2,3])**2+ np.abs(eigenVectors[1::2,3])**2

            plt.plot(x_coords, psi2, label="n=2")



            plt.xlabel("x (au)")
            plt.ylabel("$|\psi(x)|^2$")
            plt.legend(loc="upper right")
            plt.savefig("./figures/PDFs/pdf-{}.svg".format(self.time))  # With A = 0 we expect straight forward zeeman splitting
            plt.close()
            print("Plot of wave functions saved!")
            return True



    def animateWavefunction(self):

        fig = plt.figure()
        x, y = self.showWavefunction(animate=True)

        ax = plt.axes(xlim=(np.min(x), np.max(x)), ylim=(np.min(y), np.max(y) * 1.2))

        line, = ax.plot([], [], lw=2)
        plt.xlabel("x (au)")
        plt.ylabel("$|\psi(x)|^2$")
        plt.title("Animation of $n = 1$ Wave Function")

        self.gaussian_mu = 0

        # initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        def animate(i):
            self.time = i
            line.set_label('$t = {}$'.format(i))
            plt.legend(loc="upper right")

            x, y = self.showWavefunction(animate=True)
            line.set_data(x, y)

            return line,

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=1000, interval=10, blit=True)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        anim.save('wavefunction-animation.mp4', writer='ffmpeg')
        plt.close()
        print("Animation of wave function saved!")

def main():
    lattices = 100

    system1 = System("(C * k_x**2 + V(x)) * identity(2) + A * sigma_x + B * x * sigma_y", lattices)
    system1.make_system()
    system1.animateWavefunction()
    system1.time = 0
    system1.showEnergies()

    system1.showWavefunction()


if __name__ == '__main__':
    main()

# Output magnetic field profiles along the wire (mT) from mumax simulation

# Add pulse gaussian along the wire depending on time
# Compute the probability of occupation of state 2 as time varies
# Display energy levels, display potential energy
# Make code OOP


# Fix scaling to give it in-terms of eV -> give in terms of a.u.
