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
    def __init__(self, hamiltonian, pertubation, number_of_lattices=50, potential_type = 0):

        # Define units
        self.potential_type = potential_type
        self.hamiltonian = hamiltonian
        self.number_of_lattices = number_of_lattices

        # Constants in SI units
        self.hbar_SI = 1.054571817e-34
        self.e_SI = 1.602176634e-19
        self.a_0_SI = 5.2917721090380e-11
        self.total_length_SI = 0.66e-6
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.b_sl_SI = 1.16e6
        self.m_SI = 9.11e-31  # divide by ten to get the effective mass
        self.mu_B_SI= 9.2740100783e-24
        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.number_of_lattices)

        # Constants in a.u.
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)
        self.b_sl_au =self.tesla_SI_to_au(self.b_sl_SI) * (1/self.m_to_au(1))
        self.hbar_au = 1
        self.g = 2
        self.m_au = 1
        self.e_au = 1
        self.total_length_au = self.total_length_SI / self.a_0_SI # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices # The distance in atomic units spanned by a lattice point
        self.mu_B_au = .5
        self.a_0 = 1
        self.z_au =  np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.number_of_lattices)
        self.pulse_frequency_au =  (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au) # in Hz
        self.pertubation = pertubation


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

    def import_mumax3_simulations(self):

        file_name = "simulation-0-x"

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


        return B_x_mean
        # fig = plt.figure()
        #
        #
        # plt.plot(y_coords, B_x_mean, label="$B_x$")
        # # plt.errorbar(x_coords, B_x_mean, yerr=B_x_sd)
        # plt.plot(y_coords, B_y_mean, label="$B_y$")
        # # plt.errorbar(x_coords, B_y_mean, yerr=B_y_sd)
        # plt.plot(y_coords, B_z_mean, label="$B_z$")
        # # plt.errorbar(x_coords, B_z_mean, yerr=B_z_sd)
        # plt.ylabel("$B_{eff}$ $(T)$")
        # plt.xlabel("$y$ $(10^{-8}m)$")
        # plt.legend()
        # plt.savefig("./figures/{}.svg".format(file_name))
        # plt.close(fig)

    def infinite_square_well_potential(self):

        self.total_length_SI = 0.66e-6
        self.E_sl_ev = 1e-6
        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.number_of_lattices)
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.B_0_au = self.tesla_SI_to_au(self.B_0_SI)


        self.total_length_au = self.total_length_SI / self.a_0_SI  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # The distance in atomic units spanned by a lattice point

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.total_length_au) # slanted magnetic field value computed from E_sl
        self.eV_0_au = self.ev_to_hartree(10e-6) # the value in a.u. of eV_0.



    def parabolic_potential(self):
        self.total_length_SI = 0.66e-6
        self.confinement_length_SI = 0.33e-6
        self.confinement_length_au = self.m_to_au(self.confinement_length_SI)

        self.z_SI =  np.linspace(-self.total_length_SI / 2, self.total_length_SI / 2, self.number_of_lattices)
        self.B_0_SI = 5e-3 # Upon using > 250e-3 the two wavefunctions represent the n=1 and n=2 states
        self.E_sl_ev = 1e-6

        self.total_length_au = self.total_length_SI / self.a_0_SI  # Total distance of nanotube in terms of au
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # The distance in atomic units spanned by a lattice point
        self.z_au = np.linspace(-self.total_length_au / 2, self.total_length_au / 2, self.number_of_lattices)

        self.E_sl_au = self.ev_to_hartree(self.E_sl_ev)
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.confinement_length_au) # slanted magnetic field value computed from E_sl
        self.omega_0 = self.hbar_au/(self.m_au*(self.confinement_length_au)**2)
        self.eV_0_au = self.ev_to_hartree(10e-6) # the value in a.u. of eV_0.

        print("b_sl is", self.tesla_au_to_SI(self.b_sl_au) / self.au_to_m(1), "T.")
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
            total_potential = .5 * (
                    (z * self.omega_0) ** 2)  # define a parabolic potential inside the scattering region.

        total_potential += self.pertubation(time, z, self.eV_0_au, self.pulse_frequency_au, self.total_length_au)

        return total_potential

    def kwant_shape(self, site):
        """
        function to define the shape of the scattering region.
        :param site: the current site.
        :return: the a boolean saying whether the scattering site should be drawn
        """
        (z, ) = site.pos
        return (-self.total_length_au/2 <= z  < self.total_length_au/2)

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





        self.syst = kwant.Builder()

        #Add the nanotube to the system
        self.syst.fill(self.template, self.kwant_shape, (0, ))

        # kwant.plot(self.syst, file='./figures/shape.png')
        self.syst = self.syst.finalized()

        self.A_constant = -self.g * self.mu_B_au * self.B_0_au * self.hbar_au / 2
        self.B_constant = -self.g * self.mu_B_au * self.b_sl_au * self.hbar_au / 2
        self.C_constant = self.hbar_au ** 2 / (2 * self.m_au)

        self.params = dict(C=self.C_constant, V=self.potential, A=self.A_constant, B=self.B_constant)
        self.tparams = self.params.copy()
        self.tparams['time'] = 0  # the initial time
        print("System intialised.")

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


        # perturb = tkwant.onebody.kernels.PerturbationExtractor(self.syst, time_name='time', time_start=0)#, params=self.params)
        self.spin_up_state = tkwant.onebody.WaveFunction.from_kwant(syst=self.syst,
                                                              psi_init=eigenVectors[:,0],
                                                              energy=eigenValues[0],
                                                              params=self.params)
        #
        # pickle.dump((self.spin_up_state),
        #             open('state.npy', "wb")) # doesn't work

        return eigenValues, eigenVectors




    def evolve(self, time_steps = 100):
        timestr = timelib.strftime("%Y%m%d-%H%M%S")


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


        total_osc_time = self.t_pi

        # compute oscillation times.
        times_au = np.linspace(0, total_osc_time, num=time_steps)
        times_SI = np.linspace(0, self.time_au_to_SI(total_osc_time), num=time_steps)

        B_au = np.empty((time_steps,), dtype=np.double)

        # compute the spin operators at new time
        rho_sz = kwant.operator.Density(self.syst, sigma_z, sum=True)
        rho_sy = kwant.operator.Density(self.syst, sigma_y, sum=True)
        rho_sx = kwant.operator.Density(self.syst, sigma_x, sum=True)

        density_operator = kwant.operator.Density(self.syst)
        psi = self.spin_up_state
        time = times_au[0]

        data = {
            'B_0':self.B_0_au,
            'length': self.total_length_au,
            'V_0': self.eV_0_au,
            'E_sl': self.E_sl_au,
            'E_x': E_x,
            'E_z': omega_res,
            'states': []
        }

        print("Simulation starting with", time_steps, "time steps from 0.0 to",total_osc_time,"a.u.")

        # for time in times_au:
        time = times_au[0]
        print("Evolving state to time",time,"a.u.")
        psi.evolve(time)
        density = np.abs(self.spin_up_state.evaluate(density_operator))**2
        average_eff_B_field_au = np.trapz(density*eff_B_field_au, x=self.z_au) #compute the expectation of the slanted field.

        # compute the expecrations of the spin operators on the evolved state and store their values.
        spin_z = np.real(self.spin_up_state.evaluate(rho_sz))
        spin_y =  np.real(self.spin_up_state.evaluate(rho_sy))
        spin_x =  np.real(self.spin_up_state.evaluate(rho_sx))


        data['states'].append({
            'time': time,
            'pdf': density.tolist(),
            'B_x': average_eff_B_field_au,
            'rho_sx': spin_x,
            'rho_sy': spin_y,
            'rho_sy': spin_z,
        })

        json_string = json.dumps(data)

        with open('{}.json'.format(timestr), 'w') as outfile:
            outfile.write(json_string)

def main():

    lattices = 100 # number of lattice points
    def cosine_V_ac(time, z, eV_0_au, pulse_frequency_au, total_length_au):
        return ((eV_0_au * np.cos(2 * np.pi * pulse_frequency_au * time)) * z) / total_length_au
    potential = 0 # infinite square-well potential
    system = System("((C * k_z**2) + V(z, time)) * identity(2) + A * sigma_z + B  * z * sigma_x", pertubation = cosine_V_ac, number_of_lattices=lattices, potential_type=potential) # call system objecy
    system.make_system()
    # system.evolve(100)
    system.import_mumax3_simulations()


if __name__ == '__main__':
    main()

