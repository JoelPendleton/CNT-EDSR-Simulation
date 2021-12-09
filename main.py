import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
import tkwant
import scipy.sparse.linalg as sla
from numpy import vectorize

hamiltonian = "(C * k_x**2 + V(x)) * identity(2) + A * sigma_x + B * x * sigma_y"
template = kwant.continuum.discretize(hamiltonian)

# constants in SI
hbar_SI = 1.054571817e-34
e_SI = 1.602176634e-19
a_0_SI = 5.2917721090380e-11
total_length_m = 7820e-9
B_0_SI = 100e-3
b_sl_SI = 150e-3

# Constants in a.u.
hbar = 1
g = 2
m = 1
e = 1
total_length_au = total_length_m / a_0_SI
lattices = 100
lattice_size = total_length_au / lattices
mu_B = e * hbar / (2 * m)
a_0 = 1

def tesla_SI_to_au(tesla):
    """
    Function to convert the magnetic flux density from SI units to AU units
    :param tesla: the magnetic flux density in teslas.
    :return: the magnetic flux density in AU.
    """
    return tesla / (hbar_SI/(e_SI * (a_0_SI ** 2)))

B_0_au = tesla_SI_to_au(B_0_SI)
b_sl_au = tesla_SI_to_au(b_sl_SI) # in hbar/(e*(a_0)**2)

A_constant = -g * mu_B * B_0_au * hbar/2
B_constant = -g * mu_B * b_sl_au * hbar/2
C_constant = hbar**2 / (2*m)
a = 1


def make_system(length):
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
        return (0 <= x < length)


    def lead_shape(site):
        """
        function to define the shape of the leads.
        :param site: the current site.
        :return: the a boolean saying whether the lead site should be drawn
        """
        (x, ) = site.pos
        return (0 <= x < length)

    syst = kwant.Builder()

    #Add the nanotube to the system
    syst.fill(template, shape, (0, ));


    #Attach the left gate to the system - for now we assume it's parallel to the nanotube
    sym_left_lead = kwant.TranslationalSymmetry((-a, ))
    left_lead = kwant.Builder(sym_left_lead)
    left_lead.fill(template, lead_shape, (-length,))
    syst.attach_lead(left_lead)


    #Attach the right gate to the system - for now we assume it's parallel to the nanotube
    sym_right_lead = kwant.TranslationalSymmetry((a, ))
    right_lead = kwant.Builder(sym_right_lead)
    right_lead.fill(template, lead_shape, (length+a,)) # The leads have the same hamiltonian as the scattering regions
    syst.attach_lead(right_lead)

    return syst

def sorted_eigs(ev):
    """
    Function to sort eigenvalues and vectors so they're in ascending order.
    :param ev: the numpy array contianing the eigenvalues and their corresponding eigenvectors.
    :return: the sorted eigenvalues and sorted eigenvectors.
    """
    evals, evecs = ev
    evals, evecs = map(np.array, zip(*sorted(zip(evals, evecs.transpose()), key=lambda x: x[0])))
    return evals, evecs.transpose()

def gaussian(x, mu, sig):
    """
    Function to compute a gaussian pulse.
    :param x: the coordinate value
    :param mu: the centre of the gaussian
    :param sig: the standard deviation
    :return: the gaussian function values.
    """
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def potential(x):  # Infinite square well
    """
    Function to define the potential of the lead.
    :param x: the position in the system.
    :return: the potential energy.
    """
    if 0 <= x <= lattices / 2:
        return 0.001 * gaussian(x, lattices / 4, lattices / 14)
    elif lattices / 2 <= x <= lattices:
        return 0
    else:
        return 999999999


def eigenstates(syst):
    """
    Function to compute the eigenstates of the system.
    :param syst: the system object.
    :return: the sorted eigenvalues and eigenvectors.
    """
    params = dict(A=A_constant, B=B_constant, C=C_constant, V=potential)
    # Calculate the wave functions in the system.
    ham_mat = syst.hamiltonian_submatrix(sparse=True, params=params)
    return sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=2, sigma=0))

def displayPotential(syst):
    """
    Procedure to display the potential and energy levels of the system
    :param syst: the system object.
    :return:
    """
    evals, evecs = eigenstates(syst)
    x_coordinates = np.linspace(0, lattices, lattices)
    vpotential = vectorize(potential)

    y_coordinates = vpotential(x_coordinates)
    plt.figure()
    plt.plot(x_coordinates, y_coordinates, label="$V(x)$")
    E_1 = evals[0]
    E_2 = evals[1]
    plt.plot([0,100], [E_1, E_1], label="$E_1$")
    plt.plot([0,100], [E_2, E_2], label="$E_2$")
    plt.legend(loc="upper right")
    plt.savefig("energies.png")  # With A = 0 we expect straight forward zeeman splitting

def showWavefunction(syst):

    """
    Procedure to show the wave functions.
    :param syst: the system object.
    :return:
    """

    evals, evecs = eigenstates(syst)
    plt.figure()
    # print(evals) # Not sure why but eigenvalues and vectors are repeated

    x_coords = np.linspace(0, total_length_au, lattices)


    psi1 = np.abs(evecs[0::2, 0])**2 + np.abs(evecs[1::2, 0])**2  # Two eigenvectors per eigenvalue
    psi2 = np.abs(evecs[0::2, 2])**2 + np.abs(evecs[1::2, 2])**2

    plt.plot(x_coords, psi1, label="n=1")
    plt.plot(x_coords, psi2, label="n=2")

    plt.xlabel("x (au)")
    plt.ylabel("$|\psi(x)|^2$")
    plt.legend(loc="upper right")
    plt.savefig("wavefunctions.png")  # With A = 0 we expect straight forward zeeman splitting

    # Next add alternating voltage


def main():
    syst = make_system(length=lattices).finalized()

    kwant.plot(syst, file='shape.png');
    displayPotential(syst)
    # showWavefunction(syst, lattices)

if __name__ == '__main__':
    main()

# Output magnetic field profiles along the wire (mT) from mumax simulation

# Add pulse gaussian along the wire depending on time
# Compute the probability of occupation of state 2 as time varies
# Display energy levels, display potential energy
# Make code OOP


# Fix scaling to give it in-terms of eV -> give in terms of a.u.
