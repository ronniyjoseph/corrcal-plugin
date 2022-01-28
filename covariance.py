import numpy as np
from matplotlib import pyplot as plt

from skymodel import sky_moment
from telescope import beam_width



class CovarianceMatrix:

    def __init__(self, u, v, nu):
        self.u = u
        self.v = v
        self.nu = nu

        self.beam_model =  None
        self.baselines =  None
        self.matrix =  None
        self.eigenvalues = None
        self.eigenmodes =  None

        return

    def gaussian(self, aperture_diameter=2, S_low=1e-3, S_mid=1, S_high=10, gamma=0.8, nu_0 = 150e6):
        """
        Flat-sky gaussian approximation of a antenna beam

        """
        #Compute 2nd moment of the Point Source counts function thingy

        #So this point source sky model is calibrated for 150 MHz, can we scale this?
        u, v, nu =  self.return_grid()
        mu = sky_moment(2, s_low=S_low, s_mid=S_mid, s_high=S_high)


        #If we're computing the covariance matrix over multiple frequencies we still want to end up with a 2D array that
        #that has the baseline-frequency covariance

        beam_width1 = beam_width(nu[0], diameter=aperture_diameter)
        beam_width2 = beam_width(nu[1], diameter=aperture_diameter)

        sigma_nu = beam_width1 ** 2 * beam_width2 ** 2 / (beam_width1 ** 2 + beam_width2 ** 2)

        kernel = -2 * np.pi ** 2 * sigma_nu * ( (u[0] * nu[0] - u[1] * nu[1]) ** 2 +
                                                (v[0] * nu[0] - v[1] * nu[1]) ** 2) / nu_0 ** 2
        covariance = 2 * np.pi * mu * sigma_nu * (nu[0] * nu[1] / nu_0 ** 2) ** (- gamma) * np.exp(kernel)

        plt.imshow(covariance)
        plt.show()
        return

    def return_grid(self):
        uu1, uu2 = np.meshgrid(self.u, self.u)
        vv1, vv2 = np.meshgrid(self.v, self.v)

        u_grid1 = np.tile(uu1, np.array([len(self.nu), len(self.nu)]))
        u_grid2 = np.tile(uu2, np.array([len(self.nu), len(self.nu)]))

        v_grid1 = np.tile(vv1, np.array([len(self.nu), len(self.nu)]))
        v_grid2 = np.tile(vv2, np.array([len(self.nu), len(self.nu)]))

        nu_grid1 = np.zeros_like(u_grid1)
        nu_grid2 = np.zeros_like(u_grid1)

        # Create the full covariance_coordinate grid?
        for i in range(len(self.nu)):
            nu_grid1[i * len(self.u):(i + 1) * len(self.u), :] = self.nu[i]
            nu_grid2[:, i * len(self.u):(i + 1) * len(self.u)] = self.nu[i]

        return (u_grid1, u_grid2), (v_grid1, u_grid2), (nu_grid1, nu_grid2)



    def airy(self):
        """
        Flat-sky airy beam approximation of an antenna beam
        """
        return

    def actualbeam(self):
        return

    def compute(self):
        return

    def decompose(self):
        return



