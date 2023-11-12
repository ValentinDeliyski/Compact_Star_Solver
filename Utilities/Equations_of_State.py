from numpy import power, array, append, zeros
from Utilities.Units import Units

class Single_Polytrope:

    def __init__(self):

        self.K_rel     = 4.012e-4
        self.gamma_rel = 2.0

    def get_P_of_E(self, E_density_MeV_fm3):

        return self.K_rel * power(E_density_MeV_fm3, self.gamma_rel)
    
    def get_E_of_P(self, Pressure_MeV_fm3):
        
        if Pressure_MeV_fm3 < 0:

            Pressure_MeV_fm3 = 0

        return  power((Pressure_MeV_fm3 / self.K_rel), (1.0 / self.gamma_rel))
    
class Multi_Polytrope():

    def __init__(self):

        self.Units = Units()

        #============= Load magic number EOS coeffs from Daniela/Pesho =============#

        # TODO: Figure out where these come from...
        # These convert energy density in [J / m^3] to P in [dyne / cm^2]
        self.K = array([6.80110e-9, 1.06186e-6, 5.32697e+1, 3.99874e-8]) * self.Units.c_light_ms**2

        # Here I convert K to spit out P in [MeV / fm^3]
        self.K = self.Units.Jm3_to_MeVfm3(self.K / 10)

        # TODO: Figure out where these come from...
        density_cutoff_SI  = array([2.44034e7, 3.78358e11, 2.62780e12, 0.0]) * self.Units.g_cm3_to_kg_m3

        density_cutoff_4_SI = pow(10, 14.7) * self.Units.g_cm3_to_kg_m3
        density_cutoff_5_SI = pow(10, 15)   * self.Units.g_cm3_to_kg_m3
        
        self.density_cutoff_SI = append(density_cutoff_SI, [density_cutoff_4_SI, density_cutoff_5_SI, 1e99], axis = 0)

        # TODO: Figure out where these come from...
        self.Gamma = array([1.58425, 1.28733, 0.62223, 1.35692])
        
        #============================================================================#

        # The pressure at the start of the first polytrope interval... in silly units - 
        # TODO: Should eventually be read from a file...
        p_1 = pow(10, 34.495) # [ dyne / cm^2 ]... or [ 10 * Pa ]
        p_1_MeV_fm3 = self.Units.Jm3_to_MeVfm3(E_density_Jm3 = p_1 / 10)

        # TODO: This should eventually be read from a file...
        self.Gamma = append(self.Gamma, [3.446, 3.572, 2.887], axis=0)

        # Calculate K[4], based on the fit parameter logp1
        self.K = append(self.K, [p_1_MeV_fm3 / pow(self.density_cutoff_SI[4], self.Gamma[4])], axis = 0)

        # Calculate density_cutoff_SI[3], based on the fit coeffs at indecies 3 and 4
        try:
            self.density_cutoff_SI[3] = pow(self.K[3] / self.K[4], 1.0 / (self.Gamma[4] - self.Gamma[3]))

        except:
            self.density_cutoff_SI[3] = 1e14 * self.Units.g_cm3_to_kg_m3

        # Calculate the remaining two K coefficients of the polytropes
        for index in [5, 6]:
            self.K = append(self.K, [self.K[index - 1] * pow(self.density_cutoff_SI[index - 1], self.Gamma[index - 1] - self.Gamma[index])], axis=0)

        # for index in [3, 2, 1, 0]:
        #     self.K[index] = self.K[index + 1] / pow(self.density_cutoff_SI[index], self.Gamma[index] - self.Gamma[index + 1])


        # Calculate the a coefficients of the polytropes and the pressure cutoff values
        self.a_coeff = zeros(7)
        self.pressure_cutoff = zeros(7)

        for index in range(7):

            if index != 0:

                last_Energy_Density_MeVfm3 = self.Units.Jm3_to_MeVfm3(self.density_cutoff_SI[index - 1] * self.Units.c_light_ms**2)

                self.a_coeff[index] = (self.a_coeff[index-1]
                                    + self.K[index - 1] / last_Energy_Density_MeVfm3 / (self.Gamma[index - 1] - 1) * pow(self.density_cutoff_SI[index - 1], self.Gamma[index - 1])
                                    - self.K[  index  ] / last_Energy_Density_MeVfm3 / (self.Gamma[  index  ] - 1) * pow(self.density_cutoff_SI[index - 1], self.Gamma[  index  ]))

            self.pressure_cutoff[index] = self.K[index] * pow(self.density_cutoff_SI[index], self.Gamma[index])

    def get_rho_of_P(self, P_MeVfm3):

        Poly_index = 6

        for index, pressure_cutoff in enumerate(self.pressure_cutoff):

            if P_MeVfm3 < pressure_cutoff:
                Poly_index = index

                break

        # if Poly_index < 4:
        #     Poly_index = 4

        return pow(P_MeVfm3 / self.K[Poly_index], 1 / self.Gamma[Poly_index])


    def get_E_of_rho(self, rho_SI):

        Poly_index = 6

        for index, density_cutoff in enumerate(self.density_cutoff_SI):

            if rho_SI < density_cutoff:
                Poly_index = index

                break

        # if Poly_index < 4:
        #     Poly_index = 4

        return (1 + self.a_coeff[Poly_index]) * self.Units.Jm3_to_MeVfm3(rho_SI * self.Units.c_light_ms**2) + self.K[Poly_index] / (self.Gamma[Poly_index] - 1) * pow(rho_SI, self.Gamma[Poly_index])

    def get_E_of_P(self, Pressure_MeV_fm3):

        if Pressure_MeV_fm3 < 0:
            Pressure_MeV_fm3 = 0

        mass_density_SI = self.get_rho_of_P(P_MeVfm3 = Pressure_MeV_fm3)

        return self.get_E_of_rho(rho_SI = mass_density_SI)
