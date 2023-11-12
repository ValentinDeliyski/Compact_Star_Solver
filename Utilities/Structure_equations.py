class TOV:

    def __init__(self, E_density_core: float, P_core: float, P_surface: float, EOS_class, Units_class_instance):

        self.E_0   = E_density_core
        self.P_0   = P_core
        self.Units = Units_class_instance
        self.P_surface = P_surface
        self.EOS   = EOS_class


    def ODE_System(self, State_vector, Independant_var) -> list:

        """
        |========================================================================|
        | Returns the right hand side of the Toloman - Oppenheimer - Volkov      |
        | ODE system as a list. The input parameters are scaled as follows:      |
        |                                                                        |
        |  - Pressure  is dimentionless - in multiples of P_0 [MeV / fm^3]       |
        |  - E_density is dimentionless - in multiples of E_0 [MeV / fm^3]       |
        |  - Mass is dimentionless - in multiples of solar masses                |
        |  - r is dimentionless - in multiples of half Solar Schwarzschild radii |
        |========================================================================|
        
        """

        r = Independant_var

        Pressure = State_vector[0]
        mass     = State_vector[1]

        E_density = self.EOS.get_E_of_P(Pressure_MeV_fm3 = Pressure * self.P_0) / self.E_0

        dP_dr     = -1 /  r**2 * ( (self.E_0 / self.P_0) * E_density + Pressure ) * ( (self.P_0 / self.Units.E_density_norm_MeV_fm3) * Pressure * r**3 + mass ) / ( 1 - 2 * mass / r )
        dm_dr     = ( self.E_0 / self.Units.E_density_norm_MeV_fm3 ) * r**2 * E_density
        dGamma_dr = -1 / ( (self.E_0 / self.P_0) * E_density + Pressure ) * dP_dr

        return [dP_dr, dm_dr, dGamma_dr] 
