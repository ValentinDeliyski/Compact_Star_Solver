class Units:

    def __init__(self):

        self.FEMPTO = 1e-15

        self.KILO = 1e3
        self.MEGA = 1e6

        self.g_cm3_to_kg_m3 = 1000

        self.e_charge_coulumb = 1.602e-19
        self.m_electron_kg    = 9.1094e-31

        self.G_Newton   = 6.67e-11
        self.c_light_ms = 2.99792458e8
        self.h_bar      = 1.05457187e-34

        self.M_sun_kg          = 1.9e30
        self.Sun_Sch_radius_m  = self.G_Newton * self.M_sun_kg / self.c_light_ms**2
        self.Sun_Sch_radius_km = self.Sun_Sch_radius_m / self.KILO

        self.J_to_ev = 1 / self.e_charge_coulumb

        # So I dont have to import a whole module for one number...
        self.pi = 3.141592653589793238462643383279502884197

        self.E_density_norm_Jm3     = self.M_sun_kg * self.c_light_ms**2 / (4 * self.pi * self.Sun_Sch_radius_m**3)
        self.E_density_norm_MeV_fm3 = self.Jm3_to_MeVfm3(self.E_density_norm_Jm3)

    def Jm3_to_MeVfm3(self, E_density_Jm3):

        return E_density_Jm3 * (self.J_to_ev / self.MEGA) / (1 / self.FEMPTO)**3
    
    def MeVfm3_to_Jm3(self, E_density_Mevfm3):

        return E_density_Mevfm3 * (1 / self.J_to_ev * self.MEGA) / (self.FEMPTO)**3
    
