from Utilities.Equations_of_State import Single_Polytrope, Multi_Polytrope
from Utilities.Units import Units
from Utilities.Integrator import Dormand_Prince
from Utilities.Structure_equations import TOV
from numpy import tanh

import matplotlib.pyplot as plt

Units_class = Units()
EOS_class   = Multi_Polytrope()

Star_Number = 200
P_surface   = 1e-10

R_star = []
M_star = []

def f_interpolant(x, scale_param):

    return (tanh(scale_param * (x - 1)) - tanh(-scale_param)) / (-2 * tanh(-scale_param))

for core_pressure_count in range(Star_Number + 1):

    """
    |=========================== Setup Initial Conditions ===========================|

    """

    P_core = 5 + 5000 * f_interpolant(2 * core_pressure_count / Star_Number, 3.6) # [ MeV / fm^3 ]
    E_density_core = EOS_class.get_E_of_P(Pressure_MeV_fm3 = P_core)

    Init_Radial_Coord = 1e-10  # [ R_half_Schw_Sun [ m ] ]

    # The total mass in the sphere with radius r_init = Init_Radial_Coord * Sun_half_Sch_radius_m
    M_0 = 4 / 3 * Units_class.pi * (Init_Radial_Coord * Units_class.Sun_Sch_radius_m) ** 3 * E_density_core / Units_class.c_light_ms**2 / Units_class.M_sun_kg

    # The initial State Vector
    Init_State_Vector = [1, M_0, 0]

    # Init the ODE system class
    ODE_System = TOV(E_density_core, P_core, P_surface, EOS_class, Units_class)

    # Init the Integrator class
    Integrator = Dormand_Prince(init_step        = 1e-6, 
                                desired_accuracy = 1e-16, 
                                ODE_System  = ODE_System, 
                                Units_class = Units_class)

    """
    |================================================================================|

    """

    _, _, R, M = Integrator.run_integrator(Integrator_radial_stop_point = None, # [ km ], if None - stops at P(r) = P_surface
                                           Initial_State                = Init_State_Vector,
                                           Init_Independant_Var         = Init_Radial_Coord)

    R_star.append(R)
    M_star.append(M)

plt.figure()
plt.plot(R_star, M_star)
plt.xlabel("Star Radius [km]")
plt.ylabel("Star Mass [Solar Msses]")
plt.show()

