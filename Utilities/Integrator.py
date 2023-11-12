from numpy import power, dot, max, absolute, zeros

class Dormand_Prince():

    def __init__(self, init_step: float, desired_accuracy: float, ODE_System, Units_class):

        
        self.Coeff_deriv = [[       0,               0,              0,             0,            0,           0   ],
                            [    1. / 5,             0,              0,             0,            0,           0   ],
                            [    3. / 40,         9. / 40,           0,             0,            0,           0   ],
                            [   44. / 45,       -56. / 15,       32. / 9,           0,            0,           0   ],
                            [19372. / 6561,  -25360. / 2187,  64448. / 6561,  -212. / 729,        0,           0   ],
                            [ 9017. / 3168,    -355. / 33,    46732. / 5247,    49. / 176, -5103. / 18656,     0   ],
                            [   35. / 384,           0,         500. / 1113,   125. / 192, -2187. / 6784,  11. / 84]]

        self.Coeff_ind_var  = [       0,      1. / 5,    3. / 10,      4. / 5,        8. / 9,           1,         1    ]
        self.Coeff_sol      = [   35. / 384,     0,    500. / 1113,  125. / 192,  -2187. / 6784,    11. / 84,      0    ]
        self.Coeff_test_sol = [ 5179. / 57600,   0,   7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 ]

        self.RK45_size = 7

        self.step     = init_step
        self.accuracy = desired_accuracy
        self.Safety_1 = 0.8
        self.Safety_2 = 1e-20
        self.Pressure_at_surface = ODE_System.P_surface

        self.contuinue_integration = False

        self.ODE_System  = ODE_System
        self.Units_class = Units_class

        self.Propagated_State     = []
        self.Independant_variable = []

        self.Surface_flag = False

    def run_integrator_step(self, State_vector, Independant_var):

        """
        
        TODO
        
        """

        Derivatives = zeros((self.RK45_size, len(State_vector)))

        for iteration in range(self.RK45_size):

            # Propagate the Intermediate State Vector
            Inter_State_Vector = State_vector + self.step * dot(self.Coeff_deriv[iteration], 
                                                                   Derivatives[:self.RK45_size - 1])
            
            # Advance the Intermediate value of the Independant Variable
            Inter_Independant_var = Independant_var + self.Coeff_ind_var[iteration] * self.step

            # Evaluate the ODE System
            Derivatives[iteration] = self.ODE_System.ODE_System(State_vector    = Inter_State_Vector, 
                                                                 Independant_var = Inter_Independant_var)
                                                                 
        # Evaluate the candidate O(h^5) solution
        sol_O5 = State_vector + self.step * dot(self.Coeff_sol, Derivatives)
        # Evaluate the test O(h^4) solution for error estimation
        sol_O4 = State_vector + self.step * dot(self.Coeff_test_sol, Derivatives)

        Est_State_Error = max(absolute(sol_O4 - sol_O5))

        if Est_State_Error < self.accuracy:

            current_step = self.step

            self.step = self.Safety_1 * self.step * power((self.accuracy / (Est_State_Error + self.Safety_2)), 0.2)
            self.contuinue_integration = True

            return sol_O5, Independant_var + current_step

        else:

            self.step = self.Safety_1 * self.step * power((self.accuracy / (Est_State_Error + self.Safety_2)), 0.25)
            self.contuinue_integration = False

       
            return State_vector, Independant_var


    def run_integrator(self, Integrator_radial_stop_point: float, Initial_State: list, Init_Independant_Var: float):

        State_Vector         = Initial_State
        Independant_Variable = Init_Independant_Var

        while True:

            State_Vector, Independant_Variable = self.run_integrator_step(State_vector    = State_Vector, 
                                                                          Independant_var = Independant_Variable)

            if self.contuinue_integration:

                self.Propagated_State.append(State_Vector)
                self.Independant_variable.append(Independant_Variable)

                if State_Vector[0] < self.Pressure_at_surface and not self.Surface_flag:
                    R_star = Independant_Variable * self.Units_class.Sun_Sch_radius_km
                    M_star = State_Vector[1]
                    self.Surface_flag = True

                if self.eval_stop_condition(Integrator_radial_stop_point, Independant_Variable, State_Vector):

                    return self.Propagated_State, self.Independant_variable, R_star, M_star
    
    def eval_stop_condition(self, Integrator_radial_stop_point, Independant_Variable, State_Vector):

        if Integrator_radial_stop_point != None:

            return Independant_Variable * self.Units_class.Sun_Sch_radius_km > Integrator_radial_stop_point

        else:

            return State_Vector[0] < self.Pressure_at_surface and self.Surface_flag