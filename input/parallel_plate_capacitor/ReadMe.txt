###################################
#### parallel_plate_capacitor #####
###################################

  -This example models a parallel plate capacitor such that
   a voltage of 10 V is specified at Zmax boundary and 0 V
   is specified at Zmin boundary.

  -First four examples provided in this folder show different ways
   to specify these voltages.

     1) inhomo_const_dirichlet:
            This example specifies voltages directly through where
        boundary conditions are specified, e.g. dir(10) or dir(0).

     2) function_parsed_dirichlet:
            This example specifies voltages through dirichlet function
        parsers whose name is specified at the time of specifying boundaries,
        e.g. dir(Zmin) or dir(Zmax). Function parsers are specified for
        Zmin and Zmax.

     3) function_parsed_robin:
            This example specifies voltages through robin boundaries by
        creating parsers for a, b, and f coefficients of the boundary.

     4) function_parsed_dirichlet_varyingFunc:
            This example is similar to example 2 but voltage at the top plate
        is applied as a cosine function 10*cos(2 pi x/(2Lx)).

  -Fifth example, is for the case with two dielectrics, where in upper Lz/2 region there is air
   and in the lower, SiO2, epsilon_2 = 3.8.

  -On all lateral sides we can either specify Neumann boundaries
   or periodic boundaries.

For verification:
-Check that voltage in the ghost cells is as specified.
-Note:
 Charge on the plate, Q = C V_0 (where V_0 = 10 V in our case)
 Surface charge density, sigma = Q/A (A = 0.2**2 in our case)
 From Gauss' law, sigma = D (D_z in our 1-D example)
 For single dielectric, D_z = E_z / epsilon_0

 From above equations, E_z = D_z/epsilon_0  = sigma/epsilon_0 = Q/(A*epsilon_0) = (C V_0)/(A*epsilon_0)
 So, C = E_z *A*epsilon_0 / V_0

 Now from solution, we can find that |E_z| = 200 V/m
 By calculating, we obtain C = 7.04 pF (as expected from theory).

-For the two dielectric case:
 potential phi at the interface is given by,
 phi_int = epsilon_1 d_2 / (epsilon_2 d_1 + epsilon_1 d_2) V_0

 In our example, epsilon_1 = epsilon_0 and epsilon_2 = 3.8*epsilon_0, d_1 = d_2 = L_z/2
 this gives, phi_int = 2.04 V

 By probing cell-centered phi at the air-dielectric interface we find this phi_int = 2.07 V.

-E_z in the air and dielectric can also be easily estimated by taking the gradient, e.g.
 E_z1 = - (10 - 2.04) / d_1 = -317.2 (Theory: -317.34)
 E_z2 = - (2.04 - 0) / d_2 = -82.8 (Theory: -83.51)
 These values are close to what we see from the solution.



