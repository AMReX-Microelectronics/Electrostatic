###################################
#### parallel_plate_capacitor #####
###################################
 
  -This example models a parallel plate capacitor such that 
   a voltage of 10 V is specified at Zmax boundary and 0 V 
   is specified at Zmin boundary.
 
  -Three examples provided in this folder show different ways
   to specify these voltages.

     1) inhomo_const_dirichlet: 
            This example specifies voltages directly through where
        boundary conditions are specified, e.g. dir(10) or dir(0).

     2) function_parsed_dirichlet:
            This example specifies voltages through dirichlet function 
        parsers whose name is specied at the time of specifying boundaries,
        e.g. dir(Zmin) or dir(Zmax). Function parsers are specified for 
        Zmin and Zmax.

     3) function_parsed_robin:
            This example specifies voltages through robin boundaries by
        creating parsers for a, b, and f coefficients of the boundary.

  -On all lateral sides we can either specify Neumann boundaries 
   or periodic boundaries.
