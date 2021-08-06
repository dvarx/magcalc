# Magcalc TNBMNS

This repository contains code to interpolate simulation data received from Comsol Multiphysics as well as processing this data.

## Input Data

Input data are `.csv` files (`bxs.csv`,`bys.csv`,`bzs.csv`) of the following form:

```
% Model,mpi_sim_single_coils.mph
% Version,COMSOL 5.4.0.246
% Date,"Jul 30 2021, 09:29"
% Table,Table 2 - Point Evaluation 1 (mf.Bx)
% il1,il2,il3,ir1,ir2,ir3,"Magnetic flux density, x component (T), Point: (-0.03, -0.03, -0.03)","Magnetic flux density, x component (T), ...........
0,0,0,1,0,0,-1.694218193203235E-4,-1.3731514733048166E-4,
0,0,0,0,0,1,-4.015639233810901E-5,8.586807374602612E-5,
0,0,0,0,1,0,0.0013807907493288253,0.001286790760573608,
0,0,1,0,0,0,-3.411286968014771E-5,-3.022815546844687E-5,
0,1,0,0,0,0,-2.056718215456974E-4,-2.2116032333843955E-4,
1,0,0,0,0,0,-3.6135806284982603E-4,-4.06375129403028E-4,
```

