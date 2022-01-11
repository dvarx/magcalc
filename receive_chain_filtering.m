%band stop dimensioning
dcoil=3e-2;
N=100;
h=3e-2;
mu0=4*pi*1e-7;

A=(dcoil/2)^2*pi;
Lbs=N^2*mu0*A/h;

f0=10e3;
Cbs=1/((2*pi*f0)^2*Lbs);