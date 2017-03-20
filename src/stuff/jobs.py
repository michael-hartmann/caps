from numpy import linspace
from math import ceil
from convert import scaled

R          = 151.3e-6 # radius is 150µm
L_start    = 750e-9   # 750nm
L_end      = 150e-9   # 150nm
npts       = 20       # number of points
T          = 300      # in K
ωp         = 8.9      # in eV
γ          = 0.0357   # in eV

ldim_min  = 200
eta       = 7
cutoff    = 1e-9
epsrel    = 1e-6
cores     = 28

epsilon_F = 5e-6 # expected accuracy; important for numerical differentiation

cmd = "mpirun -n %d ./casimir -x %%.15g -L %%d -T %%.15g --omegap %%.15g --gamma %%.15g --cutoff %g --epsrel %g | tee %%04d_%%d.out" % (cores, cutoff, epsrel)
pfa = "./casimir_pfa %.15g %.15g %.15g %.15g"

for i,L in enumerate(linspace(L_start, L_end, npts)):
    LbyR = L/R

    # optimum value for differentiation
    hopt = LbyR*(epsilon_F/720)**(1/5)

    for j,h in enumerate((-2*hopt, -hopt, 0, hopt, 2*hopt)):
        ldim = max(int(ceil(eta/LbyR)),ldim_min)
        T_scaled, ωp_scaled, γ_scaled = scaled(R, LbyR+h, T, ωp, γ)

        print(cmd % (LbyR+h, ldim, T_scaled, ωp_scaled, γ_scaled, i,j))

    print()
