import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
outfile = "s15model.hyd"

ds = yt.load(fname)
ad = ds.all_data()
Msun = 1.989e33

MBH = float(ds.parameters["sim_eff_mass"])/Msun  #mass cut in Msun
time = float(ds.parameters["time"])  #time in s

# compute cell faces and volume
r = np.array(ad[("flash","r")])
dr = np.array(ad[("flash", "dr")])
rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)

# fields
rho = np.array(ad[("flash","dens")])
temp = np.array(ad[("flash","temp")])
velx = np.array(ad[("flash","velx")])

dmass = cvol * rho / Msun       #cell mass in Msun
Menc = MBH + np.cumsum(dmass)

Nzon   = len(r)
Rcen   = rl[0]
rhoCen = rho[0]

with open(outfile, "w") as f:
    # first line:
    # timeStart  Nzon  mass_cut  Rcen  rhoCen
    f.write(f"{time:.8e} {Nzon:d} {MBH:.8e} {Rcen:.8e} {rhoCen:.8e}\n")

    # zone lines:
    # km  dMr  r(km)  rho(km)  Tp(km)  u(km)  Mr(km)  dummy
    for i in range(Nzon):
        f.write(
            f"{i+1:d} "
            f"{dmass[i]:.8e} "
            f"{rh[i]:.8e} "
            f"{rho[i]:.8e} "
            f"{temp[i]:.8e} "
            f"{velx[i]:.8e} "
            f"{Menc[i]:.8e} "
            f"{0.0:.8e}\n"
        )

print("Wrote", outfile)
print("Nzon =", Nzon)
print("time =", time)
print("mass_cut =", MBH)
print("ejecta mass =", np.sum(dmass))
print("total mass =", Menc[-1])
print("Rcen =", Rcen)
print("Rout =", rh[-1])
