import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
outfile = "s15model_comb16.hyd"

ncombine = 16   # number of FLASH zones combined into one STELLA zone

ds = yt.load(fname)
ad = ds.all_data()
Msun = 1.989e33

MBH = float(ds.parameters["sim_eff_mass"])/Msun
time = float(ds.parameters["time"])

r = np.array(ad[("flash","r")])
dr = np.array(ad[("flash","dr")])

rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)

rho = np.array(ad[("flash","dens")])
temp = np.array(ad[("flash","temp")])
velx = np.array(ad[("flash","velx")])

mass_g = cvol * rho
dmass = mass_g / Msun

rl_new = []
rh_new = []
vol_new = []
dmass_new = []
rho_new = []
temp_new = []
velx_new = []

Nold = len(r)

for i in range(0, Nold, ncombine):
    j = min(i + ncombine, Nold)

    vol_sum = np.sum(cvol[i:j])
    mass_sum_g = np.sum(mass_g[i:j])
    mass_sum_msun = mass_sum_g / Msun

    rl_new.append(rl[i])
    rh_new.append(rh[j-1])
    vol_new.append(vol_sum)
    dmass_new.append(mass_sum_msun)

    rho_new.append(mass_sum_g / vol_sum)
    temp_new.append(np.sum(temp[i:j] * mass_g[i:j]) / mass_sum_g)
    velx_new.append(np.sum(velx[i:j] * mass_g[i:j]) / mass_sum_g)

rl = np.array(rl_new)
rh = np.array(rh_new)
dmass = np.array(dmass_new)
rho = np.array(rho_new)
temp = np.array(temp_new)
velx = np.array(velx_new)

Menc = MBH + np.cumsum(dmass)

Nzon = len(dmass)
Rcen = rl[0]
rhoCen = rho[0]

with open(outfile, "w") as f:
    f.write(f"{time:.8e} {Nzon:d} {MBH:.8e} {Rcen:.8e} {rhoCen:.8e}\n")

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
print("ncombine =", ncombine)
print("old Nzon =", Nold)
print("new Nzon =", Nzon)
print("time =", time)
print("mass_cut =", MBH)
print("ejecta mass =", np.sum(dmass))
print("total mass =", Menc[-1])
print("Rcen =", Rcen)
print("Rout =", rh[-1])
