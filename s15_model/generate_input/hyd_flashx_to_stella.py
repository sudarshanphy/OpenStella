import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
outfile = "s15model.hyd"

ds = yt.load(fname)
ad = ds.all_data()
Msun = 1.989e33

MBH0 = float(ds.parameters["sim_eff_mass"])/Msun
time = float(ds.parameters["time"])

r = np.array(ad[("flash","r")])
dr = np.array(ad[("flash", "dr")])

rho = np.array(ad[("flash","dens")])
temp = np.array(ad[("flash","temp")])
velx = np.array(ad[("flash","velx")])

rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)

dmass_all = cvol * rho / Msun
Menc_all = MBH0 + np.cumsum(dmass_all)

idx = np.argwhere(velx > 1.0e3).flatten()
if len(idx) == 0:
    raise RuntimeError("No zone found with velx > 1.0e3")

istart = idx[0]

if istart == 0:
    MBH = MBH0
else:
    MBH = Menc_all[istart-1]

r = r[istart:]
dr = dr[istart:]
rl = rl[istart:]
rh = rh[istart:]
rho = rho[istart:]
temp = temp[istart:]
velx = velx[istart:]
dmass = dmass_all[istart:]

Menc = MBH + np.cumsum(dmass)

Nzon = len(r)
Rcen = rl[0]
rhoCen = rho[0]

with open(outfile, "w") as f:
    f.write(f"{time:.8e} {Nzon:d} {MBH:.8e} {Rcen:.8e} {rhoCen:.8e}\n")

    for k in range(Nzon):
        f.write(
            f"{k+1:d} "
            f"{dmass[k]:.8e} "
            f"{rh[k]:.8e} "
            f"{rho[k]:.8e} "
            f"{temp[k]:.8e} "
            f"{velx[k]:.8e} "
            f"{Menc[k]:.8e} "
            f"{0.0:.8e}\n"
        )

print("Wrote", outfile)
print("istart =", istart)
print("Nzon =", Nzon)
print("time =", time)
print("mass_cut =", MBH)
print("ejecta mass =", np.sum(dmass))
print("total mass =", Menc[-1])
print("Rcen =", Rcen)
print("Rout =", rh[-1])
print("rhoCen =", rhoCen)
print("v first =", velx[0])
print("v min/max kept =", velx.min(), velx.max())
