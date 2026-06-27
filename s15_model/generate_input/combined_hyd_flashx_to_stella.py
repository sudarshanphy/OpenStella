import yt
import numpy as np

num1 = 1718
fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_%04d"%(num1)

ncombine = 2   # number of zones to combine; use 1 for no combining

# These are indices AFTER the istart cut.
# Inclusive ranges: (200, 210) keeps 200,201,...,210 as original FLASH cells.
preserve_ranges = [(0,32)]

outfile = "s15model_new2_%04d_comb%02d.hyd"%(num1,ncombine)


def make_blocks(Nold, ncombine, preserve_ranges):
    preserve = np.zeros(Nold, dtype=bool)

    for a, b in preserve_ranges:
        a = max(a, 0)
        b = min(b, Nold - 1)
        if b >= a:
            preserve[a:b+1] = True

    blocks = []
    i = 0

    while i < Nold:
        if preserve[i]:
            blocks.append((i, i+1))
            i += 1
        else:
            j = i
            while j < Nold and (not preserve[j]) and (j - i) < ncombine:
                j += 1

            blocks.append((i, j))
            i = j

    return blocks


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
gpot = np.array(ad[("flash","gpot")])
eint = np.array(ad[("flash","eint")])

tener = 0.5 * rho * velx * velx + rho * eint + rho * gpot

rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)

dmass_all = cvol * rho / Msun
mass_g_all = cvol * rho
Menc_all = MBH0 + np.cumsum(dmass_all)

idx = np.argwhere(tener > 0.0).flatten()
if len(idx) == 0:
    raise RuntimeError("No zone found with tener > 0.0")

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
cvol = cvol[istart:]
dmass = dmass_all[istart:]
mass_g = mass_g_all[istart:]

Nold = len(r)

rl_new = []
rh_new = []
dmass_new = []
rho_new = []
temp_new = []
velx_new = []

blocks = make_blocks(Nold, ncombine, preserve_ranges)

for i, j in blocks:
    vol_sum = np.sum(cvol[i:j])
    mass_sum_g = np.sum(mass_g[i:j])
    mass_sum_msun = np.sum(dmass[i:j])

    rl_new.append(rl[i])
    rh_new.append(rh[j-1])
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
print("ncombine =", ncombine)
print("preserve_ranges =", preserve_ranges)
print("istart =", istart)
print("old kept Nzon =", Nold)
print("new Nzon =", Nzon)
print("time =", time)
print("mass_cut =", MBH)
print("ejecta mass =", np.sum(dmass))
print("total mass =", Menc[-1])
print("Rcen =", Rcen)
print("Rout =", rh[-1])
print("rhoCen =", rhoCen)
print("v first =", velx[0])
print("v min/max kept =", velx.min(), velx.max())

print("first 20 blocks:")
for b in blocks[:20]:
    print(b, "size", b[1]-b[0])
