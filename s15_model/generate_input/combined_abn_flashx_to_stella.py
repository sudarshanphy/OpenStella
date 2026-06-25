import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
species_file = "SpeciesList.txt"
outfile = "s15model_comb16.abn"

ncombine = 16   # must match the hyd file

ds = yt.load(fname)
ad = ds.all_data()
Msun = 1.989e33

r = np.array(ad[("flash","r")])
dr = np.array(ad[("flash","dr")])
rho = np.array(ad[("flash","dens")])

rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)
mass_g = cvol * rho

species = []
with open(species_file, "r") as f:
    for line in f:
        s = line.strip()
        if s:
            species.append(s)

def elem(s):
    if s == "p" or s == "d":
        return "h"
    if s == "n":
        return "free_n"

    out = ""
    for c in s:
        if c.isalpha():
            out += c
        else:
            break
    return out

X = {}
for s in species:
    X[s] = np.array(ad[("flash", s.ljust(4))])

Nold = len(r)

H      = np.zeros(Nold)
He     = np.zeros(Nold)
C      = np.zeros(Nold)
N      = np.zeros(Nold)
O      = np.zeros(Nold)
Ne     = np.zeros(Nold)
Na     = np.zeros(Nold)
Mg     = np.zeros(Nold)
Al     = np.zeros(Nold)
Si     = np.zeros(Nold)
S      = np.zeros(Nold)
Ar     = np.zeros(Nold)
Ca     = np.zeros(Nold)
FePeak = np.zeros(Nold)
Ni58   = np.zeros(Nold)
Ni56   = np.zeros(Nold)

used = set()

for s in species:
    e = elem(s)

    if e == "h":
        H += X[s]
        used.add(s)
    elif e == "he":
        He += X[s]
        used.add(s)
    elif e == "c":
        C += X[s]
        used.add(s)
    elif e == "n":
        N += X[s]
        used.add(s)
    elif e == "o":
        O += X[s]
        used.add(s)
    elif e == "ne":
        Ne += X[s]
        used.add(s)
    elif e == "na":
        Na += X[s]
        used.add(s)
    elif e == "mg":
        Mg += X[s]
        used.add(s)
    elif e == "al":
        Al += X[s]
        used.add(s)
    elif e == "si":
        Si += X[s]
        used.add(s)
    elif e == "s":
        S += X[s]
        used.add(s)
    elif e == "ar":
        Ar += X[s]
        used.add(s)
    elif e == "ca":
        Ca += X[s]
        used.add(s)
    elif s == "ni56":
        Ni56 += X[s]
        used.add(s)
    elif s == "ni58":
        Ni58 += X[s]
        used.add(s)

for s in species:
    if s not in used:
        FePeak += X[s]

groups = [H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56]

groups_new = [[] for _ in groups]

for i in range(0, Nold, ncombine):
    j = min(i + ncombine, Nold)
    msum = np.sum(mass_g[i:j])

    for k, G in enumerate(groups):
        groups_new[k].append(np.sum(G[i:j] * mass_g[i:j]) / msum)

groups_new = [np.array(g) for g in groups_new]

H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56 = groups_new

Nzon = len(H)

with open(outfile, "w") as f:
    for i in range(Nzon):
        f.write(
            f"{i+1:d} "
            f"{0.0:.8e} {0.0:.8e} {0.0:.8e} "
            f"{H[i]:.8e} "
            f"{He[i]:.8e} "
            f"{C[i]:.8e} "
            f"{N[i]:.8e} "
            f"{O[i]:.8e} "
            f"{Ne[i]:.8e} "
            f"{Na[i]:.8e} "
            f"{Mg[i]:.8e} "
            f"{Al[i]:.8e} "
            f"{Si[i]:.8e} "
            f"{S[i]:.8e} "
            f"{Ar[i]:.8e} "
            f"{Ca[i]:.8e} "
            f"{FePeak[i]:.8e} "
            f"{Ni58[i]:.8e} "
            f"{Ni56[i]:.8e}\n"
        )

Xsum = H + He + C + N + O + Ne + Na + Mg + Al + Si + S + Ar + Ca + FePeak + Ni58 + Ni56

print("Wrote", outfile)
print("ncombine =", ncombine)
print("old Nzon =", Nold)
print("new Nzon =", Nzon)
print("Xsum min/max =", Xsum.min(), Xsum.max())
print("Ni56 sum =", Ni56.sum())
