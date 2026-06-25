import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
species_file = "SpeciesList.txt"
outfile = "s15model.abn"

ds = yt.load(fname)
ad = ds.all_data()

rho = np.array(ad[("flash","dens")])
velx = np.array(ad[("flash","velx")])
gpot = np.array(ad[("flash","gpot")])
eint = np.array(ad[("flash","eint")])

tener = 0.5 * rho * velx * velx + rho * eint + rho * gpot
idx = np.argwhere(tener > 0).flatten()

if len(idx) == 0:
    raise RuntimeError("No zone found with tener > 0.0")

istart = idx[0]

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
    print(s)
    X[s] = np.array(ad[("flash", s.ljust(4))])

Nold = len(next(iter(X.values())))

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

H      = H[istart:]
He     = He[istart:]
C      = C[istart:]
N      = N[istart:]
O      = O[istart:]
Ne     = Ne[istart:]
Na     = Na[istart:]
Mg     = Mg[istart:]
Al     = Al[istart:]
Si     = Si[istart:]
S      = S[istart:]
Ar     = Ar[istart:]
Ca     = Ca[istart:]
FePeak = FePeak[istart:]
Ni58   = Ni58[istart:]
Ni56   = Ni56[istart:]

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
print("istart =", istart)
print("old Nzon =", Nold)
print("new Nzon =", Nzon)
print("Xsum min/max =", Xsum.min(), Xsum.max())
print("Ni56 sum =", Ni56.sum())
