import yt
import numpy as np

fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_0958"
species_file = "SpeciesList.txt"
outfile = "s15model.abn"

ds = yt.load(fname)
ad = ds.all_data()

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

Nzon = len(next(iter(X.values())))

H      = np.zeros(Nzon)
He     = np.zeros(Nzon)
C      = np.zeros(Nzon)
N      = np.zeros(Nzon)
O      = np.zeros(Nzon)
Ne     = np.zeros(Nzon)
Na     = np.zeros(Nzon)
Mg     = np.zeros(Nzon)
Al     = np.zeros(Nzon)
Si     = np.zeros(Nzon)
S      = np.zeros(Nzon)
Ar     = np.zeros(Nzon)
Ca     = np.zeros(Nzon)
FePeak = np.zeros(Nzon)
Ni58   = np.zeros(Nzon)
Ni56   = np.zeros(Nzon)

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
print("Nzon =", Nzon)
print("Xsum min/max =", Xsum.min(), Xsum.max())
print("Ni56 sum =", Ni56.sum())
