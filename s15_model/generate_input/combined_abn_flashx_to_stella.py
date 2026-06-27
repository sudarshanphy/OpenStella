import yt
import numpy as np

num1 = 1718
fname = "../chk/failed_1dsph_ug_s15_18274_csm_wenoexp1_hdf5_chk_%04d"%(num1)
species_file = "SpeciesList.txt"
ncombine = 2   # must match hyd_flashx_to_stella.py

# These are indices AFTER the istart cut.
# Inclusive ranges
# Must match the .hyd script.
preserve_ranges = [(0,56)]

outfile = "s15model_new3_%04d_comb%02d.abn"%(num1, ncombine)


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

r = np.array(ad[("flash","r")])
dr = np.array(ad[("flash","dr")])
rho = np.array(ad[("flash","dens")])
velx = np.array(ad[("flash","velx")])
gpot = np.array(ad[("flash","gpot")])
eint = np.array(ad[("flash","eint")])

rl = r - 0.5 * dr
rh = r + 0.5 * dr
cvol = 4.0/3.0 * np.pi * dr * (rl*rl + rl*rh + rh*rh)
mass_g = cvol * rho

tener = 0.5 * rho * velx * velx + rho * eint + rho * gpot
#idx = np.argwhere(tener > 0.0).flatten()
idx = np.argwhere(velx > 1.0e6).flatten()

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

Nold_all = len(next(iter(X.values())))

H      = np.zeros(Nold_all)
He     = np.zeros(Nold_all)
C      = np.zeros(Nold_all)
N      = np.zeros(Nold_all)
O      = np.zeros(Nold_all)
Ne     = np.zeros(Nold_all)
Na     = np.zeros(Nold_all)
Mg     = np.zeros(Nold_all)
Al     = np.zeros(Nold_all)
Si     = np.zeros(Nold_all)
S      = np.zeros(Nold_all)
Ar     = np.zeros(Nold_all)
Ca     = np.zeros(Nold_all)
FePeak = np.zeros(Nold_all)
Ni58   = np.zeros(Nold_all)
Ni56   = np.zeros(Nold_all)

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

mass_g = mass_g[istart:]

groups = [H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56]
groups_new = [[] for _ in groups]

Nold = len(H)

blocks = make_blocks(Nold, ncombine, preserve_ranges)

for i, j in blocks:
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
print("preserve_ranges =", preserve_ranges)
print("istart =", istart)
print("old full Nzon =", Nold_all)
print("old kept Nzon =", Nold)
print("new Nzon =", Nzon)
print("Xsum min/max =", Xsum.min(), Xsum.max())
print("Ni56 sum =", Ni56.sum())

print("first 20 blocks:")
for b in blocks[:20]:
    print(b, "size", b[1]-b[0])
