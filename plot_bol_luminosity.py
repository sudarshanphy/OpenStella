import sys
import numpy as np
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "res/s15model_1718_comb08.lbol"

TIME_OFFSET_DAYS = 5.13880835648

rows = []

with open(fname, "r") as f:
    for line in f:
        parts = line.split()
        if len(parts) < 5:
            continue

        try:
            vals = [float(x.replace("D", "E").replace("d", "e")) for x in parts]
            rows.append(vals)
        except ValueError:
            continue

data = np.array(rows)

t = data[:, 0] + TIME_OFFSET_DAYS

L_ubvri = 10.0**data[:, 1]
L_bol   = 10.0**data[:, 2]
L_xeuv  = 10.0**data[:, 3]
L_ir    = 10.0**data[:, 4]

mask = np.isfinite(t) & np.isfinite(L_bol) & (L_bol > 0.0)

t = t[mask]
L_bol = L_bol[mask]
L_ubvri = L_ubvri[mask]
L_xeuv = L_xeuv[mask]
L_ir = L_ir[mask]

plt.figure(figsize=(7, 5))
plt.plot(t, L_bol, lw=1.8, label="Bolometric")
#plt.plot(t, L_ubvri, lw=1.2, label="UBVRI")
#plt.plot(t, L_xeuv, lw=1.2, label="XEUV")
#plt.plot(t, L_ir, lw=1.2, label="IR")

plt.yscale("log")
plt.xlabel("Physical time [days]")
plt.ylabel(r"Luminosity [erg s$^{-1}$]")
plt.xlim(t.min(), t.max())
plt.legend()
plt.tight_layout()

outpng = fname.split("/")[-1].replace(".lbol", "_lbol.png")
outpdf = fname.split("/")[-1].replace(".lbol", "_lbol.pdf")

plt.savefig(outpng, dpi=200)
plt.savefig(outpdf)

print("Read:", fname)
print("time range =", t.min(), t.max(), "days")
print("Lbol range =", L_bol.min(), L_bol.max(), "erg/s")
print("Saved:", outpng)
print("Saved:", outpdf)
