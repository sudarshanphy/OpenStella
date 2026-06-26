import sys
import numpy as np
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "res/s15model_comb16.lbol"

rows = []

with open(fname, "r") as f:
    for line in f:
        parts = line.split()
        if len(parts) < 2:
            continue

        try:
            vals = [float(x.replace("D", "E").replace("d", "e")) for x in parts]
            rows.append(vals)
        except ValueError:
            continue

data = np.array(rows)

t = data[:, 0]
L = data[:, 1]

# If STELLA wrote log10(L) instead of L, convert automatically.
# Real luminosity is usually ~1e35--1e45 erg/s, while log10(L) is ~35--45.
if np.nanmax(L) < 100.0:
    print("Second column looks like log10(L). Converting to erg/s.")
    L = 10.0**L

mask = np.isfinite(t) & np.isfinite(L) & (L > 0.0)
t = t[mask]
L = L[mask]

plt.figure(figsize=(7, 5))
plt.plot(t, L, lw=1.5)
plt.yscale("log")
plt.xlabel("Time [days]")
plt.ylabel(r"Bolometric luminosity [erg s$^{-1}$]")
plt.xlim(t.min(), t.max())
plt.tight_layout()

outpng = fname.split("/")[-1].replace(".lbol", "_lbol.png")
outpdf = fname.split("/")[-1].replace(".lbol", "_lbol.pdf")

plt.savefig(outpng, dpi=200)
plt.savefig(outpdf)

print("Read:", fname)
print("time range =", t.min(), t.max(), "days")
print("L range =", L.min(), L.max(), "erg/s")
print("Saved:", outpng)
print("Saved:", outpdf)
