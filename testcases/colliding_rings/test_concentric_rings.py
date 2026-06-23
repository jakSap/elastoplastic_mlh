#!/usr/bin/env python3
"""
Quick test: build a ring by stacking concentric particle circles.
Each circle at radius r gets N = round(2*pi*r / delta_p) particles,
evenly spaced angularly.  The radial step between circles equals delta_p.
"""

import numpy as np
import matplotlib.pyplot as plt

delta_p = 0.2
r_inner = 3.0
r_outer = 4.0

xs, ys = [], []

r = r_inner
while r <= r_outer + 1e-9:
    n = max(1, round(2 * np.pi * r / delta_p))
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs.extend(r * np.cos(angles))
    ys.extend(r * np.sin(angles))
    r += delta_p

xs = np.array(xs)
ys = np.array(ys)
print(f"Total particles: {len(xs)}")

fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
ax.scatter(xs, ys, s=1.5, linewidths=0)
ax.set_aspect("equal")
ax.set_title(f"Concentric-ring IC  (delta_p={delta_p}, r=[{r_inner},{r_outer}])\n"
             f"N = {len(xs)}")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.tight_layout()
plt.savefig("test_concentric_rings.png")
print("Saved test_concentric_rings.png")
plt.show()
