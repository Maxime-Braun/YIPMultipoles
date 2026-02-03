import math

def get_real_ylm(l, m, theta, phi):
    # Standard Real Spherical Harmonic
    # Y00 = 1/sqrt(4pi)
    if l == 0 and m == 0:
        return 0.5 * math.sqrt(1.0 / math.pi)
    return 0

val = get_real_ylm(0, 0, 0, 0)
weight = 1.0 * val * 4 * math.pi
print(f"Y00 value: {val}")
print(f"Integral of 1 * Y00: {weight}")
