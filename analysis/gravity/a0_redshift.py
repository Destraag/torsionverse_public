"""
a0_redshift.py — Item A unscripted #1: a0 redshift evolution prediction.

Framework prediction: a0(z) = R_s * c * H(z)
Since R_s is a geometric constant, a0 tracks H(z) exactly.

At z=1: H(z=1) = H0 * sqrt(Omega_m*(1+z)^3 + Omega_Lambda)
Using Planck 2018: Omega_m=0.315, Omega_Lambda=0.685

This is a testable prediction: galaxy rotation curves at z~1 should show
a MOND transition acceleration 76% higher than at z=0.
JWST extended rotation curve programs can test this.

Run: python analysis/a0_redshift.py
"""

import math

SEP = "=" * 58

Rs      = math.sqrt(5) / (4 * math.pi)
c       = 2.998e8          # m/s
H0      = 2.184e-18        # 1/s  (67.4 km/s/Mpc)

# Planck 2018 cosmological parameters
Omega_m      = 0.315
Omega_Lambda = 0.685

def H(z):
    """Hubble parameter at redshift z (flat LambdaCDM)."""
    return H0 * math.sqrt(Omega_m * (1 + z)**3 + Omega_Lambda)

def a0_at_z(z):
    return Rs * c * H(z)

print(SEP)
print("a0 REDSHIFT EVOLUTION — FRAMEWORK PREDICTION")
print(SEP)
print()
print(f"R_s          = {Rs:.6f}")
print(f"H0           = {H0:.4e} 1/s  (67.4 km/s/Mpc)")
print(f"Omega_m      = {Omega_m}")
print(f"Omega_Lambda = {Omega_Lambda}")
print()

a0_z0 = a0_at_z(0)
print(f"a0(z=0) = R_s * c * H0 = {a0_z0:.4e} m/s^2")
print()

redshifts = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

print(f"  {'z':>5}  {'H(z)/H0':>10}  {'a0(z) [m/s^2]':>16}  {'a0(z)/a0(0)':>12}  Note")
print(f"  {'-'*5}  {'-'*10}  {'-'*16}  {'-'*12}  ----")
for z in redshifts:
    hz  = H(z)
    a0z = a0_at_z(z)
    ratio = a0z / a0_z0
    note = ""
    if z == 0.0: note = "z=0 baseline"
    if z == 1.0: note = "<<< WHITEPAPER CLAIM"
    print(f"  {z:>5.1f}  {hz/H0:>10.4f}  {a0z:>16.4e}  {ratio:>12.4f}  {note}")

print()
Hz1  = H(1)
ratio_z1 = Hz1 / H0
print(f"H(z=1) / H0 = sqrt({Omega_m}*8 + {Omega_Lambda}) = sqrt({Omega_m*8 + Omega_Lambda:.3f}) = {ratio_z1:.4f}")
print(f"a0(z=1) = {ratio_z1:.4f} * a0(z=0) = {a0_at_z(1):.4e} m/s^2")
print()
print(f"Whitepaper claims: H(z=1)/H0 = 1.76 (script gives {ratio_z1:.4f})")
print(f"a0(z=1) = 1.76 * a0(z=0) -> {1.761 * a0_z0:.4e} m/s^2 (whitepaper)")
print(f"a0(z=1) = {ratio_z1:.4f} * a0(z=0) -> {a0_at_z(1):.4e} m/s^2 (script)")

print()
print(SEP)
print("TESTABILITY:")
print()
print("  This prediction is testable with JWST extended rotation curves.")
print("  If galaxies at z~1 show flat rotation curves transitioning at")
print(f"  a0~{a0_at_z(1):.2e} m/s^2 rather than {a0_z0:.2e} m/s^2,")
print("  that would confirm the medium is cosmologically coupled.")
print("  A constant a0 across redshift would require framework modification.")
print(SEP)
