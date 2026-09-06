import numpy as np
import json
from math import sqrt, exp
import scipy.special as sp

a0 = 1.0        # working in ratios of the Bohr's radius which is equal to (h-bar)^2/(m*k_e*e^2)

# Radial part R_{n,l}(r)
def R(n, l, r):
    rho = 2.0 * r / (n * a0)
    norm = sqrt((2.0 / (n * a0))**3 * sp.factorial(n - l - 1) / (2 * n * sp.factorial(n + l)))
    L = sp.assoc_laguerre(rho, n - l - 1, 2 * l + 1)
    return norm * np.exp(-rho / 2.0) * rho**l * L

# Full complex orbital psi(r,theta,phi)
def psi(n, l, m, r, theta, phi):
    Y = sp.sph_harm_y(l, m, theta, phi) # complex spherical harmonic
    wavefunction = R(n, l, r) * Y
    return wavefunction

# Probability density = |psi|**2 for rejection sampling
def probability(n, l, m, r, theta, phi):
    return abs(psi(n, l, m, r, theta, phi))**2

# Sample points in space
def sample_points(n, l, m, N=50000):
    points = []
    maximum_probability = 0.0

    for _ in range(N):
        r = np.random.exponential(scale=n**2)
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)
        p = probability(n, l, m, r, theta, phi)

        if (p > maximum_probability):
            maximum_probability = p

        if np.random.rand() < p / maximum_probability:
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)

            wavefunction = psi(n, l, m, r, theta, phi)
            points.append({
                "x coordinate": x,
                "y coordinate": y,
                "z coordinate": z,
                "wavefunction real component": wavefunction.real,
                "wavefunction imaginary component": wavefunction.imag
            })
    return points

def main():
    print("=== Hydrogen Orbital Generator ===\n")
    n = int(input("Enter n (1..7): "))
    l = int(input(f"Enter l (0..{n-1}): "))
    m = int(input(f"Enter m (-{l}..{l}): "))
    N = int(input("Enter sample size, N: "))

    print("\nGenerating samples... please wait...")
    points = sample_points(n, l, m, N)

    filename = f"orbital_n{n}_l{l}_m{m}.json"
    with open(filename, "w") as f:
        json.dump(points, f, indent=2)

    print(f"\nSaved {len(points)} points to {filename}")

if __name__ == "__main__":
    main()
