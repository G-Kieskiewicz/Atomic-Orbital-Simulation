import numpy as np
# import json -- This line is useless now as i dont save the data as a json
from math import sqrt
import scipy.special as sp
import pyvista as pv

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

def plot_orbital(filename,coordinates, real_parts, point_size=0.05, background="#000000", opacity=1,
                  positive_color="#ff8a4c", negative_color="#49c8ff"):
    cloud = pv.PolyData(coordinates)

    pos_rgb = pv.Color(positive_color).int_rgb  # hex -> (r, g, b) 0-255
    neg_rgb = pv.Color(negative_color).int_rgb

    colors = np.where(
        real_parts[:, None] >= 0,
        np.array(pos_rgb),
        np.array(neg_rgb)
    ).astype(np.uint8)

    cloud["colors"] = colors

    plotter = pv.Plotter()
    plotter.set_background(background)  # hex works directly here
    plotter.add_mesh(
        cloud,
        scalars="colors",
        rgb=True,
        point_size=point_size,
        opacity=opacity,
        render_points_as_spheres=False,
    )
    plotter.export_html(filename)
    plotter.show()

def main():
    print("=== Hydrogen Orbital Generator ===\n")
    n = int(input("Enter n (1..7): "))
    l = int(input(f"Enter l (0..{n-1}): "))
    m = int(input(f"Enter m (-{l}..{l}): "))
    N = int(input("Enter sample size, N: "))

    print("\nGenerating samples... please wait...")
    points = sample_points(n, l, m, N)

    filename = f"orbital_n{n}_l{l}_m{m}.html"
    """with open(filename, "w") as f:
        json.dump(points, f, indent=2)"""

    coordinates = np.array([[d["x coordinate"], d["y coordinate"], d["z coordinate"]] for d in points])
    real_parts = np.array([d["wavefunction real component"] for d in points])

    plot_orbital(filename, coordinates, real_parts)

    print(f"\nSaved simulation to {filename}")

if __name__ == "__main__":
    main()
