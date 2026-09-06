# Atomic Orbital Simulation

Generates the wavefunction probability densities for hydrogen like atoms and renders them as 3D point clouds via PyVista.

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [API](#api)
- [Contributing](#contributing)
- [License](#license)

## Background

This project builds an atomic orbital simulation from first principles through solving the radial and angular parts of the hydrogen Schrödinger equation, sampling the resulting |ψ|² probability density via rejection sampling, and visualizing the result in 3D.

`Atomic_Orbital_Sim.py` computes the wavefunction from scratch using `scipy.special` (associated Laguerre polynomials, spherical harmonics) and renders/exports the point cloud using `pyvista`.

Built for fun after a YouTube rabbit hole on atomic orbitals. I wanted to actually compute and visualize the wavefunctions myself instead of just taking the textbook pictures for granted.

## Install

### Dependencies

The Python pipeline needs:

```bash
pip install numpy scipy pyvista
```

## Usage

### CLI

Run the Python sampler directly:

```bash
python Atomic_Orbital_Sim.py
```

It will prompt for the quantum numbers and sample size:

```
=== Hydrogen Orbital Generator ===

Enter n (1..7): 2
Enter l (0..1): 1
Enter m (-1..1): 0
Enter sample size, N: 50000
```

This produces `orbital_n2_l1_m0.html` (an exported interactive PyVista scene) and opens a live render window, the live render window doesn't always work due to PyVista acting up.

## API

The Python module exposes these functions if imported rather than run as a script:

- `R(n, l, r)` — radial wavefunction `R_{n,l}(r)`
- `psi(n, l, m, r, theta, phi)` — full complex wavefunction ψ
- `probability(n, l, m, r, theta, phi)` — |ψ|² probability density
- `sample_points(n, l, m, N=50000)` — rejection-samples `N`-ish points from the probability density, returns a list of dicts with coordinates and ψ's real/imaginary components
- `plot_orbital(filename, coordinates, real_parts, ...)` — renders and exports a point cloud, colored by the sign of Re(ψ)

## Contributing

This is a personal learning project. Questions and bug reports are welcome via GitHub issues, but this project will not be maintained and will not accept pull requests or other contributions. It may see further commits from me as I add other atomic simulations.

## License

MIT © Gabriel A. Kieskiewicz