import numpy as np
import pyvista as pv
import json

def load_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)

    points = np.array([[d["x coordinate"], d["y coordinate"], d["z coordinate"]] for d in data])
    real_parts = np.array([d["wavefunction real component"] for d in data])

    return points, real_parts

def plot_orbital(points, real_parts, point_size=0.05, background="#000000", opacity=1,
                  positive_color="#ff8a4c", negative_color="#49c8ff"):
    cloud = pv.PolyData(points)

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
    plotter.export_html("orbital_n4_l3_m3.html")
    plotter.show()

points, real_parts = load_data("orbital_n4_l3_m3.json")
plot_orbital(points, real_parts)