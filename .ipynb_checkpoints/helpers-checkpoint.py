# beam_utils.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc, norm
from matplotlib.gridspec import GridSpec

__all__ = [
    "compute_hist",
    "compute_cutout",
    "plot_hist",
    "hist",            # high-level convenience
]

# --------------------------- Helpers ---------------------------

def sobol_uniform(n, d):
    sampler = qmc.Sobol(d, scramble=False).fast_forward(1)
    return sampler.random(n)
    
def sobol_normals(n, d, eps=1e-12):
    U = sobol_uniform(n, d)
    # Clamp to avoid 0 or 1, which would lead to ±inf in inverse CDF
    U = np.clip(U, eps, 1 - eps)
    Z = norm.ppf(U)  # Elementwise inverse CDF (aka "probit")
    return Z

def compute_hist(x, y, nx=80, ny=80, xlims=(-0.2, 0.2), ylims=(-0.2, 0.2)):
    """
    Compute a NaN-masked 2D histogram and 1D marginals.
    Returns a dict with hist, edges, centers, and marginals.
    """
    if xlims is None:
        xlims = (np.min(x), np.max(x))
    if ylims is None:
        ylims = (np.min(y), np.max(y))

    x_edges = np.linspace(xlims[0], xlims[1], nx + 1)
    y_edges = np.linspace(ylims[0], ylims[1], ny + 1)

    H, x_edges, y_edges = np.histogram2d(x, y, bins=[x_edges, y_edges])
    H = H.astype(float)
    H[H == 0] = np.nan  # mask empty bins like Julia's NaN convention

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    x_marginal = np.nansum(H, axis=1)  # over y
    y_marginal = np.nansum(H, axis=0)  # over x

    return dict(
        hist=H,
        x_edges=x_edges, y_edges=y_edges,
        x_centers=x_centers, y_centers=y_centers,
        x_marginal=x_marginal, y_marginal=y_marginal,
    )


def compute_cutout(H, x_edges, y_edges, h=0.05, v=0.05, center=(0.0, 0.0)):
    """
    Extract a zero-NaN-preserving rectangular cutout of H around (x0,y0)
    with half-widths (h, v). Returns the cutout array plus metadata.
    """
    x0, y0 = center

    # step sizes
    h_step = x_edges[1] - x_edges[0]
    v_step = y_edges[1] - y_edges[0]

    # index of center in bin space
    h_center = int(round((x0 - x_edges[0]) / h_step))
    v_center = int(round((y0 - y_edges[0]) / v_step))

    half_h = int(round(h / h_step))
    half_v = int(round(v / v_step))

    # clamp to array bounds (H is shape (nx, ny))
    h_start = max(h_center - half_h, 0)
    h_stop  = min(h_center + half_h, H.shape[0])
    v_start = max(v_center - half_v, 0)
    v_stop  = min(v_center + half_v, H.shape[1])

    # transpose so rows=y, cols=x for image-like use
    cut = H.T[v_start:v_stop, h_start:h_stop]

    return dict(
        cutout=cut,
        center=center,
    )


# --------------------------- Plotting ---------------------------

def plot_hist(data, title="Beam profile at Target",
              box=None, cutout=None, cmap="inferno"):
    """
    Plot 2D histogram with marginals.
    Optionally overlay a rectangular box and/or add a cutout inset.
      - box: dict with keys {'h_bounds':(xmin,xmax), 'v_bounds':(ymin,ymax)}
      - cutout: dict returned by compute_cutout (must contain 'cutout', 'h_bounds', 'v_bounds')
    Returns the Matplotlib figure.
    """
    H = data["hist"]
    x_edges = data["x_edges"]; y_edges = data["y_edges"]
    x_centers = data["x_centers"]; y_centers = data["y_centers"]
    x_marginal = data["x_marginal"]; y_marginal = data["y_marginal"]

    fig = plt.figure(figsize=(8.5, 8.5), constrained_layout=True)
    gs = GridSpec(2, 2, width_ratios=[0.8, 0.2], height_ratios=[0.2, 0.8],
                  wspace=0.05, hspace=0.05)

    # top marginal (x)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_top.fill_between(x_centers, x_marginal, alpha=0.5)
    ax_top.set_xlim((x_edges[0], x_edges[-1]))
    ax_top.set_ylabel("Count")
    ax_top.set_xticks([])

    # right marginal (y)
    ax_right = fig.add_subplot(gs[1, 1])
    ax_right.fill_betweenx(y_centers, y_marginal, alpha=0.5)
    ax_right.set_ylim((y_edges[0], y_edges[-1]))
    ax_right.set_xlabel("Count")
    ax_right.set_yticks([])

    # main heatmap
    ax_main = fig.add_subplot(gs[1, 0])
    mesh = ax_main.pcolormesh(x_edges, y_edges, H.T, shading="auto", cmap=cmap)
    ax_main.set_xlabel("x (m)")
    ax_main.set_ylabel("y (m)")
    ax_main.set_aspect("auto")

    # overlay rectangular box on main plot
    if box is not None:
        xmin, xmax = box["h_bounds"]
        ymin, ymax = box["v_bounds"]
        ax_main.axvline(xmin, color="lime", linewidth=2)
        ax_main.axvline(xmax, color="lime", linewidth=2)
        ax_main.axhline(ymin, color="lime", linewidth=2)
        ax_main.axhline(ymax, color="lime", linewidth=2)

    # inset cutout in the empty top-right cell
    if cutout is not None:
        ax_cut = fig.add_subplot(gs[0, 1])
        vmin = np.nanmin(H)
        vmax = np.nanmax(H)
        xmin, xmax = box["h_bounds"]
        ymin, ymax = box["v_bounds"]

        # imshow is convenient for arbitrary cutout sizes
        ax_cut.imshow(
            cutout["cutout"],
            extent=[xmin, xmax, ymin, ymax],
            origin="lower",
            aspect="equal",
            cmap=cmap,
            vmin=vmin, vmax=vmax,
        )
        ax_cut.set_title("Cutout")
        ax_cut.axis("off")

    plt.suptitle(title)
    plt.show()
    return fig


# --------------------------- High-level API ---------------------------

def hist(
    coords,
    nx=80, ny=80,
    xlims=(-0.2, 0.2), ylims=(-0.2, 0.2),
    make_plot=True,
    cutout=False,
    h=0.05, v=0.05,
    center=(0.0, 0.0),
    title="Beam profile at Target",
):
    """
    High-level driver:
      - If cutout=False:
          * make_plot=True  -> plot 2D + marginals; return data + figure
          * make_plot=False -> return data dict only
      - If cutout=True:
          * make_plot=False -> return cutout dict (+ some data metadata)
          * make_plot=True  -> plot 2D + marginals + box overlay + cutout inset; return everything
    """

    
    # coords should be [N_particles x N_turns x 6]
    _shape = coords.shape
    
    # Check that `coords` is 3D
    if len(_shape) != 3:
        try:
            coords = coords.reshape(int(np.prod(_shape) // 6), 1, 6)
        except Exception:
            raise ValueError("`coords` array should have shape (N_particles x 6) or (N_particles x N_turns x 6)")

    x = coords[:, -1, 0]
    y = coords[:, -1, 2]

    data = compute_hist(x, y, nx=nx, ny=ny, xlims=xlims, ylims=ylims)

    if not cutout and not make_plot:
        return data

    if cutout:
        c = compute_cutout(
            data["hist"], data["x_edges"], data["y_edges"],
            h=h, v=v, center=center
        )
        if not make_plot:
            # return cutout plus a few useful axes arrays
            return dict(
                **c,
                x_edges=data["x_edges"],
                y_edges=data["y_edges"],
                x_centers=data["x_centers"],
                y_centers=data["y_centers"],
                n_total=_shape[0]
            )
        # with plot: overlay and inset
        fig = plot_hist(
            data,
            title=title,
            box={
                "h_bounds": (center[0] - h, center[0] + h), 
                "v_bounds": (center[1] - v, center[1] + v)
            },
            cutout=c,
        )
        return dict(**data, **c, figure=fig)

    # plot, but no cutout
    fig = plot_hist(data, title=title)
    return dict(**data, figure=fig)
