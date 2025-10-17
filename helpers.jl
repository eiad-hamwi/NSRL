using Plots
using StatsBase
using Printf

# Ensure you have Plots.jl and StatsBase.jl installed:
# using Pkg
# Pkg.add("Plots")
# Pkg.add("StatsBase")

export compute_hist, compute_cutout, plot_hist, hist

"""
    compute_hist(x, y; nx=80, ny=80, xlims=(-0.2, 0.2), ylims=(-0.2, 0.2))

Compute a NaN-masked 2D histogram and 1D marginals.
Returns a Dict with hist, edges, centers, and marginals.
"""
function compute_hist(x, y; nx=80, ny=80, xlims=nothing, ylims=nothing)
    # Determine limits if not provided
    _xlims = isnothing(xlims) ? (minimum(x), maximum(x)) : xlims
    _ylims = isnothing(ylims) ? (minimum(y), maximum(y)) : ylims

    x_edges = range(_xlims[1], _xlims[2], length=nx + 1)
    y_edges = range(_ylims[1], _ylims[2], length=ny + 1)

    # Compute histogram
    h = fit(Histogram, (x, y), (x_edges, y_edges))
    H = Float64.(h.weights)
    H[H .== 0] .= NaN  # Mask empty bins

    # Compute centers and marginals
    x_centers = 0.5 .* (x_edges[1:end-1] .+ x_edges[2:end])
    y_centers = 0.5 .* (y_edges[1:end-1] .+ y_edges[2:end])
    x_marginal = vec(sum(h.weights, dims=2))
    y_marginal = vec(sum(h.weights, dims=1))

    return Dict(
        "hist" => H,
        "x_edges" => x_edges, "y_edges" => y_edges,
        "x_centers" => x_centers, "y_centers" => y_centers,
        "x_marginal" => x_marginal, "y_marginal" => y_marginal,
    )
end

"""
    compute_cutout(H, x_edges, y_edges; h=0.05, v=0.05, center=(0.0, 0.0))

Extract a zero-NaN-preserving rectangular cutout of H around `center`
with half-widths (h, v). Returns the cutout array plus metadata.
"""
function compute_cutout(H, x_edges, y_edges; h=0.05, v=0.05, center=(0.0, 0.0))
    x0, y0 = center

    # Step sizes
    h_step = x_edges[2] - x_edges[1]
    v_step = y_edges[2] - y_edges[1]

    # Index of center in bin space
    h_center_idx = round(Int, (x0 - x_edges[1]) / h_step) + 1
    v_center_idx = round(Int, (y0 - y_edges[1]) / v_step) + 1

    half_h_bins = round(Int, h / h_step)
    half_v_bins = round(Int, v / v_step)

    # Clamp to array bounds (H is shape (nx, ny))
    h_start = max(h_center_idx - half_h_bins, 1)
    h_stop  = min(h_center_idx + half_h_bins - 1, size(H, 1))
    v_start = max(v_center_idx - half_v_bins, 1)
    v_stop  = min(v_center_idx + half_v_bins - 1, size(H, 2))

    # Transpose H for image-like use (rows=y, cols=x)
    cut = H[h_start:h_stop, v_start:v_stop]'

    # Mask NaN bins
    cut[isnan.(cut)] .= 0.

    return Dict(
        "cutout" => cut,
        "center" => center,
    )
end

"""
    plot_hist(data; title="Beam profile at Target", box=nothing, cutout=nothing, cmap=:inferno)

Plot 2D histogram with marginals.
Optionally overlay a rectangular box and/or add a cutout inset.
"""
function plot_hist(data; title="Beam profile at Target", box=nothing, cutout=nothing, cmap=:inferno)
    H = data["hist"]
    x_edges, y_edges = data["x_edges"], data["y_edges"]
    x_centers, y_centers = data["x_centers"], data["y_centers"]
    x_marginal, y_marginal = data["x_marginal"], data["y_marginal"]

    # Main heatmap
    p_main = heatmap(x_centers, y_centers, H',
        xlabel="x (m)", ylabel="y (m)",
        colorbar=false, aspect_ratio=:equal,
        cmap=cmap
    )

    # Top marginal (x)
    p_top = plot(x_centers, x_marginal,
        seriestype=:line, fillrange=0, fillalpha=0.5,
        legend=false, grid=false,
        xlims=(x_edges[1], x_edges[end]),
        ylabel="Count",
        xticks=false
    )

    # Right marginal (y)
    p_right = plot(y_marginal, y_centers,
        fillrange=0, fillalpha=0.5,
        legend=false, grid=false,
        ylims=(y_edges[1], y_edges[end]),
        xlabel="Count",
        yticks=false,
    )

    # Overlay rectangular box on main plot
    if !isnothing(box)
        xmin, xmax = box["h_bounds"]
        ymin, ymax = box["v_bounds"]
        plot!(p_main, [xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin],
              linecolor=:lime, linewidth=3, label=nothing)
    end
    
    # Inset for cutout in the empty top-right cell
    if !isnothing(cutout)
        vmin, vmax = extrema(filter(!isnan, H))
        xmin, xmax = box["h_bounds"]
        ymin, ymax = box["v_bounds"]

        p_cut = heatmap(
            range(xmin, stop=xmax, length=size(cutout["cutout"], 2)),
            range(ymin, stop=ymax, length=size(cutout["cutout"], 1)),
            cutout["cutout"],
            title="Cutout",
            aspect_ratio=:equal,
            cmap=cmap,
            colorbar=false,
            clims=(vmin, vmax),
            framestyle=:none,
        )
    else
        p_cut = plot(framestyle=:none)
    end
    
    # Compose layout
    l = @layout [
        a{0.2h} b{0.2w}
        c       d{0.2w}
    ]
    
    return plot(p_top, p_cut, p_main, p_right, layout=l, size=(850, 850), plot_title=title)
end


"""
    hist(bunch; nx=80, ny=80, xlims=(-0.2, 0.2), ylims=(-0.2, 0.2),
         make_plot=true, cutout=false, h=0.05, v=0.05, center=(0.0, 0.0),
         title="Beam profile at Target")

High-level driver for creating beam profile histograms and plots.
"""
function hist(bunch;
              nx=80, ny=80,
              xlims=(-0.2, 0.2), ylims=(-0.2, 0.2),
              make_plot=true,
              cutout=false,
              h=0.05, v=0.05,
              center=(0.0, 0.0),
              title="Beam profile at Target",
              )

    # Assumes bunch is a struct/Dict with bunch.coords
    x = bunch.coords.v[:, 1]
    y = bunch.coords.v[:, 3]

    data = compute_hist(x, y; nx=nx, ny=ny, xlims=xlims, ylims=ylims)

    if !cutout && !make_plot
        return data
    end

    if cutout
        c = compute_cutout(data["hist"], data["x_edges"], data["y_edges"];
                           h=h, v=v, center=center)
        
        if !make_plot
            # Return cutout plus some useful metadata
            return merge(c, Dict(
                "x_edges" => data["x_edges"],
                "y_edges" => data["y_edges"],
                "x_centers" => data["x_centers"],
                "y_centers" => data["y_centers"],
                "n_total" => length(x)
            ))
        end

        # With plot: overlay and inset
        box_bounds = Dict(
            "h_bounds" => (center[1] - h, center[1] + h),
            "v_bounds" => (center[2] - v, center[2] + v)
        )
        
        fig = plot_hist(data, title=title, box=box_bounds, cutout=c)
        return merge(data, c, Dict("figure" => fig))
    end
    
    # Plot, but no cutout
    fig = plot_hist(data, title=title)
    return merge(data, Dict("figure" => fig))
end