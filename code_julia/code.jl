# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using LinearAlgebra
using SparseArrays

using LaTeXStrings
using CairoMakie
set_theme!(theme_latexfonts();
           fontsize = 26,
           linewidth = 3,
           markersize = 16,
           Lines = (cycle = Cycle([:color, :linestyle], covary = true),),
           Scatter = (cycle = Cycle([:color, :marker], covary = true),))



#####################################################################
# Time integration

function runge_kutta(f!, dt, tspan, u_analytical, parameters;
                     relaxation = true, relax_last_step = true,
                     dot_entropy = dot)
    # Runge-Kutta parameters
    # RDPK3SpFSAL35
    A = [
        0.0 0.0 0.0 0.0 0.0
        0.23002986245180762 0.0 0.0 0.0 0.0
        0.1028611905146702 0.30214341669482886 0.0 0.0 0.0
        0.06278872625842052 0.029432944569291738 0.8025606185416311 0.0 0.0
        0.047715578851412185 0.03637030463243372 0.20321191503846886 0.4362158943603441 0.0]
    b = [
        0.11479359710235411
        0.08933442853113313
        0.43558710250086163
        0.2473576188201451
        0.11292725304550591
    ]
    c = sum(A, dims = 2) |> vec
    s = length(b)

    # Values for the Runge-Kutta time stepping
    t = first(tspan)
    u = u_analytical(t, parameters)
    k = [similar(u) for _ in 1:s]
    utmp = similar(u)
    unew = similar(u)

    # Results gathered during the computation
    times = Vector{typeof(t)}()
    entropies = Vector{eltype(u)}()

    while t < last(tspan)
        # Avoid stepping over the final time
        if t + dt > last(tspan)
            dt = last(tspan) - t
        end

        # For relaxation
        entropy_diff = zero(t)

        # Runge-Kutta step
        for i in 1:s
            fill!(utmp, 0)
            for j in 1:(i - 1)
                @. utmp = utmp + A[i, j] * k[j]
            end
            @. utmp = u + dt * utmp
            f!(k[i], utmp, parameters, t + c[i] * dt)

            if relaxation
                entropy_rate = 2 * dot_entropy(utmp, k[i])
                entropy_diff += b[i] * dt * entropy_rate
            end
        end

        fill!(unew, 0)
        for i in 1:s
            @. unew = unew + b[i] * k[i]
        end
        @. unew = u + dt * unew
        tnew = t + dt

        if relaxation
            # For inner product norms, we have
            #   gamma = (entropy_new - entropy_old - 2 * <u, unew - u>) / |unew - u|^2
            @. utmp = unew - u
            gamma = (entropy_diff - 2 * dot_entropy(u, utmp)) / dot_entropy(utmp, utmp)
        else
            gamma = one(t)
        end

        if tnew != last(tspan)
            @. u = u + gamma * (unew - u)
            t = t + gamma * (tnew - t)
        else
            if relax_last_step
                # Use an IDP step in the last step to hit the final time
                @. u = u + gamma * (unew - u)
            else
                @. u = unew
            end
            t = tnew
        end

        # compute functionals
        push!(times, t)
        push!(entropies, dot_entropy(u, u))
    end

    return (; t, u, times, entropies, parameters)
end



#####################################################################
# Setup active flux matrices

function central_matrices(xmin, xmax, n_volumes)
    @assert n_volumes > 2
    n = 2 * n_volumes
    Δx = (xmax - xmin) / n_volumes

    D = spzeros(n, n)
    for i in 1:n_volumes
        # derivative of the volumes
        D[2i, 2i-1] = -1
        if i == n_volumes
            D[2i, 1] = 1
        else
            D[2i, 2i+1] = 1
        end

        # derivative of the point values
        if i == 1
            D[2i-1, 2n_volumes-1] = 1
            D[2i-1, 2n_volumes] = -3
            D[2i-1, 2i+0] = 3
            D[2i-1, 2i+1] = -1
        elseif i == n_volumes
            D[2i-1, 2i-3] = 1
            D[2i-1, 2i-2] = -3
            D[2i-1, 2i+0] = 3
            D[2i-1, 1] = -1
        else
            D[2i-1, 2i-3] = 1
            D[2i-1, 2i-2] = -3
            D[2i-1, 2i+0] = 3
            D[2i-1, 2i+1] = -1
        end
    end
    @. D /= Δx

    # Diagonal mass matrix
    m_v = 3 * Δx / 4
    m_p = Δx / 4
    m = zeros(n)
    for i in 1:n_volumes
        m[2i] = m_v
        m[2i-1] = m_p
    end
    M = Diagonal(m)

    M, D
end

function upwind_matrices(xmin, xmax, n_volumes)
    @assert n_volumes > 2
    n = 2 * n_volumes
    Δx = (xmax - xmin) / n_volumes

    Dm = spzeros(n, n)
    Dp = spzeros(n, n)
    for i in 1:n_volumes
        # derivative of the volumes
        Dm[2i, 2i-1] = -1
        Dp[2i, 2i-1] = -1
        if i == n_volumes
            Dm[2i, 1] = 1
            Dp[2i, 1] = 1
        else
            Dm[2i, 2i+1] = 1
            Dp[2i, 2i+1] = 1
        end

        # derivative of the point values
        if i == 1
            Dm[2i-1, 2n_volumes-1] = 2
            Dm[2i-1, 2n_volumes] = -6
            Dm[2i-1, 2i-1] = 4
            Dp[2i-1, 2i-1] = -4
            Dp[2i-1, 2i+0] = 6
            Dp[2i-1, 2i+1] = -2
        elseif i == n_volumes
            Dm[2i-1, 2i-3] = 2
            Dm[2i-1, 2i-2] = -6
            Dm[2i-1, 2i-1] = 4
            Dp[2i-1, 2i-1] = -4
            Dp[2i-1, 2i+0] = 6
            Dp[2i-1, 1] = -2
        else
            Dm[2i-1, 2i-3] = 2
            Dm[2i-1, 2i-2] = -6
            Dm[2i-1, 2i-1] = 4
            Dp[2i-1, 2i-1] = -4
            Dp[2i-1, 2i+0] = 6
            Dp[2i-1, 2i+1] = -2
        end
    end
    @. Dm /= Δx
    @. Dp /= Δx

    # Banded mass matrix
    M = spzeros(n, n)
    m_v = Δx
    m_p = 2 * m_v / 3
    m_vp = -m_v / 2
    m_pv = m_vp
    # m_vv = 0
    m_pp = m_v / 6
    for i in 1:n_volumes
        M[2i, 2i] = m_v
        M[2i-1, 2i-1] = m_p

        M[2i, 2i-1] = m_vp
        M[2i-1, 2i] = m_pv
        if i < n_volumes
            M[2i, 2i+1] = m_vp
            M[2i+1, 2i] = m_pv
        else
            M[2i, 1] = m_vp
            M[1, 2i] = m_pv
        end

        if i == 1
            # M[2i, 2n_volumes] = m_vv
            # M[2i, 2i+2] = m_vv
            M[2i-1, 2i+1] = m_pp
            M[2i-1, 2n_volumes-1] = m_pp
        elseif i == n_volumes
            # M[2i, 2i-2] = m_vv
            # M[2i, 2] = m_vv
            M[2i-1, 2i-3] = m_pp
            M[2i-1, 1] = m_pp
        else
            # M[2i, 2i-2] = m_vv
            # M[2i, 2i+2] = m_vv
            M[2i-1, 2i+1] = m_pp
            M[2i-1, 2i-3] = m_pp
        end
    end

    return M, Dm, Dp
end



#####################################################################
# Perform numerical experiments
function rhs!(du, u, parameters, t)
    mul!(du, parameters.D, u, -1, 0)
    return nothing
end

function plot_energies()
    xmin = 0.0
    xmax = 2 * pi
    n_volumes = 50
    x = range(xmin, xmax, length = 2 * n_volumes + 1)[begin:(end - 1)]
    tspan = (0.0, 80.0)
    dt = 0.5 * (xmax - xmin) / n_volumes
    u_analytical(t, parameters) = @. exp(sin(t - parameters.x))

    fig = Figure(size = (2 * 600, 450))

    ax1 = Axis(fig[1, 1];
               xlabel = L"Time $t$", ylabel = "Energy Change",
               title = "Central Version")
    ax2 = Axis(fig[1, 2];
               xlabel = L"Time $t$", ylabel = "Energy Change",
               title = "Upwind Version")

    let ax = ax1
        M, D = central_matrices(xmin, xmax, n_volumes)
        parameters = (; M = M, D = D, x = x)
        results = runge_kutta(rhs!, dt, tspan, u_analytical, parameters;
                              relaxation = true,
                              relax_last_step = true,
                              dot_entropy = (u, v) -> dot(u, M, v))
        lines!(ax, results.times, results.entropies .- first(results.entropies))
    end

    let ax = ax2
        M, Dm, Dp = upwind_matrices(xmin, xmax, n_volumes)
        parameters = (; M = M, D = Dm, x = x)
        results = runge_kutta(rhs!, dt, tspan, u_analytical, parameters;
                              relaxation = true,
                              relax_last_step = true,
                              dot_entropy = (u, v) -> dot(u, M, v))
        lines!(ax, results.times, results.entropies .- first(results.entropies))
    end

    filename = joinpath(@__DIR__, "energies.pdf")
    save(filename, fig)
    @info "results saved" filename

    return fig
end

