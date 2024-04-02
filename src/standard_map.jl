function standard_map(x, y, k)
    yprime = mod(y - (k / (2π)) * sin(2π * x), 1)
    #yprime = y - (k / (2π)) * sin(2π * x)
    return [mod(x + yprime, 1), yprime]
end

function standard_map_inverse(xp, yp, k)
    x = xp - yp
    y = yp + (k / (2π)) * sin(2π * x)
    return [mod(x, 1), mod(y, 1)]
end

function standard_map_gradient(x, y, k)
    sinx, cosx = sincos(2π * x)
    return [[(1 - k * cosx) (1)]
            [(-k * cosx) (1)]]
end

function shift_to_square_around_origin(x)
    x = mod(x, 1)
    if x > 0.5
        return x - 1
    else
        return x
    end
end

function plot_standard_map_trajectory(k, x0, y0, n; color=nothing)
    data = zeros(n, 2)
    data[1, 1] = x0
    data[1, 2] = y0
    for i in 2:n
        data[i, :] .= standard_map(data[i - 1, 1], data[i - 1, 2], k)
    end

    if isnothing(color)
        color = RGB(rand(3)...)
    end
    scatter!(
        shift_to_square_around_origin.(data[:, 1]),
        shift_to_square_around_origin.(data[:, 2]), 
        label=nothing, 
        markerstrokewidth=0,
        markersize=1,
        color=color,
        )
end

function plot_standard_map_Poincare(k, n, nx, ny)
    Random.seed!(0)
    x0s = range(-0.5, 0.5, length=nx)
    y0s = range(-0.5, 0.5, length=ny)
    s = 600
    p = plot(size=(s, s))
    for x0 in x0s
        for y0 in y0s
            plot_standard_map_trajectory(k, x0, y0, n)
        end
    end
    plot_stable_and_unstable_manifolds(k, -0.5, 0.0)
    p
end

function plot_stable_and_unstable_manifolds(k, x0, y0)
    M = standard_map_gradient(x0, y0, k)
    eigvals, eigvects = eigen(M)
    @assert prod(eigvals) ≈ 1
    ϵ = 0.01
    r = [x0, y0]
    for j in 1:2
        if eigvals[j] > 1
            color = :red
            map_function = standard_map
        else
            color = :blue
            map_function = standard_map_inverse
        end
        plot!(
            [x0, x0 + ϵ * eigvects[1, j]],
            [y0, y0 + ϵ * eigvects[2, j]],
            color=color,
            label=nothing,
            )

        for eigval_sign in [-1, 1]

            # Set how many points to plot along the manifolds and how far to trace them:
            n_points_per_eigenvalue = 50
            n_repetitions_of_eigenvalue = 10
            n_points_per_manifold = n_points_per_eigenvalue * n_repetitions_of_eigenvalue
            data = zeros(n_points_per_manifold, 2)

            # The first set of n_points_per_eigenvalue points are given by multiples of the eigenvector:
            initial_amplitude = 1e-2
            for i in 1:n_points_per_eigenvalue
                amplitude = eigval_sign * initial_amplitude * exp(log(eigvals[j]) * (i - 1) / n_points_per_eigenvalue)
                data[i, :] .= [x0, y0] + amplitude * eigvects[:, j]
            end

            # To get the remaining points, iterate the map:
            for i in (n_points_per_eigenvalue + 1):n_points_per_manifold
                data[i, :] = map_function(data[i - n_points_per_eigenvalue, :]..., k)
            end

            scatter!(
                shift_to_square_around_origin.(data[:, 1]),
                shift_to_square_around_origin.(data[:, 2]),
                color=color,
                label=nothing,
                markerstrokewidth=0,
                markersize=1,
                )

            # Connect the dots, if they don't cross the edges:
            threshold = 0.1
            for i in 1:(n_points_per_manifold - 1)
                x = shift_to_square_around_origin.(data[i:(i+1), 1])
                y = shift_to_square_around_origin.(data[i:(i+1), 2])
                if abs(x[2] - x[1]) < threshold && abs(y[2] - y[1]) < threshold
                    plot!(x, y, color=color, label=nothing)
                end
            end
        end

    end
end