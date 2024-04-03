function standard_map(x, y, k)
    yprime = mod(y - (k / (2π)) * sin(2π * x), 1)
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
            n_repetitions_of_eigenvalue = 7
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

function find_stable_unstable_manifold_intersection(k)
    x0 = 0.5
    y0 = 0

    # Initial guesses for the λ parameters. The precise values are not critical, but they
    # should be small compared to 1.
    λs = 1e-3
    λu = λs
    # At the end of variable names, "s" = stable, "u" = unstable.

    # Rough estimate of the distance between the adjacent X-points:
    estimated_distance = 1.0

    # Find the directions by which the (un)stable manifolds move out from the X point:
    M = standard_map_gradient(x0, y0, k)
    eigvals, eigvects = eigen(M)
    if eigvals[1] < eigvals[2]
        vs = eigvects[:, 1]
        vu = eigvects[:, 2]
    else
        vu = eigvects[:, 1]
        vs = eigvects[:, 2]
    end

    # Determine how many iterates of the map are needed:
    ru = [x0, y0] + λu * vu
    rs = [x0, y0] + λs * vs
    N = 0
    last_distance = 1e100
    while true
        ru = standard_map(ru..., k)
        rs = standard_map_inverse(rs..., k)
        N += 1
        dr = shift_to_square_around_origin.(ru) - shift_to_square_around_origin.(rs)
        distance = norm(dr)
        # @show N, ru, rs, dr, distance
        @show N, distance
        # # Find the distance between ru and rs, considering periodicity of the domain:
        # distance = min(
        #     norm(dr + [-1, -1]),
        #     norm(dr + [0, -1]),
        #     norm(dr + [1, -1]),
        #     norm(dr + [-1, 0]),
        #     norm(dr + [0, 0]),
        #     norm(dr + [1, 0]),
        #     norm(dr + [-1, 1]),
        #     norm(dr + [0, 1]),
        #     norm(dr + [1, 1]),
        #     )
        #@show N, distance
        if distance > last_distance
            break
        end
        last_distance = distance
    end

    function λ_to_r(λu, λs)
        ru = [x0, y0] + λu * vu
        rs = [x0, y0] + λs * vs
        # Interate the forward and backward maps N times:
        for j in 1:N
            ru = standard_map(ru..., k)
            rs = standard_map_inverse(rs..., k)
        end
        return ru, rs
    end

    function residual(λ)
        λu, λs = λ
        ru, rs = λ_to_r(λu, λs)
        @show ru, rs
        return shift_to_square_around_origin.(ru) - shift_to_square_around_origin.(rs)
    end

    # Use Newton's method with finite differences to solve for the λ parameters:
    println("Solving for 1st homoclinic point")
    solution = nlsolve(residual, [λu, λs])
    @show solution
    λ1u, λ1s = solution.zero
    ru, rs = λ_to_r(λ1u, λ1s)
    scatter!([shift_to_square_around_origin(ru[1])], [shift_to_square_around_origin(ru[2])], color=:red, markershape=:cross, label=nothing, markersize=6, markerstrokewidth=2)
    scatter!([shift_to_square_around_origin(rs[1])], [shift_to_square_around_origin(rs[2])], color=:blue, markershape=:xcross, label=nothing, markersize=6, markerstrokewidth=2)

    # Initial guess for the 2nd homoclinic point:
    factor = sqrt(eigvals[1])
    λ2u = λ1u * factor
    λ2s = λ1s / factor
    println("Solving for 2nd homoclinic point")
    solution = nlsolve(residual, [λ2u, λ2s])
    @show solution
    λ2u, λ2s = solution.zero
    ru, rs = λ_to_r(λ2u, λ2s)
    scatter!([shift_to_square_around_origin(ru[1])], [shift_to_square_around_origin(ru[2])], color=:magenta, markershape=:cross, label=nothing, markersize=6, markerstrokewidth=2)
    scatter!([shift_to_square_around_origin(rs[1])], [shift_to_square_around_origin(rs[2])], color=:darkgreen, markershape=:xcross, label=nothing, markersize=6, markerstrokewidth=2)
end