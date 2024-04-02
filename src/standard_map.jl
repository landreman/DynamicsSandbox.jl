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
    plot_fixed_point_and_eigvals(k, -0.5, 0.0)
    p
end

function plot_fixed_point_and_eigvals(k, x0, y0)
    M = standard_map_gradient(x0, y0, k)
    eigvals, eigvects = eigen(M)
    @assert prod(eigvals) ≈ 1
    ϵ = 0.01
    r = [x0, y0]
    for j in 1:2
        if eigvals[j] > 1
            color = :red
        else
            color = :blue
        end
        plot!(
            [x0, x0 + ϵ * eigvects[1, j]],
            [y0, y0 + ϵ * eigvects[2, j]],
            color=color,
            label=nothing,
            )

        # for log10δ in range(-5, -3, length=3000)
        #     δ = 10 ^ log10δ
        amplitudes = collect(range(0, log(eigvals[j]), length=100))
        for preδ in amplitudes[1:end - 1]
            δ = 1e-3 * exp(preδ)
            for eigval_sign in [-1, 1]
                plot_standard_map_trajectory(
                    k, 
                    x0 + eigval_sign * δ * eigvects[1, j],
                    y0 + eigval_sign * δ * eigvects[2, j],
                    10;
                    color=color)
            end
        end
    end
end