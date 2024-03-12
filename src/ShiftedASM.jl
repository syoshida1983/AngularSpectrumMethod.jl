export ShiftedASM
export ShiftedASM!

function CenterFrequency(x₀, S, u₊, u₋)
    if S ≤ x₀
        return (u₊ + u₋)/2
    elseif -S < x₀ < S
        return (u₊ - u₋)/2
    elseif x₀ ≤ -S
        return -(u₊ + u₋)/2
    else
        return zero(x₀)
    end
end

function BandWidth(x₀, S, u₊, u₋)
    if S ≤ x₀
        return u₊ - u₋
    elseif -S < x₀ < S
        return u₊ + u₋
    elseif x₀ ≤ -S
        return u₋ - u₊
    else
        return zero(x₀)
    end
end

"""
    ShiftedASM(u, λ, Δx, Δy, z, x₀, y₀; expand=true)

return shifted diffraction field with the shift distance ``x_{0}`` and ``y_{0}`` by the shifted ASM (see Ref. 1).

> 1. [Kyoji Matsushima, "Shifted angular spectrum method for off-axis numerical propagation," Opt. Express **18**, 18453-18463 (2010)](https://doi.org/10.1364/OE.18.018453)
"""
function ShiftedASM(u, λ, Δx, Δy, z, x₀, y₀; expand=true)
    N = ifelse(expand, size(u).*2, size(u)) # row and column directions are x- and y-axis, respectively
    r₀ = [y₀, x₀]
    S = N.*[Δy, Δx]./2                      # size of source sampling window
    ν = fftfreq.(N, inv.([Δy, Δx]))         # spatial frequencies (DC corner)
    ν₊ = @. inv(λ*√(z^2/(r₀ + S)^2 + 1))    # upper limit of bandwidth
    ν₋ = @. inv(λ*√(z^2/(r₀ - S)^2 + 1))    # lower limit of bandwidth
    ν₀ = CenterFrequency.(r₀, S, ν₊, ν₋)
    νw = BandWidth.(r₀, S, ν₊, ν₋)
    H = @. exp(2π*im*((r₀[1]*ν[1] + r₀[2]*ν[2]') + z*√(1/λ^2 - ν[1]^2 - ν[2]'^2 + 0im)))    # transfer function
    W = @. rect((ν[1] - ν₀[1])/νw[1])*rect((ν[2]' - ν₀[2])/νw[2])   # window function
    ũ = select_region_view(u, new_size=N)
    û = fftshift(ifft(fft(ifftshift(ũ)).*H.*W))

    return select_region_view(û, new_size=size(u))
end

"""
    ShiftedASM!(u, λ, Δx, Δy, z, x₀, y₀; expand=true)

Same as ShiftedASM, but operates in-place on u, which must be an array of complex floating-point numbers.
"""
function ShiftedASM!(u, λ, Δx, Δy, z, x₀, y₀; expand=true)
    u[:,:] = ShiftedASM(u, λ, Δx, Δy, z, x₀, y₀; expand)
end
