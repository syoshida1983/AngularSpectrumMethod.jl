export ShiftedASM

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
    û = ifelse(expand, ifftshift(fft(fftshift(padzeros(u)))), ifftshift(fft(fftshift(u))))
    N₁, N₂ = size(û)
    r₀ = [x₀, y₀]
    S = [N₁*Δx/2, N₂*Δy/2]              # size of source sampling window
    v₊ = @. 1/(λ*√(z^2/(r₀ + S)^2 + 1)) # upper limit of bandwidth
    v₋ = @. 1/(λ*√(z^2/(r₀ - S)^2 + 1)) # lower limit of bandwidth
    v₀ = CenterFrequency.(r₀, S, v₊, v₋)
    vw = BandWidth.(r₀, S, v₊, v₋)

    @fastmath @inbounds for j ∈ 1:N₂, i ∈ 1:N₁
        v = [(i - 1 - N₁/2)/(N₁*Δx), (j - 1 - N₂/2)/(N₂*Δy)]    # spatial frequency u, v
        w = √(1/λ^2 - v⋅v + 0im)                                # spatial frequency w
        û[i, j] *= exp(2π*im*(r₀⋅v + z*w))*prod(rect.((v .- v₀)./vw))
    end

    return ifelse(expand, crop(ifftshift(ifft(fftshift(û)))), ifftshift(ifft(fftshift(û))))
end
