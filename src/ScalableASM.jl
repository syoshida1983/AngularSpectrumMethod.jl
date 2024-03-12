export ScalableASM
export ScalableASM!

function BandLimit(λ, z, L, ν, νz)
    return ifelse(abs(ν/νz - λ*ν) ≤ L/(2z), one(λ), zero(λ))
end

function DistanceLimit(λ, L, N)
    R = L/N/λ
    return L*inv(2*abs(inv(4R) - inv(√(16R^2 + 2))))
end

function MinimumDistance(λ, L, N)
    Δ = L/N
    return L*Δ/λ
end

"""
    ScalableASM(u, λ, Δx, Δy, z; expand=true)

return automatically scaled diffraction field by the scaled scalable ASM (see Ref. 1).
The sampling pitch in the destination plane ``\\Delta_{d}`` is ``\\Delta_{d}=\\frac{\\lambda z}{pN\\Delta_{s}}``,
where ``\\Delta_{s}`` is the sampling pitch in the source plane, ``N`` is the number of pixels in the source or destination plane, and ``p=2`` is the padding factor.

> 1. [Rainer Heintzmann, Lars Loetgering, and Felix Wechsler, "Scalable angular spectrum propagation," Optica **10**, 1407-1416 (2023)](https://doi.org/10.1364/OPTICA.497809)
"""
function ScalableASM(u, λ, Δx, Δy, z; expand=true)
    N = ifelse(expand, size(u).*2, size(u)) # row and column directions are x- and y-axis, respectively
    L = N.*[Δy, Δx]                         # computational domain sizes
    z ≥ minimum(DistanceLimit.(λ, L, N)) &&  @warn "Propagated field might be affected by vignetting"
    z ≤ maximum(MinimumDistance.(λ, L, N)) &&  @warn "Magnification might be less than one"
    ν = fftfreq.(N, inv.([Δy, Δx]))         # spatial frequencies (DC corner)
    ν² = @. ν[1]^2 + ν[2]'^2
    νz = @. √(1/λ^2 - ν² + 0im)
    H = @. exp(2π*im*z*(νz - (1/λ - λ*ν²/2)))   # transfer function
    W = @. BandLimit(λ, z, L[1], ν[1], νz)*BandLimit(λ, z, L[2], ν[2]', νz) # window function
    ũ = select_region_view(u, new_size=N)
    û = ifft(fft(ifftshift(ũ)).*H.*W)
    r₁ = fftfreq.(N, L)                             # coordinates in the source plane
    r₂ = fftfreq.(N, λ*z*inv.([Δy, Δx]))            # coordinates in the destination plane
    Q₁ = @. exp(π*im/(λ*z)*(r₁[1]^2 + r₁[2]'^2))    # Fresnel kernel in the source plane
    Q₂ = @. exp(2π*im*z/λ)/(im*λ*z)*exp(π*im/(λ*z)*(r₂[1]^2 + r₂[2]'^2))    # Fresnel kernel in the destination plane
    û = fftshift(Q₂.*fft(û.*Q₁))

    return select_region_view(û, new_size=size(u))
end

"""
    ScalableASM!(u, λ, Δx, Δy, z; expand=true)

Same as ScalableASM, but operates in-place on u, which must be an array of complex floating-point numbers.
"""
function ScalableASM!(u, λ, Δx, Δy, z; expand=true)
    u[:,:] = ScalableASM(u, λ, Δx, Δy, z; expand)
end
