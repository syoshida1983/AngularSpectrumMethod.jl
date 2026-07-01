export ScalableASM
export ScalableASM!

"""
    ScalableASM(u, λ, Δx, Δy, z; expand=true)

return automatically scaled diffraction field by the scalable ASM (see Ref. 1).
The sampling interval in the destination plane ``\\Delta_{d}`` is ``\\Delta_{d}=\\dfrac{\\lambda z}{pN\\Delta_{s}}``,
where ``\\Delta_{s}`` is the sampling interval in the source plane, ``N`` is the number of pixels in the source or destination plane, and ``p=2`` is the padding factor.
If the diffraction distance ``z`` is negative, backpropagation is performed,
i.e. diffracting at positive ``z`` and then diffracting at negative ``z``, the field returns to its original state.

> 1. [Rainer Heintzmann, Lars Loetgering, and Felix Wechsler, "Scalable angular spectrum propagation," Optica **10**, 1407-1416 (2023)](https://doi.org/10.1364/OPTICA.497809)
"""
function ScalableASM(u, λ, Δx, Δy, z; expand=true)
    Δ = [Δy, Δx]
    L = Δ.*size(u)
    R = Δ./λ
    zₗᵢₘᵢₜ  = @. L*inv(abs(inv(4R) - inv(√(16R^2 + 2))))
    zₘᵢₙ   = 4L.*R
    reduce(|, abs(z) .≥ zₗᵢₘᵢₜ) && @warn "Propagated field might be affected by vignetting"
    reduce(|, abs(z) .< zₘᵢₙ) && @warn "Magnification might be less than one"

    N   = ifelse(expand, size(u).*2, size(u))       # row and column directions are x- and y-axis, respectively
    Lₚ  = N.*Δ                                      # computational domain sizes
    ν   = fftfreq.(N, inv.(Δ))                      # spatial frequencies (DC corner)
    ν²  = @. ν[1]^2 + ν[2]'^2
    νz  = @. √(1/λ^2 - ν² + 0im)                    # spatial frequencies in the z-axis
    H   = @. exp(2π*im*z*(real(νz) - (1/λ - λ*ν²/2)))*exp(-2π*abs(z)*imag(νz))                      # transfer function
    W   = @. (abs(ν[1]/νz - λ*ν[1]) ≤ abs(Lₚ[1]/(2z)))*(abs(ν[2]'/νz - λ*ν[2]') ≤ abs(Lₚ[2]/(2z)))  # window function
    r₁  = fftfreq.(N, Lₚ)                           # coordinates in the source plane
    r₂  = fftfreq.(N, λ*z*inv.(Δ))                  # coordinates in the destination plane
    Q₁  = @. exp(π*im/(λ*z)*(r₁[1]^2 + r₁[2]'^2))   # Fresnel kernel in the source plane
    Q₂  = @. exp(2π*im*z/λ)*exp(π*im/(λ*z)*(r₂[1]^2 + r₂[2]'^2))
    û::Matrix{ComplexF64} = select_region_view(u, new_size=N)

    if ~signbit(z)  # forward
        û = ifft(fft(ifftshift(û)).*H.*W)
        û = fftshift(Q₂.*fft(Q₁.*û))/(im*√length(û))
    else            # backward
        û = Q₂.*conj(fft(conj(Q₁.*ifftshift(û))))/(-im*√length(û))
        û = fftshift(ifft(fft(û).*H.*W))
    end

    return select_region(û, new_size=size(u))
end

"""
    ScalableASM!(u, λ, Δx, Δy, z; expand=true)

Same as ScalableASM, but operates in-place on `u`, which must be an array of complex floating-point numbers.
"""
function ScalableASM!(u, λ, Δx, Δy, z; expand=true)
    u[:,:] = ScalableASM(u, λ, Δx, Δy, z; expand)
end
