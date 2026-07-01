export BandLimitedASM
export BandLimitedASM!

"""
    BandLimitedASM(u, λ, Δx, Δy, z; expand=true)

return diffracted field by the band-limited ASM (see Ref. 1).

> 1. [Kyoji Matsushima and Tomoyoshi Shimobaba, "Band-Limited Angular Spectrum Method for Numerical Simulation of Free-Space Propagation in Far and Near Fields," Opt. Express **17**, 19662-19673 (2009)](https://doi.org/10.1364/OE.17.019662)
"""
function BandLimitedASM(u, λ, Δx, Δy, z; expand=true)
    N   = ifelse(expand, size(u).*2, size(u))                           # row and column directions are x- and y-axis, respectively
    ν   = fftfreq.(N, inv.([Δy, Δx]))                                   # spatial frequencies (DC corner)
    νz  = @. √(1/λ^2 - ν[1]^2 - ν[2]'^2 + 0im)                          # spatial frequencies in the z-axis
    Δν  = inv.(N.*[Δy, Δx])                                             # sampling interval
    νₗᵢₘᵢₜ  = @. 1/(λ*√((2Δν*z)^2 + 1))                                 # bandwidth limit
    H   = @. exp(2π*im*z*real(νz))*exp(-2π*abs(z)*imag(νz))             # transfer function
    W   = @. (abs(ν[1]/(2νₗᵢₘᵢₜ[1])) < 1/2)*(abs(ν[2]'/(2νₗᵢₘᵢₜ[2])) < 1/2) # window function
    ũ   = select_region_view(u, new_size=N)
    û   = fftshift(ifft(fft(ifftshift(ũ)).*H.*W))

    return select_region(û, new_size=size(u))
end

"""
    BandLimitedASM!(u, λ, Δx, Δy, z; expand=true)

Same as BandLimitedASM, but operates in-place on `u`, which must be an array of complex floating-point numbers.
"""
function BandLimitedASM!(u, λ, Δx, Δy, z; expand=true)
    u[:,:] = BandLimitedASM(u, λ, Δx, Δy, z; expand)
end
