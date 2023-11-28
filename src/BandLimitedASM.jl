export BandLimitedASM

"""
    BandLimitedASM(u, λ, Δx, Δy, z; expand=true)

return diffracted field by the band-limited ASM (see Ref. 1).

> 1. [Kyoji Matsushima and Tomoyoshi Shimobaba, "Band-Limited Angular Spectrum Method for Numerical Simulation of Free-Space Propagation in Far and Near Fields," Opt. Express **17**, 19662-19673 (2009)](https://doi.org/10.1364/OE.17.019662)
"""
function BandLimitedASM(u, λ, Δx, Δy, z; expand=true)
    û = ifelse(expand, ifftshift(fft(fftshift(padzeros(u)))), ifftshift(fft(fftshift(u))))
    Ny, Nx = size(û)                # row and column directions are x- and y-axis, respectively
    Δv = [1/(Ny*Δy), 1/(Nx*Δx)]     # sampling interval
    vₗ = @. 1/(λ*√((2Δv*z)^2 + 1))  # bandwidth limit

    @fastmath @inbounds for j ∈ 1:Nx, i ∈ 1:Ny
        v = [(i - 1 - Ny/2)/(Ny*Δy), (j - 1 - Nx/2)/(Nx*Δx)]    # spatial frequency u, v
        w = √(1/λ^2 - v⋅v + 0im)                                # spatial frequency w
        û[i, j] *= exp(2π*im*z*w)*prod(rect.(v./(2 .*vₗ)))
    end

    return ifelse(expand, crop(ifftshift(ifft(fftshift(û)))), ifftshift(ifft(fftshift(û))))
end
