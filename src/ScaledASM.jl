export ScaledASM
export ScaledASM!

"""
    ScaledASM(u, λ, Δx, Δy, z, R; expand=true)

return scaled diffraction field according to the scale factor ``R`` by the scaled ASM (see Ref. 1).

> 1. [Tomoyoshi Shimobaba, Kyoji Matsushima, Takashi Kakue, Nobuyuki Masuda, and Tomoyoshi Ito, "Scaled angular spectrum method," Opt. Lett. **37**, 4128-4130 (2012)](https://doi.org/10.1364/OL.37.004128)
"""
function ScaledASM(u, λ, Δx, Δy, z, R; expand=true)
    N = ifelse(expand, size(u).*2, size(u)) # row and column directions are x- and y-axis, respectively
    f = @. fftshift(fftfreq(N))             # DFT sample frequencies (DC-centered)
    S = ifelse(abs(R) ≤ 1, R, 1/R)          # scale factor
    k = [repeat(f[1].*S, N[2])';            # nodes for NFFT
        repeat(f[2].*S, inner=N[1])']
    J = ifelse(abs(R) ≤ 1, abs(R)/(N[1]*N[2]), 1/abs(R))    # Jacobian for energy conservation
    H = @. exp(2π*im*z*√(1/λ^2 - (f[1]/Δy)^2 - (f[2]'/Δx)^2 + 0im)) # transfer function
    û::Matrix{ComplexF64} = select_region_view(u, new_size=N)

    if abs(R) ≤ 1   # scale down
        û = ifftshift(fft(fftshift(û)))
        û = @. û*H*J
        û = reshape(conj.(nfft(k, conj.(û))::Vector{ComplexF64}), N)
    else    # scale up
        û = conj.(nfft_adjoint(k, N, conj.(û[:]))::Matrix{ComplexF64})
        û = @. û*H*J
        û = ifftshift(ifft(fftshift(û)))
    end

    return select_region_view(û, new_size=size(u))
end

"""
    ScaledASM!(u, λ, Δx, Δy, z, R; expand=true)

Same as ScaledASM, but operates in-place on u, which must be an array of complex floating-point numbers.
"""
function ScaledASM!(u, λ, Δx, Δy, z, R; expand=true)
    u[:,:] = ScaledASM(u, λ, Δx, Δy, z, R; expand)
end
