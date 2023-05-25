export ScaledASM

"""
    ScaledASM(u, λ, Δx, Δy, z, R; expand=true)

return scaled diffraction field according to the scale factor ``R`` by the scaled ASM (see Ref. 1).

> 1. [Tomoyoshi Shimobaba, Kyoji Matsushima, Takashi Kakue, Nobuyuki Masuda, and Tomoyoshi Ito, "Scaled angular spectrum method," Opt. Lett. **37**, 4128-4130 (2012)](https://doi.org/10.1364/OL.37.004128)
"""
function ScaledASM(u, λ, Δx, Δy, z, R; expand=true)
    û::Matrix{ComplexF64} = ifelse(expand, padzeros(u), u)
    N₁, N₂ = size(û)
    S = ifelse(abs(R) ≤ 1, R, 1/R)                                              # scale factor
    k = vcat(                                                                   # nodes for NFFT
        repeat(Vector(range(-1/2, 1/2, length=N₁)).*S, N₂)',                        # x
        repeat(Vector(range(-1/2, 1/2, length=N₂)).*S, inner=N₁)')                  # y
    J = ifelse(abs(R) ≤ 1, abs(R)/(N₁*N₂), 1/abs(R))                            # Jacobian for energy conservation
    û = ifelse(abs(R) ≤ 1,                                                      # F[u]
            ifftshift(fft(fftshift(û))),                                            # scale down
            conj.(nfft_adjoint(k, size(û), conj.(û[:]))::Matrix{ComplexF64}))       # scale up

    @fastmath @inbounds for j ∈ 1:N₂, i ∈ 1:N₁
        v = [(i - 1 - N₁/2)/(N₁*Δx), (j - 1 - N₂/2)/(N₂*Δy)]    # spatial frequency u, v
        w = √(1/λ^2 - v⋅v + 0im)                                # spatial frequency w
        û[i, j] *= J*exp(2π*im*z*w)
    end

    return ifelse(expand,
            crop(ifelse(abs(R) ≤ 1,                                             # with 4× expansion
                reshape(conj.(nfft(k, conj.(û))::Vector{ComplexF64}), size(û)),     # scale down
                ifftshift(ifft(fftshift(û))))),                                     # scale up
            ifelse(abs(R) ≤ 1,                                                  # without expansion
                reshape(conj.(nfft(k, conj.(û))::Vector{ComplexF64}), size(û)),     # scale down
                ifftshift(ifft(fftshift(û)))))                                      # scale up
end
