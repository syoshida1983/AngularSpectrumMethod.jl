export ASM

"""
    ASM(u, λ, Δx, Δy, z; expand=true)

return diffracted field by the angular spectrum method (ASM).
Evanescent waves are not eliminated but attenuated as ``\\exp(-2{\\pi}wz)``.
Without attenuation, the total energy ``\\iint|u|\\mathrm{d}x\\mathrm{d}y`` is conserved.

# Arguments
- `u`: input field.
- `λ`: wavelength.
- `Δx`: sampling interval the in x-axis.
- `Δy`: sampling interval the in y-axis.
- `z`: diffraction distance.
- `expand=true`: if true (default), perform 4× expansion and zero padding for aliasing suppression.

!!! note
    The x-axis is the horizontal direction, and the y-axis is the vertical.
"""
function ASM(u, λ, Δx, Δy, z; expand=true)
    û = ifelse(expand, ifftshift(fft(fftshift(padzeros(u)))), ifftshift(fft(fftshift(u))))
    Ny, Nx = size(û)    # row and column directions are x- and y-axis, respectively

    @fastmath @inbounds for j ∈ 1:Nx, i ∈ 1:Ny
        v = [(i - 1 - Ny/2)/(Ny*Δy), (j - 1 - Nx/2)/(Nx*Δx)]    # spatial frequency u, v
        w = √(1/λ^2 - v⋅v + 0im)                                # spatial frequency w
        û[i, j] *= exp(2π*im*z*w)
    end

    return ifelse(expand, crop(ifftshift(ifft(fftshift(û)))), ifftshift(ifft(fftshift(û))))
end
