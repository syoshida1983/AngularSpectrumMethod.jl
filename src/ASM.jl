export ASM
export ASM!

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
    N = ifelse(expand, size(u).*2, size(u)) # row and column directions are x- and y-axis, respectively
    ν = fftfreq.(N, inv.([Δy, Δx]))         # spatial frequencies (DC corner)
    H = @. exp(2π*im*z*√(1/λ^2 - ν[1]^2 - ν[2]'^2 + 0im))   # transfer function
    ũ = select_region_view(u, new_size=N)
    û = fftshift(ifft(fft(ifftshift(ũ)).*H))

    return select_region_view(û, new_size=size(u))
end

"""
    ASM!(u, λ, Δx, Δy, z; expand=true)

Same as ASM, but operates in-place on u, which must be an array of complex floating-point numbers.
"""
function ASM!(u, λ, Δx, Δy, z; expand=true)
    u[:,:] = ASM(u, λ, Δx, Δy, z; expand)
end
