export TiltedASM
export TiltedASM!

"""
    TiltedASM(u, λ, Δx, Δy, T; expand=true, weight=false)

return tilted diffraction field for a rotation matrix ``T`` by the tilted ASM (see Ref. 1, 2).
If `weight=true`, a diagonal weighting matrix is used as the Jacobian determinant (default `false`).
In this case, the energy conservation improves, but the computational cost is high (see Ref. 3).

!!! note
    [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) is helpful in generating rotation matrices.

> 1. [Kyoji Matsushima, Hagen Schimmel, and Frank Wyrowski, "Fast calculation method for optical diffraction on tilted planes by use of the angular spectrum of plane waves," J. Opt. Soc. Am. A **20**, 1755-1762 (2003)](https://doi.org/10.1364/JOSAA.20.001755)
> 2. [Kyoji Matsushima, "Formulation of the rotational transformation of wave fields and their application to digital holography," Appl. Opt. **47**, D110-D116 (2008)](https://doi.org/10.1364/AO.47.00D110)
> 3. [James G. Pipe and Padmanabhan Menon, "Sampling density compensation in MRI: Rationale and an iterative numerical solution," Magn. Reson. Med. **41**, 179-186 (1999)](https://doi.org/10.1002/(sici)1522-2594(199901)41:1%3C179::aid-mrm25%3E3.0.co;2-v)
"""
function TiltedASM(u, λ, Δx, Δy, T; expand=true, weight=false)
    N = ifelse(expand, size(u').*2, size(u'))   # Note that the array is transposed for consistency with the rotation matrix
    Δ = [Δx, Δy]
    ν = fftfreq.(N, inv.(Δ))            # spatial frequencies νx, νy (DC corner) in the source plane
    νz² = @. 1/λ^2 - ν[1]^2 - ν[2]'^2   # squared spatial frequency νz²
    ν̂₀ = T*[0, 0, 1/λ]                  # carrier frequency in the reference plane
    ũ = select_region(transpose(u), new_size=N)
    û = fft(ifftshift(ũ))
    i = findall(==(true), νz²[:] .> 0)  # valid indices
    f̂ = û[i]                            # spectrum data
    ν̃ = @. [ν[1][(i - 1)%N[1] + 1]'; ν[2][(i - 1)÷N[1] + 1]'; √(νz²[i])']   # spatial frequencies in the source plane
    ν̂ = T*ν̃ .- ν̂₀               # spatial frequencies in the reference plane
    k̂ = ν̂[1:2,:].*Δ             # frequency node
    k̂ = @. k̂ - floor(k̂ + 1/2)   # periodic boundary [-1/2, 1/2)
    p = plan_nfft(k̂, N)::NFFTPlan{Float64, 2, 1}

    if weight
        û = adjoint(p)*(f̂.*sqrt.(sdc(p, iters=10)./length(f̂)))
    else
        # Jacobian determinant
        ν̂[3,:] .= @. √(1/λ^2 - ν̂[1,:]^2 - ν̂[2,:]^2)
        j = [T[4]*T[8] - T[7]*T[5],
            T[7]*T[2] - T[1]*T[8],
            T[1]*T[5] - T[4]*T[2]]
        J = .√([ν̂[1,:]./ν̂[3,:] ν̂[2,:]./ν̂[3,:] ones(eltype(ν̂), length(f̂), 1)]*j)
        û = adjoint(p)*(f̂.*J./length(f̂))
    end

    # superposition of carrier wave (positive directions of the coordinates are right and up)
    f = select_region(transpose(û), new_size=size(u))
    r = fftshift.(fftfreq.(size(u), size(u).*[Δy, -Δx]))

    return @. f*exp(2π*im*(ν̂₀[2]*r[1] + ν̂₀[1]*r[2]'))
end

"""
    TiltedASM!(u, λ, Δx, Δy, T; expand=true, weight=false)

Same as TiltedASM, but operates in-place on u, which must be an array of complex floating-point numbers.
"""
function TiltedASM!(u, λ, Δx, Δy, T; expand=true, weight=false)
    u[:,:] = TiltedASM(u, λ, Δx, Δy, T; expand, weight)
end
