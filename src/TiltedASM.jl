export TiltedASM

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
    û = ifelse(expand, ifftshift(fft(fftshift(padzeros(u)))), ifftshift(fft(fftshift(u))))
    N₁, N₂ = size(û)
    v̂₀ = (T*[0, 0, 1/λ])[1:2]                   # carrier frequency û₀, v̂₀ in the reference plane
    k̂ = Matrix{Float64}(undef, 2, length(û))    # frequency nodes in the reference plane
    f̂ = Vector{ComplexF64}(undef, length(û))    # spectrum data
    J = Vector{Float64}(undef, length(û))       # Jacobian determinant
    a = [T[4]*T[8] - T[7]*T[5],                 # for Jacobian calculation
        T[7]*T[2] - T[1]*T[8],
        T[1]*T[5] - T[4]*T[2]]
    n = 0

    @fastmath @inbounds for j ∈ 1:N₂, i ∈ 1:N₁
        v = [(i - 1 - N₁/2)/(N₁*Δx), (j - 1 - N₂/2)/(N₂*Δy)]    # spatial frequency u, v in the source plane
        w² = 1/λ^2 - v⋅v                                        # squared spatial frequency w² in the source plane

        # ignore evanescent wave
        if w² < 0
            continue
        end

        append!(v, √(w²))
        v̂ = (T*v)[1:2] - v̂₀    # spatial frequency û, v̂ in the reference plane
        t̂ = v̂.*[Δx, Δy]        # frequency node

        # bounds checking of nodes
        if true ∈ (@. abs(t̂) > 1/2)
            continue
        end

        n += 1
        append!(v̂, √(1/λ^2 - v̂⋅v̂))
        f̂[n] = û[i, j]
        J[n] = √(a⋅[v̂[1]/v̂[3], v̂[2]/v̂[3], 1])   # Jacobian determinant
        k̂[:, n] = t̂
    end

    p = plan_nfft((@view k̂[:, 1:n]), size(û))::NFFTPlan{Float64, 2, 1}
    if weight
        f = adjoint(p)*((@view f̂[1:n]).*sqrt.(sdc(p, iters=10)./n))
    else
        f = adjoint(p)*((@view f̂[1:n]).*(@view J[1:n])./n)
    end

    # superposition of carrier wave
    @inbounds @fastmath for j ∈ 1:N₂, i ∈ 1:N₁
        r = [(i - 1 - N₁/2)*Δx, (j - 1 - N₂/2)*Δy]
        f[i, j] *= exp(2π*im*v̂₀⋅r)
    end

    return ifelse(expand, crop(f), f)
end
