module AngularSpectrumMethod

using FFTW
using LinearAlgebra
using NFFT
using NFFTTools

# Rectangular function
function rect(x)
    return ifelse(abs(x) < 1/2, one(x), zero(x))
end

# expand and pad zeros
function padzeros(x)
    y = zeros(eltype(x), 2 .*size(x))
    u, l = round.(Int, size(x)./2, RoundNearestTiesAway) .+ 1
    d, r = round.(Int, size(x)./2 .*3, RoundNearestTiesAway)
    y[u:d, l:r] = @view x[:, :]
    return y
end

# crop center
function crop(x)
    u, l = round.(Int, size(x)./4, RoundNearestTiesAway) .+ 1
    d, r = round.(Int, size(x)./4 .*3, RoundNearestTiesAway)
    return x[u:d, l:r]
end

include("ASM.jl")
include("BandLimitedASM.jl")
include("ScaledASM.jl")
include("ShiftedASM.jl")
include("TiltedASM.jl")

end
