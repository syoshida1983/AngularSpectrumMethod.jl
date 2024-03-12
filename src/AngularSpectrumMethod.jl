module AngularSpectrumMethod

using FFTW
using LinearAlgebra
using NDTools
using NFFT
using NFFTTools

# Rectangular function
function rect(x)
    return ifelse(abs(x) < 1/2, one(x), zero(x))
end

include("ASM.jl")
include("BandLimitedASM.jl")
include("ScalableASM.jl")
include("ScaledASM.jl")
include("ShiftedASM.jl")
include("TiltedASM.jl")

end
