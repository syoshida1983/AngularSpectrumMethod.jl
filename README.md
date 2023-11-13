# AngularSpectrumMethod

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://syoshida1983.github.io/AngularSpectrumMethod.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://syoshida1983.github.io/AngularSpectrumMethod.jl/dev/)
[![Build Status](https://github.com/syoshida1983/AngularSpectrumMethod.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/syoshida1983/AngularSpectrumMethod.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This package provides the functions for diffraction calculations based on the angular spectrum method (ASM).
In addition to the standard diffraction calculation with ASM, various diffraction calculations with band-limited, scaled, shifted, and tilted ASM are available.

## ASM

The function `ASM` returns the diffracted field by the angular spectrum method (ASM).
Evanescent waves are not eliminated but attenuated as $\exp(-2{\pi}wz)$.
Without attenuation, the total energy $\iint|u|\mathrm{d}x\mathrm{d}y$ is conserved.

```julia
v = ASM(u, λ, Δx, Δy, z; expand=true)
```

### Arguments
- `u`: input field.
- `λ`: wavelength.
- `Δx`: sampling interval in the x-axis.
- `Δy`: sampling interval in the y-axis.
- `z`: diffraction distance.
- `expand=true`: if true (default), perform 4× expansion and zero padding for aliasing suppression.

> [!NOTE]
> The x-axis is the horizontal direction, and the y-axis is the vertical.

<p align="center">
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/rect.jpg" width="250px">
    &emsp;&emsp;
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ASM50mm.jpg" width="250px">
    <br>
    &emsp;&emsp;&emsp;
    input field
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    diffracted field by ASM
</p>

## Band-limited ASM

The function `BandLimitedASM` returns a diffracted field by the band-limited ASM (see Ref. 1).

```julia
v = BandLimitedASM(u, λ, Δx, Δy, z; expand=true)
```

> 1. [Kyoji Matsushima and Tomoyoshi Shimobaba, "Band-Limited Angular Spectrum Method for Numerical Simulation of Free-Space Propagation in Far and Near Fields," Opt. Express **17**, 19662-19673 (2009) ](https://doi.org/10.1364/OE.17.019662)

<p align="center">
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ASM5m.jpg" width="250px">
    &emsp;&emsp;
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/BandLimitedASM.jpg" width="250px">
    <br>
    ASM with a long propagation distance
    &emsp;&emsp;&emsp;&emsp;&emsp;
    Band-limited ASM
    &emsp;&emsp;&emsp;&emsp;&emsp;
</p>

## Scaled ASM

The function `ScaledASM` returns a scaled diffraction field according to the scale factor $R$ by the scaled ASM (see Ref. 2).

```julia
v = ScaledASM(u, λ, Δx, Δy, z, R; expand=true)
```

> 2. [Tomoyoshi Shimobaba, Kyoji Matsushima, Takashi Kakue, Nobuyuki Masuda, and Tomoyoshi Ito, "Scaled angular spectrum method," Opt. Lett. **37**, 4128-4130 (2012)](https://doi.org/10.1364/OL.37.004128)

<p align="center">
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ScaledASMx2.jpg" width="250px">
    &emsp;&emsp;
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ScaledASMx0.5.jpg" width="250px">
    <br>
    Scaled ASM with R = 2.0
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    R = 0.5
    &emsp;&emsp;&emsp;
</p>

## Shifted ASM

The function `ShiftedASM` returns a shifted diffraction field with the shift distance $x_{0}$ and $y_{0}$ by the shifted ASM (see Ref. 3).

```julia
v = ShiftedASM(u, λ, Δx, Δy, z, x₀, y₀; expand=true)
```

> 3. [Kyoji Matsushima, "Shifted angular spectrum method for off-axis numerical propagation," Opt. Express **18**, 18453-18463 (2010)](https://doi.org/10.1364/OE.18.018453)

<p align="center">
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ASM50mm.jpg" width="250px">
    &emsp;&emsp;
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/ShiftedASM.jpg" width="250px">
    <br>
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    ASM
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    Shifted ASM with x direction shift
</p>

## Tilted ASM

The function `TiltedASM` returns a tilted diffraction field for a rotation matrix $T$ by the tilted ASM (see Ref. 4, 5).
If `weight=true`, a diagonal weighting matrix is used as the Jacobian determinant (default `false`).
In this case, the energy conservation improves, but the computational cost is high (see Ref. 6).

```julia
v = TiltedASM(u, λ, Δx, Δy, T; expand=true, weight=false)
```

> 4. [Kyoji Matsushima, Hagen Schimmel, and Frank Wyrowski, "Fast calculation method for optical diffraction on tilted planes by use of the angular spectrum of plane waves," J. Opt. Soc. Am. A **20**, 1755-1762 (2003)](https://doi.org/10.1364/JOSAA.20.001755)
> 5. [Kyoji Matsushima, "Formulation of the rotational transformation of wave fields and their application to digital holography," Appl. Opt. **47**, D110-D116 (2008)](https://doi.org/10.1364/AO.47.00D110)
> 6. [James G. Pipe and Padmanabhan Menon, "Sampling density compensation in MRI: Rationale and an iterative numerical solution," Magn. Reson. Med. **41**, 179-186 (1999)](https://doi.org/10.1002/(sici)1522-2594(199901)41:1%3C179::aid-mrm25%3E3.0.co;2-v)

> [!NOTE]
> [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) is helpful in generating rotation matrices.

<p align="center">
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/rect.jpg" width="250px">
    &emsp;&emsp;
    <img src="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/images/TiltedASM.jpg" width="250px">
    <br>
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    input field
    &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
    Tilted ASM with rotation around the x-axis
</p>
