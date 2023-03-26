# MonolithicFEMVLFS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oriolcg.github.io/MonolithicFEMVLFS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oriolcg.github.io/MonolithicFEMVLFS.jl/dev)
[![Build Status](https://github.com/oriolcg/MonolithicFEMVLFS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/oriolcg/MonolithicFEMVLFS.jl/actions/workflows/CI.yml?query=branch%3Amain)

A monolithic Finite Element formulation for the hydroelastic analysis of Very Large Floating Structures

This repository contains all the tests performed in the manuscript:
*A Monolithic Finite Element formulation for Very Large Floating Structures* by Oriol Colomés, Francesc Verdugo and Ido Akkerman. If you use this formulation, please cite:
```
@article{colomes2023monolithic,
  title={A monolithic finite element formulation for the hydroelastic analysis of very large floating structures},
  author={Colom{\'e}s, Oriol and Verdugo, Francesc and Akkerman, Ido},
  journal={International Journal for Numerical Methods in Engineering},
  volume={124},
  number={3},
  pages={714--751},
  year={2023},
  publisher={Wiley Online Library}
}
```
and
```
@software{Colomes_MonolithicFEMVLFS_2022,
    author = {Colomes, Oriol},
    doi = {10.4121/19601419},
    month = {4},
    title = {{MonolithicFEMVLFS.jl}},
    url = {https://github.com/oriolcg/MonolithicFEMVLFS.jl},
    year = {2022},
    version = {0.1.0}
}
```
### Abstract
In this work we present a novel monolithic Finite Element Method (FEM) for the hydroelastic analysis of Very Large Floating Structures (VLFS) with arbitrary shapes that is stable, energy conserving and overcomes the need of an iterative algorithm. The new formulation enables a fully monolithic solution of the linear free-surface flow, described by linear potential flow, coupled with floating thin structures, described by the Euler-Bernoulli beam or Poisson-Kirchhoff plate equations. 

The formulation presented in this work is general in the sense that solutions can be found in the frequency and time domains, it overcomes the need of using elements with $ C^1 $ continuity by employing a continuous/discontinuous Galerkin (C/DG) approach, and it is suitable for finite elements of arbitrary order.

We show that the proposed approach can accurately describe the hydroelastic phenomena of VLFS with a variety of tests, including structures with elastic joints, variable bathymetry and arbitrary strucutral shapes.

## Installation
`MonolithicFEMVLFS` is a package registered in the official [Julia package registry](https://github.com/JuliaRegistries/General).  Thus, the installation of this package is straight forward using the [Julia's package manager](https://julialang.github.io/Pkg.jl/v1/). Open the Julia REPL, type `]` to enter package mode, and install as follows
```julia
pkg> add MonolithicFEMVLFS
```

## Usage
To run all the test cases in the paper do:
```julia
using MonolithicFEMVLFS
run_tests("all")
```

To run only a specific test, for example the Khabakhpasheva test in frequency domain, do:
```julia
using MonolithicFEMVLFS
run_tests("5-2-1-Khabakhpasheva-freq-domain.jl")
```
Note that the numbers in front of the script indicate the section in the manuscript.

After execution, the data will be stored in the respective folder `data/<Section-number>-<test-name>`. If the flag to generate VTK files is active, the VTK output will be stored in `data/VTKOutput/<Section-number>-<test-name>`. The plots shown in the manuscript are stored in `plots/<Section-number>-<test-name>`.

This repository uses DrWatson package, the data will only be generated the first time the tests are executed. If the data is already stored, the scripts will only regenerate the figures.

The code snipped appearing in Figure 3 of the manuscript can be found in `src/lst_periodic_beam.jl`.
