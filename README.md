# ReACT
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19668133.svg)](https://doi.org/10.5281/zenodo.19668133)

## Overview
ReACT is a reachability algorithm for linear time invarient systems. 


It introduces variable time steps and a new discretization technique found [Here.](ReACTDiscretize.jl)

## Usage
To use ReACT install a julia interpreter found [Here](https://julialang.org/downloads/) and install the relevant packages by running [The Setup File](setupLibraries.jl)

To start use the [Base Example](baseExample.jl) which produces a plot using ReACT. Here you can also adjust time step sizes and which model to plot.

To run the SLICOT benchmarks go to [Our Test Suite](testSuite.jl).
Additionally you can also run benchmarks for [LGG](LGG09Test.jl) and [BFFPSV](BFFPSV18Test.jl).

Importantly we refer to our approach (Algorithm 3) as ReACT found [Here.](ReACT.jl)

## Authors
This repository has been made by:
* Daniel Hilo Hansen
* Mikkel Bjørn
* Grace Melchiors


## License
MIT License

Copyright (c) 2026 Daniel Hilo Hansen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
