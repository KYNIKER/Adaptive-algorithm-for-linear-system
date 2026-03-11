# ReACT
## Overview
ReACT is a reachability algorithm for linear time invarient systems. 

It introduces variable time steps and a new discretization technique found [Here.](ReACTDiscretize.jl)

## Usage
To use ReACT install a julia interpreter and install the relevant packages.

To start use the [base example](baseExample.jl) which produces a plot using ReACT.

To run the SLICOT benchmarks go to [testSuite](testSuite.jl), which can be adapted further.
Additionally you can also run benchmarks for [LGG](LGG09Test.jl) and [BFFPSV](BFFPSV18Test.jl).

Notably there exists ReACT and ReACTv2. Where ReACTv2 uses support functions more than ReACT. 

## Authors
This repository has been made by:
* Daniel Hilo Hansen
* Mikkel Bjørn
* Grace Melchiors

## License
MIT License

Copyright (c) 2025 Daniel Hilo Hansen

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
