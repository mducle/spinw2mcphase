# SpinW2McPhase

A collection of matlab scripts which convert a [SpinW](https://spinw.org) object
to a set of input files for [McPhase](http://mcphase.de) and runs the calculation,
returning a SpinW-compatible `spec` structure for plotting.

SpinW requires Matlab, and can be installed as a Matlab "Add-on" 
([instructions](https://spinw.org/Ispra2019/#/install2)) or downloaded as a 
[zip](https://github.com/SpinW/spinw/releases/tag/v3.1) and manually installed.

To obtain a binary distribution of McPhase, please email <service@mcphase.de>.
Alternatively, you can obtain the source code
[here](https://github.com/mducle/mcphase) and compile it yourself. 
This requires the GNU C, C++ and Fortran compilers and GNU Make 
(e.g. use `msys` on Windows).
Thereafter you should be able to compile McPhase with `make` or `make fast=1`.

## Examples

The [examples](examples) folder contains examples m-files.
For them to work, you need to add the [spinw2mcphase](spinw2mcphase) folder
to your Matlab path (e.g. `addpath('/path/to/spinw2mcphase')`).
This contains two files: `mcphase_sqw.m` and `sw_plotspec.m`.
The `sw_plotspec.m` file is required because the original SpinW version
assumes that there are 2 * _N_ modes, where _N_ is the number of magnetic ions,
but for RPA calculcations the number of modes can exceed this.

Some examples require [Horace](https://github.com/pace-neutrons/Horace/) to be installed.
