## Introduction

This codebase contains implementation of various Neo-hookean strain energies discussed in the [Stable Neo-Hookean Flesh Simulation](https://graphics.pixar.com/library/StableElasticity/paper.pdf) paper. We also provide modifiers for changing various settings of the simulation such as poisson's ratio, young modulus, force multiplier and mesh (bunny or armadillo).  

[![Stable Neo-Hookean Flesh Simulation](https://github.com/yashkant/stable-neo-hookean-simulation/blob/master/neo-hookean-flesh.gif)](https://m.youtube.com/watch?v=aUB7hpQMkM0&amp;feature=youtu.be)

### Prerequisite installation

On all the platforms, we will assume you have installed cmake and a modern c++
compiler on Mac OS X[¹](#¹macusers), Linux[²](#²linuxusers), or
Windows[³](#³windowsusers).

We also assume that you have cloned this repository using the `--recursive`
flag (if not then issue `git submodule update --init --recursive`). 

### Layout


The `CMakeLists.txt` file setups up the cmake build routine for this
repository.

The `main.cpp` file will include the headers in the `include/` directory and
link to the functions compiled in the `src/` directory. This file contains the
`main` function that is executed when the program is run from the command line.

The `include/` directory contains one file for each function that you will
implement as part of the assignment.

The `src/` directory contains _ implementations_ of the functions
specified in the `include/` directory.

The `data/` directory contains _sample_ input data for your program.

## Compilation for Debugging

This and all following assignments will follow a typical cmake/make build
routine. Starting in this directory, issue:

    mkdir build
    cd build
    cmake ..

If you are using Mac or Linux, then issue:

    make

## Compilation for Testing

Compiling the code in the above manner will yield working, but very slow executables. To run the code at full speed, you should compile it in release mode. Starting in the **build directory**, do the following:

    cmake .. -DCMAKE_BUILD_TYPE=Release
    
Followed by:

    make 
  
Your code should now run significantly (sometimes as much as ten times) faster. 

If you are using Windows, then running `cmake ..` should have created a Visual Studio solution file
called `a3-finite-elements-3d.sln` that you can open and build from there. Building the project will generate an .exe file.


## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./neo-hookean-flesh-simulation
This will start a simulation with armadillo mesh using Smith et al. energy model.
You can use the options below to change energy model, Poisson's ratio, force multiplier and Younge Modulus:
```
add --poisson_rate 0.4 or -mu 0.4 to change Poisson's rate
add --younge_modulus 6000 or -ym 6000 to change Younge Modulus
add --force_multiplier 10 or -fm 10 to change force multiplier
add --energy_model smith_14 or -em smith_14 to change energy model. You can use 'smith_14', 'smith_13', 'ogden', 'wang' or 'bower'.
add --mesh bunny or -m bunny to change mesh to bunny.
add --sim_volume true or -sm true to run volume preservation simulation test.

```

