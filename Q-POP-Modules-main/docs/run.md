# Build and Run

## Build
In the base directory 'Q-POP-Modules', for Ubuntu-20.04 run
```sh
source $HOME/fenics/dolfin/share/dolfin/dolfin.conf
cmake -B build -S .
cmake --build build
```
and if using Ubuntu-22.04, run
```sh
source $HOME/fenics/dolfin/share/dolfin/dolfin.conf
cmake -DBOOST_ROOT=$HOME/boost -B build -S .
cmake --build build
```

The above can also be run to similar effect by combining the first two commands. This builds Q-POP-IMT without adding Dolfin's paths to the environment. For example, on Ubuntu-20.04:
```sh
cmake -DDOLFIN_DIR=$HOME/fenics/dolfin/share/dolfin/cmake -B build -S .
cmake --build build
```

## Run

To run the simulation using Python, navigate to Q-POP-Modules/imt and then use
```sh
mpirun -np 2 python3 ./qpop-imt.py
```

To run the simulation using C++, there are two equivalent methods:
1. You can change directory (cd) to the `build/imt/imt-cpp` folder, create an input file (see 'input.md' for an example), and then run the `qpop-imt` executable directly from the build directory using
```sh
mpirun -np 2 qpop-imt
```
where the number following the option 'np' is the number of processors to use.

2. Or you can add the path `build/imt/imt-cpp` folder to the PATH environment variable using `PATH="<ADD PATH TO Q-POP-MODULES HERE>/build/imt/imt-cpp${PATH:+:${PATH}}"`, which will allow the use of `qpop-imt` as a command from any folder with an input file.
