# Environment Setup

:::{caution}
For a hassle-free installation, it is recommended to use Ubuntu 20.04.
:::

## Installing FEniCS

### GNU 9
It is necessary to use GNU Compiler Collection 9. On systems that default to other versions of the GNU compiler, manually switching to GNU-9 will be required.
:::{caution}
Always exercise caution while running commands using sudo. If you do not have root access, please contact your system administrator.
:::
```sh
sudo apt install gcc-9 g++-9 gfortran-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90 --slave /usr/bin/gcov gcov /usr/bin/gcov-9
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 90
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-9 90
sudo update-alternatives --config gcc
sudo update-alternatives --config g++
sudo update-alternatives --config gfortran
```

### OpenMPI-3
All testing was done with OpenMPI-3.1.6. Newer versions of OpenMPI may not work due to API changes that are not supported by FEniCS.
Since it is no longer under active development, it can only be obtained from source. To avoid conflicts with existing MPI installations on your system, install it in your home directory.
Please ensure that the GNU-9 compilers are the default while installing OpenMPI.
:::{Note}
mpi4py, if already installed, will need to uninstalled and reinstalled after the following steps.
:::
```sh
wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz
tar xvzf openmpi-3.1.6.tar.gz
cd openmpi-3.1.6
./configure --prefix=$HOME/openmpi3
make all
make install
MPI_DIR=$HOME/openmpi3
echo "export LD_LIBRARY_PATH=$MPI_DIR/lib:$LD_LIBRARY_PATH; export PATH=$MPI_DIR/bin:$PATH" >> ~/.bashrc
pip3 install -U mpi4py 
```
To switch back to your pre-existing MPI installations, remove/decrease the priority of the newly added values to the LD_LIBRARY_PATH and PATH environment variables.

### Boost 1.71
In newer versions of Boost, the following error may be encountered:
```sh
dolfin/io/VTKWriter.cpp:27:10: fatal error: boost/detail/endian.hpp: No such file or directory
   27 | #include <boost/detail/endian.hpp>
      |          ^~~~~~~~~~~~~~~~~~~~~~~~~
compilation terminated.
```
All testing was done using Boost-1.71.0, which is fully compatible with FEniCS 2019. On systems running Ubuntu 20.04, `sudo apt install libboost-dev` may be used to install Boost-1.71.0. On newer versions of Ubuntu, Boost will need to be built from source:
```sh
wget https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz
tar -xzf boost_1_71_0.tar.gz
cd boost_1_71_0/
./bootstrap.sh --prefix=$HOME/boost
./b2 install
```
This will install Boost-1.71.0 to the system's home directory. To install elsewhere, the prefix option in the above command must be set to a desired path.

### PETSc
The linear matrix solvers used by Q-POP-IMT are provided by PETSc. All testing was done with PETSc-3.15.1, though any version newer than 3.7 is compatible with FEniCS. To install, please refer to the PETSc (https://petsc.org/release/install/download/) and petsc4py (https://pypi.org/project/petsc4py/) documentation. A minimal configuration spec is 
```sh
--download-metis --download-parmetis --download-ptscotch --download-suitesparse --download-mumps --download-scalapack --download-hypre --with-mpi-dir=$MPI_DIR --with-petsc4py
```

### Dolfin
```sh
sudo apt install git libeigen3-dev python3-pip cmake 
pip3 install fenics-ffc --upgrade
echo 'export PATH=$PATH:$HOME/.local/bin' >> ~/.bashrc
git clone --branch=2019.1.0.post0 https://bitbucket.org/fenics-project/dolfin.git
cd dolfin
git describe --tags --abbrev=0
mkdir build
cd build
```
Then, on Ubuntu 20.04, use 
```sh
cmake -DCMAKE_INSTALL_PREFIX=$HOME/fenics/dolfin ..
```
Or, on systems with Ubuntu 22.04, use
```sh
cmake -DCMAKE_INSTALL_PREFIX=$HOME/fenics/dolfin -DBOOST_ROOT=$HOME/boost ..
```
To install FEniCS C++ and Python both,
```sh
make
make install
cd ../python
pip3 install .
echo "source $HOME/fenics/dolfin/share/dolfin/dolfin.conf" >> ~/.bashrc
```
:::{note}
The Python interface for FEniCS has various dependencies. Refer the FEniCS (https://fenics.readthedocs.io/en/latest/installation.html) page.
:::
A new terminal must be opened for changes to take effect. Alternatively, `source ~/.bashrc` can be used.

## Installing Sphinx 

:::{attention}
Optional: Build documentation locally
:::

```sh
pip install Sphinx==6.2.1
pip install myst-parser sphinxcontrib-bibtex sphinx_design sphinx_copybutton sphinxcontrib-mermaid sphinx_multiversion sphinx-book-theme
```

```sh
cd docs
make html
```