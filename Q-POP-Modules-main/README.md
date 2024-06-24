![GitHub release version](https://img.shields.io/github/v/release/DOE-COMMS/Q-POP-Modules?color=%2350C878&include_prereleases)
![License](https://img.shields.io/github/license/DOE-COMMS/Q-POP-Modules)
![GitHub Size](https://img.shields.io/github/repo-size/DOE-COMMS/Q-POP-Modules)
![HitCount](https://hits.dwyl.com/DOE-COMMS/Q-POP-Modules.svg?style=flat-square&show=unique)
![HitCount](https://img.shields.io/endpoint?url=https%3A%2F%2Fhits.dwyl.com%2FDOE-COMMS%2FQ-POP-Modules.json&label=total%20hits&color=pink)

# Q-POP Modules
The core modules of Q-POP (**Q**uantum **P**hase-field **O**pen-source **P**ackage) aim to deliver flexible, scalable, and extendable phase-field solvers that will enable researchers to study quantum materials and phase transitions with ease. The software will be made open-source in three stages, with new physics capabilities added with release. The module for insulator-metal transitions (Q-POP-IMT) was the first to be released in April 2024. The superconductor and the dynamical phase-field modules will follow soon.

## Acknowledgements
The development of this open-source software is supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award No. DE-SC0020145 as part of the Computational Materials Sciences Program.

## Folder Structure
```sh
├── stash
├── external
│   └── pugixml
│       ├── pugiconfig.hpp
│       ├── pugixml.cpp
│       └── pugixml.hpp
└── imt
    ├── VO2_GSBlockATP.py
    ├── qpop-imt.py
    └── examples
        └── input.xml
```

## Q-POP-IMT
This module is an open-source phase-field package for simulating mesoscale, nonequilibrium, inhomogeneous processes of insulator-metal phase transitions and their associated electrical and structural responses.

### Setup Development Environment
This module uses FEniCS C++ library and its Python interface for defining and solving finite-element partial differential equations, and is capable of exascale simulations. The FEniCS C++ library and its Python interface of version 2019.1.0.post0 must be installed. Note that this version of FEniCS is compatible only with openMPI v3.1 or older versions, due to the dependence of FEniCS on mpi4py, which itself does not have a stable OpenMPI 4+ implementation currently. 

### Running the Program
The program supports parallel computing. One can run the Python interface of the module without any compilation. To run it on, say, 8 processors, run the below command in your desired directory:
```
mpirun -np 8 python directory-to-python-script/qpop-imt.py
```
The program requires an input file for specifying parameters. The output files will be generated in the current directory.

### Input File
Users can provide various parameters in the `input.xml` file to control the simulation. There are several sections in `input.xml`:
Section            | Explanation
----------         | ------------
`internal`         | Internal parameters of the material, e.g., Landau potential coefficients
`external`         | External parameters, e.g., system sizes and ambient temperature
`time`             | Parameters related to the time, e.g., simulation end time
`initialization`   | Parameters for initialization, e.g., initial temperature
`solverparameters` | Parameters for finite-element solver, e.g., linear solver choice

Currently, the program considers one structural order parameter $\eta$ and one electronic order parameter $\psi$, and the Landau potential has the form 
```math
\begin{aligned}
f_L=& \frac{a_1(T-T_1)}{2T_c}\eta^2 + \frac{a_2}{4}\eta^4 + \frac{a_3}{6}\eta^6  \\
    & + \frac{b_1(T-T_2)}{2T_c}\psi^2 + \frac{b_2}{4}\psi^4 + \frac{b_3}{6}\psi^6  \\
    & + c_1\eta\psi - \frac{c_2}{2}\eta^2\psi^2 + \frac{c_3}{2}\eta^3\psi ,
\end{aligned}
```
where $T$ is temperature. The energy gap form is approximated by the symmetry-allowed lowest order term on the electronic order parameter, $E_g = E_{g0}\psi^2$, where $E_{g0}$ is the gap coefficient. The detail of the phase-field model of insulator-metal transitions can be found in [Y. Shi and L.-Q. Chen, 2019](https://doi.org/10.1103/PhysRevApplied.11.014059 "Current-Driven Insulator-To-Metal Transition in Strongly Correlated VO2"). The Landau coefficients, gap coefficient, etc., can be input in the `internal` section like
```
<internal>
 <a1 unit='kB*Tc/u.c.'>1.0</a1>
 <T1 unit='K'>200.0</T1>
 <gapcoeff unit='eV'>0.5</gapcoeff>
 ...
</internal>
```
The units are fixed in the program and here they are only to remind the users (u.c. means unit cell volume). Other internal parameters are:
Name | Explanation
---- | -----------
`ucellvol` | Unit cell volume [nm<sup>-3</sup>]
`SPTtimeconst` | Time constant of the structural phase transition [ns]
`IMTtimeconst` | Time constant of the insulator-metal transition [ns]
`excitdensity` | Excited electron density inducing the insulator-metal transition [u.c.<sup>-1</sup>]; attribute of `IMTtimeconst`
`SOPgradcoeff` | Gradient energy coefficient for the structural order parameter [eV/nm]
`EOPgradcoeff` | Gradient energy coefficient for the electronic order parameter [eV/nm]
`conductbandeffDOS` | Effective density of states of the conduction band [u.c.<sup>-1</sup>]
`valencebandeffDOS` | Effective density of states of the valence band [u.c.<sup>-1</sup>]
`elecmobilityx` | Electron mobility in the $\hat{x}$ direction [cm<sup>2</sup>/(Vs)]
`elecmobilityy` | Electron mobility in the $\hat{y}$ direction [cm<sup>2</sup>/(Vs)]
`holemobilityx` | hole mobility in the $\hat{x}$ direction [cm<sup>2</sup>/(Vs)]
`holemobilityy` | hole mobility in the $\hat{y}$ direction [cm<sup>2</sup>/(Vs)]
`ehrecrate` | Electron-hole recombination rate constant [u.c./ns]
`relativepermitt` | Relative permittivity
`volheatcapacity` | Volumetric heat capacity at constant pressure [J/(m<sup>3</sup>K)]
`thermalconduct` | Thermal conductivity [W/(mK)]
`heatdiss` | Heat dissipation coefficient [W/(m<sup>2</sup>K)]

If not explicitly specified, the internal parameters will default to those of VO<sub>2</sub>, a prototypical strongly correlated material exhibiting an insulator-metal transition near room temperature.

The example below simulates a rectangular VO<sub>2</sub> device, supplied with a direct voltage through a series resistor.
```
<?xml version="1.0"?>
<input>
 <external>
  <temperature unit='K'>300.1</temperature>
  <voltage unit='V'>5.6</voltage>
  <resistor unit='Ohm'>2.5e5</resistor>
  <capacitor unit='nF'>0.0</capacitor>
  <heatdiss unit='W/m^2K'>3e6</heatdiss>
  <Lx unit='nm' mesh='100'>100.0</Lx>
  <Ly unit='nm' mesh='40'>40.0</Ly>
  <Lz unit='nm'>20.0</Lz>
 </external>
 <time>
  <endtime unit='ns'>2e3</endtime>
  <savemethod>auto</savemethod>
  <saveperiod>20</saveperiod>
 </time>
 <initialization>
  <temperature unit='K'>300.0</temperature>
  <SOP>1.119</SOP>
  <EOP>-1.293</EOP>
  <Tcvariance method='nucleus1'>
   <Tcshift unit='K'>-5.0</Tcshift>
   <radius unit='nm'>3.0</radius>
  </Tcvariance>
 </initialization>
 <solverparameters>
  <Newtonabsolutetolerance>1e-5</Newtonabsolutetolerance>
  <Newtonrelativetolerance>1e-3</Newtonrelativetolerance>
  <Newtonmaxiteration>15</Newtonmaxiteration>
  <timesteptolerance>1e-2</timesteptolerance>
  <directsolver>pastix</directsolver>
  <useGSblockpreconditioner>0</useGSblockpreconditioner>
  <loglevel>INFO</loglevel>
 </solverparameters>
</input>
```
Almost all the parameters are self-explanatory. For the `external` section: 
Name          | Explanation
------------- | -------------------
`temperature` | Ambient temperature
`voltage`     | Direct voltage applied
`resistor`    | Resistance of the series resistor
`capacitor`   | The capacitance representing the parasitic capacitance or the external capacitor parallelly connected with the series resistor
`heatdiss`    | Heat transfer coefficient from the device to the environment
`Lx`          | Width of the device
`Ly`          | Length of the device, along the electric field
`Lz`          | Thickness of the device
`mesh`        | Number of mesh vertices along a given dimension

For the `time` section:
Name         | Explanation
------------ | -----------
`endtime`    | Simulation end time; beginning time is zero
`savemethod` | Method of saving time-dependent solutions, used with `saveperiod` parameter. `auto` means to save every `saveperiod`-step solution; `fixed` means to save every `saveperiod`-nanosecond solution
`saveperiod` | See `savemethod`

For the `initialization` section:
Name          | Explanation
------------- | -----------
`temperature` | Initial temperature of the device
`SOP`         | Initial value for the structural order parameter
`EOP`         | Initial value for the electronic order parameter
`Tcvariance`  | `method`: how to set up a nucleus of the high-temperature phase. `nucleus1` means to set up a semicircle with a radius of `radius` and a transition temperature shift of `Tcshift`, located at the midpoint of the $y = 0$ edge

The `solverparameters` section defines the parameters for both the Newton's iteration solver (for nonlinear differential equations) and the linear solver (for solving linear equations in each iteration of the Newton's method):
Name                      | Explanation
------------------------- | -----------
`Newtonabsolutetolerance` | Absolute tolerance for the Newton iteration
`Newtonrelativetolerance` | Relative tolerance for the Newton iteration. Meeting either the absolute tolerance or the relative tolerance is considered converged
`Newtonmaxiteration`      | Limit of the iteration times for the Newton's method
`timesteptolerance`       | Relative tolerance for the adaptive time stepping error
`directsolver`            | Which direct solver to use for solving the linear problem
`useGSblockpreconditioner` | Whether to use our developed block Gauss-Seidel preconditioned GMRES iterative solver to solve the linear problem
`loglevel`                | Log level; see [FEniCS manual](https://fenics.readthedocs.io/projects/dolfin/en/2017.2.0/apis/api_log.html "FEniCS log level")

This setup produces the intrinsic voltage self-oscillation in VO<sub>2</sub> thin films, which was published in Physical Review Applied; see [Y. Shi and L.-Q. Chen, 2022](https://doi.org/10.1103/PhysRevApplied.17.014042 "Intrinsic voltage self-oscillation"). In our test runs, the simulation took 2 hours to complete using 16 processors of an AMD EPYC 7742 CPU.

### Visualization of solutions
The program generates solutions in pvd format that can be read and plotted by [ParaView](https://www.paraview.org "ParaView website"). The solution files are:
Name       | Explanation
---------- | -----------
`eta.pvd`  | Time-dependent structural order parameter field
`psi.pvd`  | Time-dependent electronic order parameter field
`phi.pvd`  | Time-dependent electric potential
`T.pvd`    | Time-dependent temperature field
`n.pvd`    | Time-dependent electron density field
`p.pvd`    | Time-dependent hole density field

Each solution file listed above links to data at different moments. The data at a moment in turn links to data parallelly computed on distributed processors. The user can just use ParaView to read the `.pvd` files and then follow the guide of [ParaView](https://docs.paraview.org/en/latest/UsersGuide/index.html "ParaView user's guide") to plot spatiotemporal fields.

The program also generates and updates on the fly a log.txt file that summarizes the temporal evolution of several physical quantities and diagnostics. The columns are summarized in the table below.
Column name    | Explanation
-------------- | -----------
`#Step`        | Number of time steps
`Time`         | Time in nanosecond
`Time step`    | Current time step size
`Tfail`        | Number of times the adaptive time stepping shrinks the time step
`Nfail`        | Number of times the Newton's iterative solver diverges
`Other fail`   | Number of times the linear solver diverges
`Av. EOP norm` | Spatially averaged electronic order parameter norm
`Av. T (K)`    | Spatially averaged temperature
`V (V)`        | Voltage drop across the VO<sub>2</sub> film
`R (Ohm)`      | Resistance of the VO<sub>2</sub> film

The time-dependent voltage drop across the VO<sub>2</sub> film (Column `V (V)` in `log.txt`) generated in the example is shown in the figure below.

<p align="center">
<img src="https://github.com/DOE-COMMS/Q-POP-Modules/files/12205719/sch_osc.pdf" alt="Temporal evolution of the voltage drop across the film." width="500">
</p>

You can see the self-oscillation of the voltage output. The `psi.pvd` generated in the example is shown as several snapshots below.

<p align="center">
<img src="https://github.com/DOE-COMMS/Q-POP-Modules/files/12197852/morph.pdf" alt="Spatiotemporal evolution of the electronic order parameter" width="500">
</p>

$\psi$ represents the electronic phases: $\psi=0$ means metal while $|\psi|\sim 1$ means insulator. You will see a metallic filament automatically growing and shrinking back and forth, generating the oscillating voltage output across the VO<sub>2</sub> film.
