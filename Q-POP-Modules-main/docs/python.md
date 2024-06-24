# Python Program for IMT

## Environment setup
The setup of FEniCS with Python is easier than for C++, you can follow the official installation guide [here](https://fenicsproject.org/download/archive/).

## Run the code 
```sh
python xxx.py
mpirun -n 4 python xxx.py  # running the code on 4 cores in parallel
```

## Brief explanation of the code
### File VO2_GSBlockATP.py 
-   For the block gauss seidel preconditioner code, we group the variables in the Allen-Cahn equation, and the variables in the PNP equations. We then give the following order: the Allen-Cahn equation, temperature equation and the PNP equation. The Gauss Seidel block preconditioner uses the lower triangular blocks of the three system of equations. 
-   Non-linear Solver: Newton solver 
-   Linear solver: Preconditioner GMRES, the block preconditioner is inverted using a direct solver at this stage. 

- The code allows the following automation:
    1. Choose adaptive time stepping by setting the varibale:  adaptive_time = True 
    2. User can input time schedule using the list variable: time_schedule_custom
    3. User can load time schedule from file: load_schedule, time_schedule_file
            - the time_schedule_file is read by the function read_in_time_schedule(filename)
    4. Initialize solution with read-in solutions: load_initialization

### File VO2_Nit_EfficAdap.py 
	1. An efficient adaptive time-stepping scheme is used

## Results Visualization
The program generates solutions in pvd format that can be read and plot by [ParaView](https://www.paraview.org "ParaView website"). The solution files are:
Name       | Explanation
---------- | -----------
`eta.pvd`  | Time-dependent structural order parameter field
`mu.pvd`  | Time-dependent electronic order parameter field
`phi.pvd`  | Time-dependent electric potential
`T.pvd`    | Time-dependent temperature field
`n.pvd`    | Time-dependent electron density field
`p.pvd`    | Time-dependent hole density field
`gamma_e.pvd` | Time-dependent transformed variable that corresponds to energy level of electrons
`gamma_h.pvd` | Time-dependent transformed variable that corresponds to energy level of holes 

Each solution file listed above links to data at different moments. The data at a moment in turn links to data parallelly computed on distributed processors. The user can just use ParaView to read the `.pvd` files and then follow the guide of [ParaView](https://docs.paraview.org/en/latest/UsersGuide/index.html "ParaView user's guide") to plot spatiotemporal fields.
