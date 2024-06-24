# Insulator-Metal Transition module

1. Enviroment configuration: 

2. how to run the code 
- python xxx.py
- mpirun -n 4 python xxx.py  # running the code on 4 cores in parallel

3. Explanation for the file VO2_GSBlockATP.py 
-    For the block gauss seidel preconditioner code, we group the variables in the Allen-Cahn equation, and the variables in the PNP equations. We thee give the following order: the Allen-Cahn equation, temperature equation and the PNP equation. The Gauss Seidel block preconditioner uses the lower triangular blocks of the three system of equations. 
-   Nonlinear Solver: Newton solver 
-   linear solver: Preconditioner GMRES, the block preconditioner is inverted using a direct solver at this stage. 

- The code allows the following automation:
    1. Choose adaptive time stepping by setting the varibale:  adaptive_time = True 
    2. User can input time schedule using the list variable: time_schedule_custom
    3. User can load time schedule from file: load_schedule, time_schedule_file
            - the time_schedule_file is read by the function read_in_time_schedule(filename)
    4. Initialize solution with read-in solutions: load_initialization

4. Explanation for the file VO2_Nit_EfficAdap.py 
	1. an efficient adaptive time stepping scheme is used

## Visualization of solutions
The program generates solutions in pvd format that can be read and plot by [ParaView](https://www.paraview.org "ParaView website"). The solution files are:
Name       | Explanation
---------- | -----------
`eta.pvd`  | Time-dependent structural order parameter field
`mu.pvd`  | Time-dependent electronic order parameter field
`phi.pvd`  | Time-dependent electric potential
`T.pvd`    | Time-dependent temperature field
`n.pvd`    | Time-dependent electron density field
`p.pvd`    | Time-dependent hole density field

Each solution file listed above links to data at different moments. The data at a moment in turn links to data parallelly computed on distributed processors. The user can just use ParaView to read the `.pvd` files and then follow the guide of [ParaView](https://docs.paraview.org/en/latest/UsersGuide/index.html "ParaView user's guide") to plot spatiotemporal fields.
