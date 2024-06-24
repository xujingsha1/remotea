# Visualization

## Visualization of solutions

Q-POP-IMT generates solutions in the .pvd format, which can be read and visualized by [ParaView](https://www.paraview.org "ParaView website"). The solution files are:

Name       | Explanation
---------- | -----------
`eta.pvd`  | Time-dependent structural order parameter field
`psi.pvd`  | Time-dependent electronic order parameter field
`phi.pvd`  | Time-dependent electric potential
`T.pvd`    | Time-dependent temperature field
`n.pvd`    | Time-dependent electron density field
`p.pvd`    | Time-dependent hole density field
`niov.pvd` | Time-dependent ionized oxygen vacancy concentration field
`nnov.pvd` | Time-dependent neutral oxygen vacancy concentration field

Each solution file listed above links to data at different moments. The data at a moment in turn links to data parallelly computed on distributed processors. The user can use ParaView to read the `.pvd` files and then follow the guide of [ParaView](https://docs.paraview.org/en/latest/UsersGuide/index.html "ParaView user's guide") to plot spatiotemporal fields.