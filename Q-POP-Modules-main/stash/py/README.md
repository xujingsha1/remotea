## A brief description of the codes
 
boundaryUtilities.py:
> This module handles the set up of Dirichlet and Periodic Boundary Conditions. I have some other codes that I have written for robin and neumann conditions, but I am still working those into this module in a user-friendly manner. Once I get that done I will build those into this module. The goal of this module is to merely provide an interface for a user to say “hey I want a Dirichlet BC at this point of the mesh with this value” and then the module handles all of the necessary fenics code for creating that boundary condition. Same goes for periodic conditions.
 
meshUtilities.py:
> This module handles the setup of simple fenics meshes. At the moment I just have rectangular meshs, but eventually we will want to add in the possibility for a user to import their mesh (Dan I think you already do this, so perhaps that is something you can handle). Again, a pretty simple interface where the user just provides some basic information and the module handles the necessary meshing.
 
Solvers.py:
> This is my pride and joy and is still under a bit of construction as I move some of my ugly solvers into the format I am going for here. Basically (as you will see in some of the *test.py files) this module allows a very easy interface for the user/developer to call a specific solver and a solver object/class is returned. Then we can simply call a solve routine to solve the equation, which also handles updating the function fields (think AC, CH, or transient problems). I have a vectorized Allen-Cahn/TDGL solver built in a separate code along with elastodynamics and polarization dynamics solvers but I need to convert them into this style which is more modular and simpler. There are some other pieces that need to eventually get built into here (like allowing the user/developer to change solver parameters and making the default solver parameters to be the best we can find via testing)
 
tensorUtilities.py:
> This handles some common tensor and vector operations that we do. Think converting from voigt into tensor notation and vice versa as well as doing stress/strain calculations using both Lame parameters and the typical stiffness tensor. It can also handle tensor assignment which I am working on an example to demonstrate that capability.
 
terminalPrinting.py:
> As one might expect, this handles a lot of stuff that is printed to the screen. Pretty self-explanatory.
 
*test.py:
> This are some simple test codes I have put together to test the different solvers I have put together. I used examples from the fenics website as it would be a good demonstration of the efficiency and ease-of-use I tried to build into the code.
 
ac_ener.py:
> This really small piece of code is used with ac_2.test.py. Basically what I am currently working on doing, is creating a robust way so a user/developer/us can create different energy functionals to be solved without the need for hard-coding directly into QPOP. For example, consider $\frac{\partial \eta}{\partial t} = -t\frac{\delta F}{\delta \eta}$ We know that the driving force is $\frac{\delta F}{\delta \eta}$ which can be calculated quite easily most of the time. Thus, the idea I have is instead of hard coding a MIT-driving force, an SC-driving force, [insert process here]-driving force, to just treat this as we would most other things, as simply another input. Therefore, the code can be quickly repurposed for many different projects/research. Gradient effects can be written into these driving forces as well pretty easily using the .dx() operation with a FEniCS function object. Obviously still have some huddles here to include different order parameters or fields (like elastic, electric, or thermal), but shouldn’t be impossible to get included.

> This modular idea also allows the user/developer/etc. to specify their initial conditions as well, and thus we could also include an “energy” term if we wanted to track different energy contributions during a calculation. This has been done a little bit in one of my one-off codes which I am working to get included into this format.
 
## Future Steps/Current Work In Progress:
 
Input/output:
> I have a large chunk of this complete as well but need to get it converted into the current style. As we start to include different meshes or things like “ac_ener.py” we will need a way to fix that, but those are relatively small potatoes. I will get this out of the jupyter notebooks I was working in and into the style over the next few weeks. Output is relatively simple with FEniCS because everything is done in-house with vtk/xdmf formatting but I just want a qpop interface again for simplicity.
 
Timestepping:
> As you can see in ch.test.py, ac.test.py, and heat.test.py the time stepping is handled in a loop in the main code. I would like to eventually create a separate module that can perform this and just make it a simple call. Obviously once we get into more complex problems, where we need to solve for magnetic, electrostatic, thermal, chemical, and elastic equilibrium for example, we need to method for handling that, but I have seen some interesting code snippets on the web that do something similar, as well as I have thoughts on doing this from some different coding projects I have done. Essentially, we would have a list of solves and then just loop through the list doing the solves in the order provided. We would need to work through updating various components of different solvers but overall this should work well.
 
Adaptivity:
> If we can eventually develop an AMR (adaptive mesh refinement) scheme into the code that would be a boon to our performance and also provide a capability nothing else in the phase field market currently has (as far as I know). In addition, some adaptive time stepping schemes could also be helpful and that is something I can handle as well.
 
Implicit/Explicit/Crank-Nicholson Schemes:
> For the time being, I am using implicit time stepping for AC and Heat/Diffusion problems, and a CN scheme (where the user can change theta as they see fit) for CH problems. I would like to eventually make them all CN schemes but just been lazy and hasn’t included that yet. It wouldn’t be a wholy difficult process but would take a fair amount of testing.
 
Convergence, Speed, and Accuracy:
> We need to test the solvers across different preconditioners and linear solvers to determine the most efficient for each particular type of problem. Then we can make these default parameters for the solvers which will be good and easy. But still allow the user to change the parameters if they want to. I have done this already in some jupyter notebooks but need to do a more thorough test and then implement the parameters within the solvers themselves.
 
## More Examples:
I am working on finishing a REALLY basic MIT problem I got from Yin, which will serve as our proof of concept for this work. Dan whenever you can if you can share a BASIC SC problem I can work on trying to get included as well (speaking of which I know your complex TDGL solver is going to be different as well so I can work on implementing that as well). I am also putting together some examples from the NIST phase field benchmarks site that is run by a team out of Northwestern so we can provide even more examples to this library. Essentially we want a phase-field package that is focused on solving quantum and SC problems, but I think this library can also be a good alternative to MOOSE and others because those interfaces and installations are much different and don’t quite allow the rapid prototyping and research that an architecture like ours can provide.