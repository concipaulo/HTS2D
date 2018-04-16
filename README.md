# Heat Transfer Solver 2D

---
Heat Transfer Solver 2D is a code to solve linear systems that comes from
discretization of physical domains. It currently solves by using GMRES
precondiotioned with the ILU family of preconditioners (ILU0, ILUT and MILU0).
This code is being used for my undergraduate thesis, so is constantly subject
to changes. So sometimes the code in the latest commit could not be working. I
recommend use the code from
[__releases__](https://github.com/concipaulo/HTS2D/releases) tab, there it will
be working for the recent case of study that I'm working on (and will be
working correctly).
---

## About the code

This code has been developed for some time and has the main objective to solve
systems of linear equations, with high rate of convergence, fast and modular.
It uses GMRES (_Generalized Minimal Residual_) that is a iterative method for
solve linear equations based on _Krylov_ subspace. 

It's wrote in Fortran 95/03, so in theory must be easy to add another
subroutines wrote in previous Fortran convections like (f90, f77). If you
intent to request a merge, keep this logic.
---

## Using the code

If you use the release version, it should be pretty straight forward, the
__compile__ file has and example of how to compile the code and run it. If you
make it executable (using chmod +x, for linux users) it can be run by simple
type `./compile` in the shell. 

### Changing the code

__Remember__ this is still a source code so it needs a little understanding of
_how_ to do it and _where_. How is pretty simple if you are familiar with
Fortran or any programming language, just follow the syntax and all be good. In
relation to where, here are some places you should take a look.

> Subroutines
  + bound_cond: this subroutine is responsible to set the boundary conditions
    of your problem. Most of the time, you will have a value in the boundary
that you want to set, like a temperature in some wall. These subroutine uses
*temp_dist* to set those values base in what coordinate you are. 
  +temp_dist: A little subroutine that has the function to set a temperature
value in some coordinate. Small and simple.
  + temp_exact: This subroutine is used to generate the exact solution for the
    problem if exists. It's quite simple, and it only require the mesh and the
function that represents the solution of the problem. 

> Modules
  + Constants: This module is responsible to set the constants of the linear
    system. Here we define if is going to be solve with second or fourth order.
Usually the variable represents where is going to be applied and what constant.

Here are some drawbacks from the recent update, there is no simple way to set
these constants, once they depend on mesh parameters, so __be aware__ that you
need to change this in every new mesh. A solution is coming.

  + parmesh: define the mesh name, precision, if is steady state or transient.
    It's well commented so just follow the instructions. 

--- 

## Contact

Feel free to open an issue or send me an email (pauloagconci@gmail.com).
