## Finite Element Method for Stokes problem

This repo contains `MATLAB` implementation for Weak Galerkin FEM solver that can be used to solve the steady-state Stokes problem in 2D
```math
\begin{align*}
    -\Delta \mathbf{u} + \nabla p &= \mathbf{f},\quad \mathrm{in\ } \Omega\\
    \nabla \cdot \mathbf{u} &= 0,\quad \mathrm{in\ } \Omega\\
    \mathbf{u} &= \mathbf{g},\quad \mathrm{in\ } \partial\Omega\,
\end{align*}
```
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, $\mathbf{g}$ defines the Dirichlet type boundary condition (usually constant) and $\mathbf{f}$ is the load function. The solver was implemented as a project for course *MS-E1653 Finite Element Method* in Aalto University, but it is still nuder work.

The mesh generation part of the solver is based on the 2D FEM solver implementation that was given on the course.[1] The unmodified version of that solver can be found from the [MyCourses page](https://mycourses.aalto.fi/course/view.php?id=36259&section=4). The Weak Galerkin method with piecewise constant basis functions and their assembly was adopted from paper *"A SIMPLE FINITE ELEMENT METHOD FOR THE STOKES EQUATIONS"* written by Lin Mu and Xiu Ye in 2017.[2]

### References

[1] Hannukainen A., "*Finite element assembly for 2D meshes (util.zip)*", MS-E1653 Finite Eleme Method course, Aalto University, 2023

[2] Mu, L., Ye, X. *“A simple finite element method for the Stokes equations”*, Advances in Computational Mathematics, 43, 1305–1324, 2017, [https://doi.org/10.1007/s10444-017-9526-z](https://doi.org/10.1007/s10444-017-9526-z).