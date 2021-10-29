# HM-RK
Implementation of explicit Runge-Kutta time-marching schemes for the coupled hydro-mechanical problem at pseudostationary conditions (\textbf{u}-p_w formulation) 

Three types of triangular elements are considered: T3T3 employs linear interpolants for both displacement and water pressure, T6T3 uses a quadratic shape functions for displacements whereas water pressure is discretized with linear functions and T6T6 describes both field with quadratic interpolants.

Constitutive models: linear elasticity, non-linear elasticity (pressure-dependent hypo-elastic law) and Subloading Modified Cam Clay [(Hashiguchi, 2017)](https://doi.org/10.1007/978-3-319-48821-9)

### Matlab versions

The code has been tested in Matlab R2019a and R2020a.

### Submodules

Some of the constitutive models (elastoplastic) are developed in another [repository](https://github.com/lluis-mv/ExplicitStressIntegration)

## Sign convention
Geomechanics sign conventions is used in this implementation (compressions are positive).

## License and citation

This code is distributed under BSD 3-Clause License. 


If you find these codes useful for your work, please cite as follows:
- Monforte, L., Carbonell, J.M., Arroyo, M. and Gens, A (2022) An unconditionally stable explicit stabilized finite element for the coupled hydromechanical formulation in soil mechanics in pseudo-stationary conditions [https://doi.org/10.1680/jgeot.19.p.354](https://doi.org/10.1680/jgeot.19.p.354)
