# HM-RK

 [![License][license-image]][license]
 
 [license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/lluis-mv/HM-RK/blob/main/LICENSE

Matlab codes to simulate the coupled hydro-mechanical problem at pseudo-static conditions (**u**-p<sub>w</sub> formulation) using the Finite Element method and node-based smoothed finite element method.

In addition to a fully implicit time integration approach, explicit Runge-Kutta time-marching schemes are also considered.


Several types of elements are considered: T3T3 employs linear interpolants for both displacement and water pressure, T6T3 uses a quadratic shape functions for displacements whereas water pressure is discretized with linear functions and T6T6 describes both field with quadratic interpolants. Q4Q4 employs bilinear shape functions for both nodal variables whereas Q8Q4 uses quadratic shape functions for displacements and bilinar interpolants for water pressure. 

A three field element has displacement, volumetric strain and water pressure as nodal variables, which are discretized with linear shape functions (T3T3T3) in triangles and bilinear interpolants (Q4Q4Q4) in quadrilater elements. These elements have two stabilization terms.

Additionally, node-based smoothed finite elements are considered. This implementation is only valid for linear triangular elements and bilinear quadrilateral elements.


Constitutive models: linear elasticity, non-linear elasticity (pressure-dependent hypo-elastic law) and Subloading Modified Cam Clay [(Hashiguchi, 2017)](https://doi.org/10.1007/978-3-319-48821-9)

### Matlab versions

The code has been tested in Matlab R2019a, R2020a and R2022b.

### Submodules

Some of the constitutive models (elastoplastic) are available in another [repository](https://github.com/lluis-mv/ExplicitStressIntegration)

## Sign convention
Computational geomechanics sign conventions is used in this implementation: in the water pressure field compression pressures are positive values whereas negative values of the effective Cauchy stress tensor denote compression.


## License and citation

These codes are distributed under BSD 3-Clause License. 


If you find these codes useful for your work, please cite as follows:
- Monforte, L., Carbonell, J.M., Arroyo, M. and Gens, A. (2022) An unconditionally stable explicit stabilized finite element for the coupled hydromechanical formulation in soil mechanics in pseudo-stationary conditions. International Journal for Numerical Methods in Engineering [https://doi.org/10.1002/nme.7064](https://doi.org/10.1002/nme.7064)

- Monforte, L., Collico, S., Carbonell, J.M., Arroyo, M. and Gens, A. (2023) Exploring the numerical performance of node-based smoothed finite elements in coupled hydro-mechanical problems. Computers and Geotechnics [https://doi.org/10.1016/j.compgeo.2023.105547](https://doi.org/10.1016/j.compgeo.2023.105547)
