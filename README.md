# Interface Reconstruction Library (IRL)

The Interface Reconstruction Library (IRL) is a library of computational algorithms geared towards enabling easy and fast implementations of geometric Volume of Fluid (VOF) schemes for simulating multiphase flows. The methods contained in IRL, however, are general and could be used to compute volume-volume intersections, interface reconstructions, or perform optimizations. IRL was also specifically designed to be agnostic of the underlying meshes, making it useable for structured, unstructured, and adaptively refined meshes.

Additional documentation and examples are currently being added. If you would like to become involved in its development or inquire about potentially using it for your application, please contact the developers.

## Developers
- [Robert Chiodi] (principal developer, maintainer)
- [Fabien Evrard](mailto:fabien.evrard@cornell.edu)

# License and scope
IRL is open-sourced under the Mozilla Public License 2 (MPL2). It has already been used in incompressible and compressible CFD codes on structured, unstructured, and AMReX based flow solvers, with some images of simulations performed using IRL shown below.

# Start guide
IRL is written in C++17 with a Fortran 2008 interface. Below are links to additional information on how to install and use IRL:

- [Compilation & Installation](docs/markdown/install_main_page.md)

- [C/Fortran Interface Usage & Description](docs/markdown/interface_main_page.md)

- [Brief Introduction](docs/reference_powerpoint.pdf)

- [Examples](docs/markdown/examples_main_page.md)

## How to cite us
If you use IRL for your scientific work, please consider citing the [paper](10.1016/j.jcp.2021.110787) introducing the half-edge structure and its implementation

    R. Chiodi and O. Desjardins, General, robust, and efficient polyhedron intersection in the interface reconstruction library, Journal of Computational Physics 449 (2022), 110787. 

If you use IRL for paraboloid-polytope intersection, please consider citing the [preprint](https://arxiv.org/abs/2210.07772) introducing the closed-form expressions of the moments and their implementation

    F. Evrard, R. Chiodi, A. Han, B. van Wachem and O. Desjardins, First moments of a polyhedron clipped by a paraboloid, ArXiv:2210.07772 (2022).

## Acknowledgements
The development of paraboloid-polytope intersection in IRL has directly benefitted from research funding provided by the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 101026017.