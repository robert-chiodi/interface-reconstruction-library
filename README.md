# Interface Reconstruction Library (IRL)
The Interface Reconstruction Library (IRL) is a library of computational algorithms geared towards enabling easy and fast implementations of geometric Volume of Fluid (VOF) schemes for simulating multiphase flows. The methods contained in IRL, however, are general and could be used to compute volume-volume intersections, interface reconstructions, or perform optimizations. IRL was also specifically designed to be agnostic of the underlying meshes, making it useable for structured, unstructured, and adaptively refined meshes.

IRL is open-sourced under the Mozilla Public License 2 (MPL2) and hosted on GitLab, available here. It has already been used in incompressible and compressible CFD codes on structured, unstructured, and AMReX based flow solvers, with some images of simulations performed using IRL shown below.

I am currently in the process of writing additional documentation and examples. If you would like to become involved in its development or inquire about potentially using it for your application, please contact me.


IRL is written in C++14 with a Fortran 2008 interface. Below are links to additional information on how to install and use IRL:

- [Compilation & Installation](docs/markdown/install_main_page.md)

- [C/Fortran Interface Usage & Description](docs/markdown/interface_main_page.md)

- [Brief Introduction](docs/reference_powerpoint.pdf)

- [Examples](docs/markdown/examples_main_page.md)

