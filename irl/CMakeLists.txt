# List of directories which will add their own files.
set(IRL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(data_structures)
add_subdirectory(distributions)
add_subdirectory(generic_cutting)
add_subdirectory(geometry)
add_subdirectory(graphs)
add_subdirectory(helpers)
add_subdirectory(interface_reconstruction_methods)
add_subdirectory(moments)
add_subdirectory(optimization)
add_subdirectory(parameters)
add_subdirectory(planar_reconstruction)
add_subdirectory(paraboloid_reconstruction)
add_subdirectory(surface_mesher)

# Other directories for the C/Fortran interface.

add_subdirectory(c_interface)
add_subdirectory(fortran_interface)
