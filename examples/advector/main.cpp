#include <chrono>
#include <iostream>
#include <string>

#include "examples/advector/circle_rotation_2d.h"
#include "examples/advector/deformation_2d.h"
#include "examples/advector/deformation_3d.h"
#include "examples/advector/reconstruction_types.h"
#include "examples/advector/solver.h"
#include "examples/advector/vof_advection.h"

static int startSimulation(const std::string& a_simulation_type,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency,
                           const std::string& a_viz_type);

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout
        << "Simulation to run. Options: Deformation2D, Deformation3D, CircleRotation2D\n";
    std::cout << "Advection method. Options: SemiLagrangian, "
                 "SemiLagrangianCorrected, FullLagrangian\n";
    std::cout << "Reconstruction method. Options: ELVIRA2D, LVIRA2D, MOF2D, "
                 "AdvectedNormals, R2P2D\n";
    std::cout << "Time step size, dt (double)\n";
    std::cout << "Number of cycles (integer)\n";
    std::cout
        << "Amount of time steps between visualization output (integer)\n";
    std::cout
        << "Type of visualization files. Options: VTK, VIZ (to use visualize.m)\n";
    std::exit(-1);
  }

  std::string simulation_type = argv[1];
  std::string advection_method = argv[2];
  std::string reconstruction_method = argv[3];
  double time_step_size = std::stod(argv[4]);
  int n_cycles = std::atoi(argv[5]);
  int viz_frequency = atoi(argv[6]);
  std::string viz_type = argv[7];

  double time_duration;
  if (simulation_type == "Deformation2D") {
    time_duration = static_cast<double>(n_cycles) * 8.0;
  } else if (simulation_type == "Deformation3D") {
    time_duration = static_cast<double>(n_cycles) * 3.0;
  } else if (simulation_type == "CircleRotation2D") {
    time_duration = static_cast<double>(n_cycles) * 1.0;

  } else {
    std::cout << "Unknown simulation type of " << simulation_type << std::endl;
    return -1;
  }

  if (viz_type != "VTK" && viz_type != "VIZ") {
    std::cout << "Unknown visualization type of " << viz_type << std::endl;
    return -1;
  }

  auto start = std::chrono::system_clock::now();
  startSimulation(simulation_type, advection_method, reconstruction_method,
                  time_step_size, time_duration, viz_frequency, viz_type);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;
  printf("Total run time: %20f \n\n", runtime.count());

  return 0;
}

static int startSimulation(const std::string& a_simulation_type,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency,
                           const std::string& a_viz_type) {
  if (a_simulation_type == "CircleRotation2D") {
    return runSimulation<CircleRotation2D>(
        a_advection_method, a_reconstruction_method, a_time_step_size,
        a_time_duration, a_viz_frequency, a_viz_type);
  } else if (a_simulation_type == "Deformation2D") {
    return runSimulation<Deformation2D>(
        a_advection_method, a_reconstruction_method, a_time_step_size,
        a_time_duration, a_viz_frequency, a_viz_type);
  } else if (a_simulation_type == "Deformation3D") {
    return runSimulation<Deformation3D>(
        a_advection_method, a_reconstruction_method, a_time_step_size,
        a_time_duration, a_viz_frequency, a_viz_type);
  } else {
    std::cout << "Unknown simulation type of : " << a_simulation_type << '\n';
    std::cout << "Value entries are: CircleRotation2D, Deformation2D. \n";
    std::exit(-1);
  }
  return -1;
}
