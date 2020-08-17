#include <chrono>
#include <iostream>
#include <string>

#include "examples/advector/circle_rotation_2d.h"
#include "examples/advector/deformation_2d.h"
#include "examples/advector/reconstruction_types.h"
#include "examples/advector/solver.h"
#include "examples/advector/vof_advection.h"

static int startSimulation(const std::string& a_simulation_type,
                           const std::string& a_advection_method,
                           const std::string& a_reconstruction_method,
                           const double a_time_step_size,
                           const double a_time_duration,
                           const int a_viz_frequency);

int main(int argc, char* argv[]) {
  if (argc != 7) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout
        << "Simulation to run. Options: Deformation2D, CircleRotation2D\n";
    std::cout << "Advection method. Options: SemiLagrangian, "
                 "SemiLagrangianCorrected, FullLagrangian\n";
    std::cout << "Reconstruction method. Options: ELVIRA2D, LVIRA2D, MOF2D, "
                 "AdvectedNormals, R2P2D\n";
    std::cout << "Time step size, dt (double)\n";
    std::cout << "Simulation duration(double)\n";
    std::cout
        << "Amount of time steps between visualization output (integer)\n";
    std::exit(-1);
  }

  std::string simulation_type = argv[1];
  std::string advection_method = argv[2];
  std::string reconstruction_method = argv[3];
  double time_step_size = std::stod(argv[4]);
  double time_duration = std::stod(argv[5]);
  int viz_frequency = atoi(argv[6]);

  auto start = std::chrono::system_clock::now();
  startSimulation(simulation_type, advection_method, reconstruction_method,
                  time_step_size, time_duration, viz_frequency);
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
                           const int a_viz_frequency) {
  if (a_simulation_type == "CircleRotation2D") {
    return runSimulation<CircleRotation2D>(
        a_advection_method, a_reconstruction_method, a_time_step_size,
        a_time_duration, a_viz_frequency);
  } else if (a_simulation_type == "Deformation2D") {
    return runSimulation<Deformation2D>(
        a_advection_method, a_reconstruction_method, a_time_step_size,
        a_time_duration, a_viz_frequency);
  } else {
    std::cout << "Unknown simulation type of : " << a_simulation_type << '\n';
    std::cout << "Value entries are: CircleRotation2D, Deformation2D. \n";
    std::exit(-1);
  }
  return -1;
}
