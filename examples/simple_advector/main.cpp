#include <chrono>
#include <iostream>
#include <string>

#include "examples/simple_advector/solver.h"

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cout << "Incorrect amount of command line arguments supplied. \n";
    std::cout << "Arguments should be \n";
    std::cout << "Number of cells per dimension, N(integer) \n";
    std::cout << "Amount of circle rotations(integer)\n";
    std::cout << "Number of visualization output per rotation(integer)\n";
    std::exit(-1);
  }

  int number_of_cells = atoi(argv[1]);
  int number_of_rotations = atoi(argv[2]);
  int viz_output_freq = atoi(argv[3]);

  auto start = std::chrono::system_clock::now();
  runSimulation(number_of_cells, number_of_rotations, viz_output_freq);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = end - start;
  printf("Total run time: %20f \n\n", runtime.count());

  return 0;
}
