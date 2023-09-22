#include "declare.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"
#include "eigen-3.4.0/Eigen/Core"
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>
#include <list>

int main(int argc, char* argv[])
{ 

  int num_modes;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

  // check for correct number of arguments and correct extension
  assert(argc == 3 && "Program needs a .geom and .txt file as arguments.");

  std::list<std::string> extensions = {".dat", ".txt"};
   
  for (int i = 1; i < argc; ++i) {
      std::string argument(argv[i]);
      std::string desiredExtension = extensions.front();
      extensions.pop_front();
      if (endsWithExtension(argument, desiredExtension)) {
          continue;
      } else {
          std::cout << "Argument" << argument << " does not have the correct extension." << std::endl;
          exit(1);
      }
  }

  // read in geometry and hessian
  InputGeometry mol(argv[1]);
  Hessian hess(argv[2], mol.zvals);

  // write to output file
  num_modes = 3*hess.hess_natom - 6;

  Matrix eigen_mat(hess.hess_natom*3, hess.hess_natom*3);

  for (int i=0; i < hess.hess_natom*3; i++){
    for (int j=0; j < hess.hess_natom*3; j++){
      eigen_mat(i,j) = hess.hessian_mat[i][j];
    }
  }

  Eigen::SelfAdjointEigenSolver<Matrix> es(eigen_mat);
  Vector evals = es.eigenvalues();
  Matrix evecs = es.eigenvectors();
  std::ofstream myfile;
  myfile.open ("freqs.txt");
  myfile << "Frequencies in cm^-1\n";
  for (int i=0; i < num_modes; i++){
    std::cout << evals << "\n";
    myfile << sqrt(evals(evals.size()-1-i))*EigenValueInCm() << "\n";
  }
  myfile.close();
  
  return 0;
}