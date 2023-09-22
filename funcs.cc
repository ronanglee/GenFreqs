#include "declare.h"
#include "AtomicMasses.h"
#include "constants.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"
#include "eigen-3.4.0/Eigen/Core"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <string>


InputGeometry::InputGeometry(const char *filename)
{

  zvals = new int[natom];
  std::string line;
  int num = 0;
  int firstInt;

   // open filename
  std::ifstream is(filename);

  assert(is.good());

  // read the number of atoms from filename
  is >> natom;

  while (std::getline(is, line)) {
    std::istringstream iss(line);

    if (iss >> firstInt) {
        // Process the first integer of the line here
        zvals[num] = firstInt;
    } else {
        // Handle the case where there is no integer on the line
        continue;
    }
    ++num;
  }
  is.close();
}

InputGeometry::~InputGeometry()
{
  delete[] zvals;
}

Hessian::Hessian(const char *filename, int *zvals)
{
  // open filename
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

  std::ifstream is(filename);
  assert(is.good());

  // read the number of atoms from filename
  is >> hess_natom;

  // allocate space the array of pointers
  hessian_mat = new double* [hess_natom*3];

  for(int i=0; i < hess_natom*3; i++){
    hessian_mat[i] = new double[hess_natom*3];
}

  // read the data from filename
  for( int i=0; i < hess_natom*3; i++){
    for(int j=0; j < hess_natom*3; j++){ 
    is >> hessian_mat[i][j];
        }
  }
  is.close();

  float amasses[hess_natom];

// get atomic masses from AtomicMasses.h
  for (int i=0; i < hess_natom; i++){
    amasses[i] = atomicMasses[zvals[i]];
      }

// mass-weight hessian
  for(int i=0; i < hess_natom*3; i++){
    for(int j=0; j < hess_natom*3; j++){
        hessian_mat[i][j] /= sqrt(amasses[i/3]*amasses[j/3]);
            }
  }
}
  
Hessian::~Hessian()
{
  for(int i=0; i < hess_natom; i++)
    delete[] hessian_mat[i];
  delete[] hessian_mat;
}

bool endsWithExtension(const std::string& str, const std::string& extension) {
    // Check if str ends with the given extension
    if (str.length() >= extension.length()) {
        return (str.compare(str.length() - extension.length(), extension.length(), extension) == 0);
    } else {
        return false; // The string is shorter than the extension
    }
}

double EigenValueInCm(){

  double conversion;

  conversion = (1/(2*M_PI*Constants::SpeedOfLight))*sqrt(Constants::AvogadroNumber*Constants::HartreeToKilojoules)*sqrt(1/pow(Constants::BohrToCMeters,2));
  return conversion;

}
  