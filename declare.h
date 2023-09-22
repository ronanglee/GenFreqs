#include <string>
#include <utility>

using namespace std;
 
class InputGeometry {
  public:

    int *zvals;
    int natom;

    InputGeometry(const char *filename);
    ~InputGeometry();
};

class Hessian {
  public:

    int hess_natom;
    double **hessian_mat;

    Hessian(const char *filename, int *zvals);
    ~Hessian();

};

bool endsWithExtension(const std::string& str, const std::string& extension);
double EigenValueInCm();
