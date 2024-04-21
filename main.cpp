#include "Matrix.h"
#include "chrono.h"

#include <iostream>

using namespace algebra;
using namespace Timings;

int main() {

    Matrix<double, typeOrder::rowWise> matrix = read<double, typeOrder::rowWise>("lnsp_131.mtx");
    std::vector<double> vector(131, 0.);
    vector[0] = 1.;

    Chrono time;

    time.start();
    std::vector<double> result = matrix * vector;
    std::cout << "Elapsed time " << time.wallTimeNow() << " microseconds" << std::endl;

    matrix.compress();
    result = matrix * vector;
    time.stop();
    std::cout << time;

    return 0;
}