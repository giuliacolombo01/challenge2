#include "Matrix.hpp"
#include "Matrix_more.hpp"
#include "chrono.hpp"

#include <iostream>

using namespace algebra;
using namespace Timings;

int main() {

    Matrix<double, typeOrder::rowWise> matrix(0, 0);
    matrix.read("lnsp_131.mtx");
    std::vector<double> vector(131, 0.);
    vector[0] = 1.;

    Chrono time;

    time.start();
    std::vector<double> result = matrix * vector;
    std::cout << "Elapsed time for uncompressed matrix" << time.wallTimeNow() << " microseconds" << std::endl;

    std::cout << "The result is: [";
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << "]^T" << std::endl;

    matrix.compress();
    result = matrix * vector;
    time.stop();
    std::cout << "Elapsed time for compressed matrix" << time << " microseconds" << std::endl;

    std::cout << "The result is: [";
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << "]^T" << std::endl;

    return 0;
}