#include "Matrix.hpp"
#include "Matrix_more.hpp"
#include "chrono.hpp"

#include <iostream>

using namespace algebra;
using namespace Timings;

int main() {

    //Create a matrix by using the constructor and then initialize it by the function read
    Matrix<double, typeOrder::rowWise> matrix(0, 0);
    matrix.read("lnsp_131.mtx");

    //Create the vector to multiply
    std::vector<double> vector(131, 0.);
    vector[1] = 1.;

    Chrono time;

    //Do the product both with uncompressed and compressed matrixes and measure the time
    time.start();
    std::vector<double> result = matrix * vector;
    std::cout << "Elapsed time for uncompressed matrix: " << time.wallTimeNow() << " microseconds" << std::endl;

    std::cout << "The result is: [";
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << "]'" << std::endl;

    std::cout << std::endl;

    matrix.compress();
    result = matrix * vector;
    time.stop();
    std::cout << "Elapsed time for compressed matrix: " << time << " microseconds" << std::endl;

    std::cout << "The result is: [";
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << "]'" << std::endl;

    std::cout << std::endl;

    //Compute the three norms
    std::cout << "1-norm of the matrix is: " << matrix.norm<typeNorm::one>() << std::endl;
    std::cout << "Infinity norm of the matrix is: " << matrix.norm<typeNorm::infinty>() << std::endl;
    std::cout << "Frobenius norm of the matrix is: " << matrix.norm<typeNorm::frobenius>() << std::endl;

    return 0;
}