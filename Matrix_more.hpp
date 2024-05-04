//
// Created by Giulia Colombo on 02/05/24.
//

#ifndef CHALLENGE2_MATRIX_MORE_HPP
#define CHALLENGE2_MATRIX_MORE_HPP

#include "Matrix.hpp"

namespace algebra {

    /*!
     * @brief It changes the size of the matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @param n_r is the new number of rows
     * @param n_c is the new number of columns
     */
    template <typename T, typeOrder StorageOrder>
    void Matrix<T, StorageOrder>::resize (std::size_t n_r, std::size_t n_c) {

        if (compressed == 0) {

            if (n_r >= n_rows && n_c >= n_cols) {
                //Increase the number of elements of the matrix (they're automatically 0 since they're not put into the map)
                n_rows = n_r;
                n_cols = n_c;
            } else {
                n_rows = n_r;
                n_cols = n_c;

                //Delete the elements of the matrix which exceed the new size
                if constexpr (StorageOrder == typeOrder::rowWise) {
                    for (const auto& element : data) {
                        if (element.first[0] >= n_rows || element.first[1] >= n_cols) {
                            data.erase(element.first);
                            //Decrement the number of non zero
                            --nz;
                        }
                    }
                } else if constexpr (StorageOrder == typeOrder::columnWise) {
                    for (const auto& element : data) {
                        if (element.first[1] >= n_rows || element.first[0] >= n_cols) {
                            data.erase(element.first);
                            //Decrement the number of non zero
                            --nz;
                        }
                    }
                } else {
                    std::cerr << "Type of ordering not recognised!" << std::endl;
                }
            }
        } else {
            std::cerr << "Matrix compressed!" << std::endl;
        }
    }

    /*!
     * @brief It extract an element of the matrix (read only)
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @param row is the row of the element that has to be extracted
     * @param col is the column of the element that has to be extracted
     * @return the element required
     */
    template <typename T, typeOrder StorageOrder>
    T Matrix<T, StorageOrder>::operator() (std::size_t row, std::size_t col) const {

        T value = 0;

        if (row >= n_rows || col >= n_cols) {
            std::cerr << "Indexes out of range!" << std::endl;
            return value;
        }

        //Search the element with the input indexes and, if there exist, return its value otherwise return 0
        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (compressed == 0) {
                if (data.find({row, col}) != data.cend()) {
                    return data.at({row, col});;
                }
            } else {
                std::size_t i = inner[row];
                while (i >= inner[row] && i < inner[row + 1]) {
                    if (outer[i] == col) {
                        return values[i];
                    }
                    ++i;
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            if (compressed == 0) {
                if (data.find({col, row}) != data.cend()) {
                    return data.at({col, row});;
                }
            } else {
                std::size_t i = inner[col];
                while (i >= inner[col] && i < inner[col + 1]) {
                    if (outer[i] == row) {
                        return values[i];
                    }
                    ++i;
                }
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }

        return value;
    }

    /*!
     * @brief It extract an element of the matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @param row is the row of the element that has to be extracted
     * @param col is the column of the element that has to be extracted
     * @return the element required
     */
    template <typename T, typeOrder StorageOrder>
    T& Matrix<T, StorageOrder>::operator() (std::size_t row, std::size_t col) {

        if (row >= n_rows || col >= n_cols) {
            std::cerr << "Indexes out of range!" << std::endl;
            return data.end()->second;
        }

        //Search the element with the input indexes and, if there exist, return its reference to be able to change its value
        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (compressed == 0) {
                if (data.find({row, col}) == data.cend()) {
                    //Increment the number of non zero elements
                    ++nz;
                    return data[{row, col}];
                } else {
                    std::cout << "Element already present, this changes its value" << std::endl;
                    return data[{row, col}];
                }
            } else {
                std::size_t i = inner[row];
                while (i >= inner[row] && i < inner[row + 1]) {
                    if (outer[i] == col) {
                        return values[i];
                    }
                    ++i;
                }
                //If there's not this element it is added by uncompressing the matrix
                uncompress();
                operator()(row, col);
                compress();
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            if (compressed == 0) {
                if (data.find({col, row}) == data.cend()) {
                    //Increment the number of non zero elements
                    ++nz;
                    return data[{col, row}];
                } else {
                    std::cout << "Element already present, this changes its value" << std::endl;
                    return data[{col, row}];
                }
            } else {
                std::size_t i = inner[col];
                while (i >= inner[col] && i < inner[col + 1]) {
                    if (outer[i] == row) {
                        return values[i];
                    }
                    ++i;
                }
                //If there's not this element it is added by uncompressing the matrix
                uncompress();
                operator()(row, col);
                compress();
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }

        return data.end()->second;
    }

    /*!
     * @brief It compresses the matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     */
    template <typename T, typeOrder StorageOrder>
    void Matrix<T, StorageOrder>::compress() {

        if (compressed == 1) {
            std::cout << "Matrix already compressed!" << std::endl;
            return;
        }

        compressed = true;

        if constexpr (StorageOrder == typeOrder::rowWise) {
            //Resize the vectors needed
            inner.resize(n_rows + 1, 0);
            outer.resize(nz, 0);
            values.resize(nz, 0);
            std::size_t first = 0;

            for (std::size_t i = 0; i < inner.size(); ++i) { 
                inner[i] = first;
                for (std::size_t j = 0; j < n_cols; ++j) {
                    if (data.find({i, j}) != data.cend()) {
                        outer[first] = j;
                        values[first] = data[{i, j}];
                        ++first;
                    }
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            //Resize the vectors needed
            inner.resize(n_cols + 1, 0);
            outer.resize(nz, 0);
            values.resize(nz, 0);
            std::size_t first = 0;

            for (std::size_t i = 0; i < inner.size(); ++i) {
                inner[i] = first;
                for (std::size_t j = 0; j < n_rows; ++j) {
                    if (data.find({i, j}) != data.cend()) {
                        outer[first] = j;
                        values[first] = data[{i, j}];
                        ++first;
                    }
                }
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }
        data.clear();
    }

    /*!
     * @brief It uncompresses the matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     */
    template <typename T, typeOrder StorageOrder>
    void Matrix<T, StorageOrder>::uncompress() {

        if (compressed == 0) {
            std::cout << "Matrix already uncompressed!" << std::endl;
            return;
        }

        compressed = false;
        std::size_t idx = 0;

        if constexpr (StorageOrder == typeOrder::rowWise || StorageOrder == typeOrder::columnWise) {

            for (std::size_t i = 0; i < outer.size(); ++i) {
                if (i <= inner[idx]) {
                    data[{idx, i}] = values[i];
                } else {
                    ++idx;
                    data[{idx, i}] = values[i];
                }
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }
        inner.clear();
        outer.clear();
        values.clear();
    }

    /*!
     * @brief It tells if the matrix is compressed or not
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @return a boolean (0 if the matrix is uncompresse, 1 otherwise)
     */
    template <typename T, typeOrder StorageOrder>
    bool Matrix<T, StorageOrder>::is_compressed() const {
        return compressed;
    }

    /*!
     * @brief It does the matrix - vector product
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @param M is the matrix that I want to multiply
     * @param v is the vector that I want to multiply
     * @return the product between the two inputs
     */
    template <typename T, typeOrder StorageOrder>
    std::vector<T> operator* (const Matrix<T, StorageOrder>& M, const std::vector<T>& v) {

        std::vector<T> result(M.n_rows, 0);

        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (M.compressed == 0) {
                for (const auto& element : M.data) {
                    result[element.first[0]] += element.second * v[element.first[1]];
                }
            } else {
                for (std::size_t i = 0; i < M.n_rows; ++i) {
                    for (std::size_t j = M.inner[i]; j < M.inner[i + 1]; ++j) {
                        result[i] += M.values[j] * v[M.outer[j]];
                    }
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            if (M.compressed == 0) {
                for (const auto& element : M.data) {
                    result[element.first[1]] += element.second * v[element.first[0]];
                }
            } else {
                for (std::size_t i = 0; i < M.n_rows; ++i) {
                    for (std::size_t j = 0; j < M.n_cols; ++j) {
                        for (std::size_t k = M.inner[j]; k < M.inner[j + 1]; ++k) {
                            if (k == i) {
                                result[i] += M.values[k] * v[j];
                            }
                        }
                    }
                }
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }

        return result;
    }

    /*!
     * @brief It reads a file containing the data of a matrix and it transforms them into a sparse matrix of class Matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @param filename is the name of the file that the function has to read
     */
    template <typename T, typeOrder StorageOrder>
    void Matrix<T, StorageOrder>::read (const std::string& filename) {

        //Create the file
        std::ifstream infile(filename);

        if (!infile.is_open()) {
            std::cerr << "File not open!" << std::endl;
            return;
        }

        std::string s;
        getline(infile, s);

        //Verify that it's a MatrixMarket file by checking thee fist line
        if (s != "%%MatrixMarket matrix coordinate real general") {
            std::cerr << "Matrix not in MatrixMarket" << std::endl;
            return;
        }

        //Go to the fist line containing the data
        while (getline(infile, s)) {
            if (s[0] != '%') {
                break;
            }
        }

        std::size_t r;
        std::size_t c;
        T v;

        if constexpr (StorageOrder == typeOrder::rowWise) {

            std::istringstream ss(s);
            //Extract the number of rows and columns from the file
            ss >> n_rows >> n_cols >> nz;

            for (std::size_t i = 0; i < nz; i++) {
                getline(infile, s);
                std::istringstream sss(s);
                //Extract the values from the file
                sss >> r >> c >> v;
                data[{r - 1, c - 1}] = v;  //-1 to scale the indexes of the file (they start from 1, not 0)
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            
            std::istringstream ss(s);
            //Extract the number of rows and columns from the file
            ss >> n_cols >> n_rows >> nz;

            for (std::size_t i = 0; i < nz; i++) {
                 getline(infile, s);
                std::istringstream sss(s);
                //Extract the values from the file
                sss >> c >> r >> v;
                data[{c - 1, r - 1}] = v;  //-1 to scale the indexes of the file (they start from 1, not 0)
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }
        
        infile.close();
    }

    /*!
     * @brief It evaluates the norm of a matrix
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     * @tparam Norm is the type of norm that has to be evaluated
     * @return the norm
     */
    template <typename T, typeOrder StorageOrder>
    template<typeNorm Norm>
    T Matrix<T, StorageOrder>::norm() {

        T result;
        T sum;

        if constexpr (Norm == typeNorm::one) {

            if (compressed == 0) {

                if constexpr (StorageOrder == typeOrder::rowWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_cols; ++i) {
                        sum = 0;
                        for (std::size_t j = 0; j < n_rows; ++j) {
                            if (data.find({j, i}) != data.cend()) {
                                sum += abs(data[{j, i}]);
                            }
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }

                } else if constexpr (StorageOrder == typeOrder::columnWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_cols; ++i) {
                        sum = 0;
                        for (std::size_t j = 0; j < n_rows; ++j) {
                            if (data.find({i, j}) != data.cend()) {
                                sum += abs(data[{i, j}]);
                            }
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }

                } else {
                    std::cerr << "Type of ordering not recognised!" << std::endl;
                }

            } else {
                if constexpr (StorageOrder == typeOrder::rowWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_rows; ++i) {
                        for (std::size_t j = 0; j < n_cols; ++j) {
                            sum = 0;
                            for (std::size_t k = inner[j]; k < inner[j + 1]; ++k) {
                                if (k == i) {
                                    sum += abs(values[k]);
                                }
                            }

                            if (sum > result) {
                                result = sum;
                            }
                        }
                    }
                } else if constexpr (StorageOrder == typeOrder::columnWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_cols; ++i) {
                        sum = 0;
                        for (std::size_t j = inner[i]; j < inner[i + 1]; ++j) {
                            sum += abs(values[j]);
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }
                } else {
                    std::cerr << "Type of ordering not recognised!" << std::endl;
                }
            }

        } else if constexpr (Norm == typeNorm::infinty) {

            if (compressed == 0) {

                if constexpr (StorageOrder == typeOrder::rowWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_rows; ++i) {
                        sum = 0;
                        for (std::size_t j = 0; j < n_cols; ++j) {
                            if (data.find({i, j}) != data.cend()) {
                                sum += abs(data[{i, j}]);
                            }
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }
                } else if constexpr (StorageOrder == typeOrder::columnWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_rows; ++i) {
                        sum = 0;
                        for (std::size_t j = 0; j < n_cols; ++j) {
                            if (data.find({j, i}) != data.cend()) {
                                sum += abs(data[{j, i}]);
                            }
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }
                } else {
                    std::cerr << "Type of ordering not recognised!" << std::endl;
                }

            } else {
                if constexpr (StorageOrder == typeOrder::rowWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_rows; ++i) {
                        sum = 0;
                        for (std::size_t j = inner[i]; j < inner[i + 1]; ++j) {
                            sum += abs(values[j]);
                        }
                        if (sum > result) {
                            result = sum;
                        }
                    }
                } else if constexpr (StorageOrder == typeOrder::columnWise) {
                    result = 0;

                    for (std::size_t i = 0; i < n_rows; ++i) {
                        for (std::size_t j = 0; j < n_cols; ++j) {
                            sum = 0;
                            for (std::size_t k = inner[j]; k < inner[j + 1]; ++k) {
                                if (k == i) {
                                    sum += abs(values[k]);
                                }
                            }

                            if (sum > result) {
                                result = sum;
                            }
                        }
                    }
                } else {
                    std::cerr << "Type of ordering not recognised!" << std::endl;
                }
            }

        } else if constexpr (Norm == typeNorm::frobenius) {

            if (compressed == 0) {
                sum = 0;

                for (auto& element : data) {
                    sum += abs(element.second) * abs(element.second);
                }
                result = sqrt(sum);

            } else {
                sum = 0;

                for (std::size_t i = 0; i < values.size(); ++i) {
                    sum += abs(values[i]) * abs(values[i]);
                }
                result = sqrt(sum);
            }
        } else {
            std::cerr << "Type of norm not recognised!" << std::endl;
        }

        return result;
    }

}

#endif //CHALLENGE2_MATRIX_MORE_HPP
