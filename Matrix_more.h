//
// Created by Giulia Colombo on 02/05/24.
//

#ifndef CHALLENGE2_MATRIX_MORE_H
#define CHALLENGE2_MATRIX_MORE_H

#include "Matrix.h"

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
                n_rows = n_r;
                n_cols = n_c;
            } else {
                n_rows = n_r;
                n_cols = n_c;
                for (const auto& element : data) {
                    if (element.first[0] >= n_rows || element.first >= n_cols) {
                        data.remove(element.first);
                    }
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

        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (compressed == 0) {
                if (data.find({row, col}) != data.cend()) {
                    return data.at({row, col});;
                } else if (row >= n_rows || col >= n_cols) {
                    std::cerr << "Indexes out of range!" << std::endl;
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
                } else if (row >= n_rows || col >= n_cols) {
                    std::cerr << "Indexes out of range!" << std::endl;
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
        }

        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (compressed == 0) {
                if (data.find({row, col}) == data.cend()) {
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
                //se non c'era già
                uncompress();
                operator()(row, col);
                compress();
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            if (compressed == 0) {
                if (data.find({col, row}) == data.cend()) {
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
                //se non c'era già
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
            inner.resize(n_rows + 1, 0);
            std::size_t first = 0;

            for (std::size_t i = 0; i < inner.size(); ++i) {  //se una riga ha tutti 0?
                inner[i] = first;
                for (std::size_t j = 0; j < n_cols; ++j) {
                    if (data.find({i, j}) != data.cend()) {
                        outer.push_back(j);
                        values.push_back(data[{i, j}]);
                        ++first;
                    }
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            inner.resize(n_cols + 1, 0);
            std::size_t first = 0;

            for (std::size_t i = 0; i < inner.size(); ++i) {
                inner[i] = first;
                for (std::size_t j = 0; j < n_rows; ++j) {
                    if (data.find({i, j}) != data.cend()) {
                        outer.push_back(j);
                        values.push_back(data[{i, j}]);
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

        if constexpr (StorageOrder == typeOrder::rowWise) {
            std::size_t row = 0;

            for (std::size_t i = 0; i < outer.size(); ++i) {
                if (i <= inner[row]) {
                    data[{row, i}] = values[i];
                } else {
                    ++row;
                    data[{row, i}] = values[i];
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {  //è uguale, non faccio if?
            std::size_t col = 0;

            for (std::size_t i = 0; i < outer.size(); ++i) {
                if (i <= inner[col]) {
                    data[{col, i}] = values[i];
                } else {
                    ++col;
                    data[{col, i}] = values[i];
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

        std::vector<T> result(v.size(), 0);

        if constexpr (StorageOrder == typeOrder::rowWise) {
            if (M.compressed == 0) {
                for (const auto& element : M.data) {
                    result[element.first[0]] += element.second * element.first[1];
                }
            } else {
                for (std::size_t i = 0; i < v.size(); ++i) {
                    for (std::size_t j = 0; j < M.n_cols; ++j) {
                        result[i] += M(i, j) * v[j];
                    }
                }
            }
        } else if constexpr (StorageOrder == typeOrder::columnWise) {
            if (M.compressed == 0) {
                for (const auto& element : M.data) {
                    result[element.first[1]] += element.second * element.first[0];
                }
            } else {
                for (std::size_t i = 0; i < v.size(); ++i) {
                    for (std::size_t j = 0; j < M.n_cols; ++j) {
                        result[i] += M(i, j) * v[j];
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

        std::ifstream infile(filename);

        if (!infile.is_open()) {
            std::cerr << "File not open!" << std::endl;
            return;
        }

        std::string s;
        getline(infile, s);

        if (s != "%%MatrixMarket matrix coordinate real general") {
            std::cerr << "Matrix not in MatrixMarket" << std::endl;
            return;
        }

        /*while (getline(infile, s)) {
            if (s[0] != '%') {
                break;
            }
        }*/

        //Go to the second line
        getline(infile, s);

        std::istringstream ss(s);
        ss >> n_rows >> n_cols >> nz;

        std::size_t r;
        std::size_t c;
        T v;

        for (std::size_t i = 0; i < nz; i++) {
            ss >> r >> c >> v;
            data[{r, c}] = v;
        }

        infile.close();
    }

}

#endif //CHALLENGE2_MATRIX_MORE_H
