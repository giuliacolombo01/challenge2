//
// Created by Giulia Colombo on 19/04/24.
//

#ifndef CHALLENGE2_MATRIX_H
#define CHALLENGE2_MATRIX_H

#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <fstream>
#include <cstdio>
#include "mmio.h"

namespace algebra {

    enum class typeOrder{rowWise, columnWise};

    template <typename T, typeOrder StorageOrder = typeOrder::rowWise> class Matrix {

    private:
        std::size_t n_rows = 0;
        std::size_t n_cols = 0;
        std::map<std::array<std::size_t, 2>, T> data;
        std::vector<std::size_t> inner;
        std::vector<std::size_t> outer;
        std::vector<T> values;
        bool compressed = false;

    public:
        Matrix(std::size_t n_r, std::size_t n_c): n_rows(n_r), n_cols(n_c) {};

        void resize (std::size_t n_r, std::size_t n_c) {

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

        const T operator() (std::size_t row, std::size_t col) {

            T value = 0;

            if constexpr (StorageOrder == typeOrder::rowWise) {
                if (compressed == 0) {
                    if (data.find({row, col}) != data.cend()) {
                        return data[{row, col}];
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
                        return data[{col, row}];
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

        void operator() (std::size_t row, std::size_t col, T value) {

            if (row >= n_rows || col >= n_cols) {
                std::cerr << "Indexes out of range!" << std::endl;
                return;
            }

            if constexpr (StorageOrder == typeOrder::rowWise) {
                if (compressed == 0) {
                    if (data.find({row, col}) == data.cend()) {
                        data[{row, col}] = value;
                        return;
                    } else {
                        std::cout << "Element already present, this changes its value" << std::endl;
                        data[{row, col}] = value;
                        return;
                    }
                } else {
                    std::size_t i = inner[row];
                    while (i >= inner[row] && i < inner[row + 1]) {
                        if (outer[i] == col) {
                            values[i] = value;
                            return;
                        }
                        ++i;
                    }
                    //se non c'era già
                    uncompress();
                    operator()(row, col, value);
                    compress();
                }
            } else if constexpr (StorageOrder == typeOrder::columnWise) {
                if (compressed == 0) {
                    if (data.find({col, row}) == data.cend()) {
                        data[{col, row}] = value;
                        return;
                    } else {
                        std::cout << "Element already present, this changes its value" << std::endl;
                        data[{col, row}] = value;
                        return;
                    }
                } else {
                    std::size_t i = inner[col];
                    while (i >= inner[col] && i < inner[col + 1]) {
                        if (outer[i] == row) {
                            values[i] = value;
                            return;
                        }
                        ++i;
                    }
                    //se non c'era già
                    uncompress();
                    operator()(row, col, value);
                    compress();
                }
            } else {
                std::cerr << "Type of ordering not recognised!" << std::endl;
            }
        }

        void compress() {

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
            delete data;
        }

        void uncompress() {

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

        const bool is_compressed() {
            return compressed;
        }

        friend std::vector<T> operator* (const Matrix& M, const std::vector<T>& v);

        friend Matrix read (const std::string& filename);
    };

    template <typename T, typeOrder StorageOrder = typeOrder::rowWise>
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
                        result[i] += M.operator()(i, j) * v[j];
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
                        result[i] += M.operator()(i, j) * v[j];
                    }
                }
            }
        } else {
            std::cerr << "Type of ordering not recognised!" << std::endl;
        }

        return result;
    }

    template <typename T, typeOrder StorageOrder = typeOrder::rowWise>
    Matrix<T, StorageOrder> read(const std::string& filename) {

        MM_typecode matcode;
        std::ifstream infile(filename);
        std::size_t n_rows, n_cols, nz;
        std::vector<std::size_t> idx_rows, idx_cols;
        std::vector<T> val;

        if (mm_read_banner(&infile, &matcode) != 0) {
            throw std::runtime_error("Could not process Matrix Market banner.");
        }

        if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
            throw std::runtime_error("Sorry, this function does not support Matrix Market type: " + std::string(mm_typecode_to_str(matcode)));
        }

        if ((mm_read_mtx_crd_size(&infile, &n_rows, &n_cols, &nz)) != 0) {
            throw std::runtime_error("Error reading matrix size.");
        }

        idx_rows.resize(nz);
        idx_cols.resize(nz);
        val.resize(nz);

        Matrix<T, StorageOrder> matrix(n_rows, n_cols);

        for (std::size_t i = 0; i < nz; i++) {
            infile >> idx_rows[i] >> idx_cols[i] >> val[i];
            idx_rows[i]--;
            idx_cols[i]--;
            matrix.operator()(idx_rows[i], idx_cols[i], val[i]);
        }

        if (infile.is_open()) {
            infile.close();
        }

        return matrix;
    }
}


#endif //CHALLENGE2_MATRIX_H
