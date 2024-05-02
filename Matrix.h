//
// Created by Giulia Colombo on 19/04/24.
//

#ifndef CHALLENGE2_MATRIX_H
#define CHALLENGE2_MATRIX_H

#include "mmio.h"

#include <map>
#include <array>
#include <vector>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>

namespace algebra {

    /*!
     * @brief It contains all the possible type of storage ordering
     */
    enum class typeOrder{rowWise, columnWise};

    /*!
     * @brief It represents sparse matrices in uncompressed and compressed way
     * @tparam T is the type of elements
     * @tparam StorageOrder is the storage order of the matrix
     */
    template <typename T, typeOrder StorageOrder = typeOrder::rowWise> class Matrix {

    private:
        std::size_t n_rows = 0;
        std::size_t n_cols = 0;
        std::size_t nz = 0;
        std::map<std::array<std::size_t, 2>, T> data;
        std::vector<std::size_t> inner;
        std::vector<std::size_t> outer;
        std::vector<T> values;
        bool compressed = false;

    public:
        Matrix(std::size_t n_r, std::size_t n_c): n_rows(n_r), n_cols(n_c) {};

        void resize (std::size_t n_r, std::size_t n_c);

        T operator() (std::size_t row, std::size_t col) const;

        T& operator() (std::size_t row, std::size_t col);

        void compress();

        void uncompress();

        bool is_compressed() const;

        template <typename t, typeOrder so>
        friend std::vector<t> operator* (const Matrix<t, so>& M, const std::vector<t>& v);

        //Matrix<T, StorageOrder> read (const std::string& filename);

        void read (const std::string& filename);
    };
}

#endif //CHALLENGE2_MATRIX_H