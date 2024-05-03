
# Sparse matrix
## Description
The code implements a class to store a sparse matrix in some convenient ways, not considering zero elements.
In particular we use COOmap for the uncompressed storing and CSR or CSC for the compressed one. Both of them are implemented row-wise and column-wise and the choice of which of these two methods to use is already given in the code by the programmer.

For row-wise implementation, COOmap uses a map whose key is the pair (row, column), which is the position of the non-zero element, and whose value is the value of this element.
For column-wise implementation it's the same but now the key is the pair (column, row) of the element.

CSR uses three vector: the first one (inner) contains the starting index for the elements of each row, the second one (outer) contains the corresponding column index and the last one contains the values of the non-zero elements. This implementation is again row-wise.

CSC it's the same but it's the respective column-wise: now inner contains the starting index for the elements of each column and outer ontains the corresponding row index.

It's important to take care that indexes of the matrix start from 0 both for rows and columns!

## Implementation
The class, called Matrix, contains as private members the number of rows and columns of the matrix, the number of non-zero elements, the map containing the values for COOmap implementation, the three vectors for CSR and CSC implementation and a boolean variable which tells if the matrix is compressed or not.

As public members there's a constructor, which takes as input the size of the matrix, and some functions. The first one resizes the matrix given the new size.
Then there are two similar operators which return the value and a pointer to an element of the matrix respectively, given its indexes.
Then there is a function which compresses the matrix,another one which decompresses it and one which returns a boolean which tells if the matrix is compressed or not.
Then there's a friend function which computes the product between a matrix and a vector given as input, and at last there's one which takes as input a MatrixMarket file and creates a matrix of this class with its content.
## Example
All the functions needed are contained into some .hpp files (one for the class, another one for the implementation of its functions and the last one for chrono's functions).
Then there's a Makefile which has to be used to compile the code by using the command make (or make doc for doxygen).
By running the code with the command ./main it's possible to see an example with a 131x131 matrix of double from MatrixMarket: this matrix is used for fluid flow modelling. Pay attention that the indexes in this file starts from 1, so they have to be scaled to construct the sparse matrix.