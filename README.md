
# Sparse matrix

The code implements a class to store a sparse matrix in some convenient ways, not considering zero elements.
In particular we use COOmap for the uncompressed storing and CSR or CSC for the compressed one. Both of them are implemented row-wise and column-wise.

For row-wise implementation, COOmap uses a map whose key is the pair (row, column), which is the position of the non-zero element, and whose value is the value of this element.
For column-wise implementation it's the same but now the key is the pair (column, row) of the element.

CSR uses 3 vector: the first one (inner) contains the starting index for the elements of each row, the second one (outer) contains the corresponding column index and the last one contains the values of the non-zero elements. This implementation is again row-wise.

CSC it's the same but it's the respective column-wise: now inner contains the starting index for the elements of each column and outer ontains the corresponding row index.llenges