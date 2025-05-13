// $Source$
//----------------------------------------------------------------------
// matrix
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file matrix.cpp
 * @brief Implements methods with matrices and vectors.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\matrix.hpp"

/**
 * @brief Construct a new Matrix:: Matrix object
 * 
 */
Matrix::Matrix(){
	this->n_row=0;
	this->n_column=0;
	this->data = nullptr;
}
/**
 * @brief Construct a new Matrix:: Matrix object
 * 
 * @param n_row The number of rows of the matrix
 * @param n_column The number of columns of the matrix
 */

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column1\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

/**
 * @brief Construct a new Matrix:: vector object
 * 
 * @param v_size the size of the vector
 */

Matrix::Matrix(const int v_size){
	if (v_size < 0) {
		cout << "Vector create: error in v_size1\n";
        exit(EXIT_FAILURE);
	}
	this->n_row=1;
	this->n_column=v_size;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
	if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	this->data[0] = (double *) calloc(v_size,sizeof(double));
}

/**
 * @brief Access an element of the matrix by row and column
 * 
 * @param row The row index (1-based)
 * @param column The column index (1-based)
 * @return double& Reference to the element at the specified position
 */
double& Matrix::operator () (const int row, const int column) {

    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
        cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
    }

    return this->data[row - 1][column - 1];
}

/**
 * @brief Access an element of the matrix by vector size
 * 
 * @param v_size The index of the element in a 1D representation (1-based)
 * @return double& Reference to the element at the specified position
 */
double& Matrix::operator () (const int v_size) {

    if (v_size <= 0 || v_size > this->n_column*this->n_row) {
        cout << "Matrix get: error in column\n";
        exit(EXIT_FAILURE);
    }

return this->data[(v_size-1)/this->n_column][(v_size-1)%this->n_column];
}

/**
 * @brief Overload the addition operator for matrix addition
 * 
 * @param m The matrix to add
 * @return Matrix& Reference to the resulting matrix after addition
 */
Matrix& Matrix::operator + (Matrix &m) {
    if (this->n_row != m.n_row || this->n_column != m.n_column) {
        cout << "Matrix sum: error in n_row/n_column2\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i,j) = (*this)(i,j) + m(i,j);
        }
    }
    
    return *m_aux;
}

/**
 * @brief Overloads the subtraction operator for the Matrix class.
 * 
 * @param m The matrix to subtract from the current matrix.
 * @return Matrix& A reference to a new matrix containing the result of 
 * the subtraction.
 */
Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column3\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

/**
 * @brief Overloads the assignment operator for the Matrix class.
 * 
 * @param m Reference to the Matrix object to be assigned.
 * @return Reference to the current Matrix object after assignment.
 */
Matrix& Matrix::operator = (Matrix &m) {
	this->n_row = m.n_row;
	this->n_column = m.n_column;

	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix assignment: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < m.n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
		for (int j = 0; j < this->n_column; j++) {
			this->data[i][j]=m.data[i][j];
		}
	}
	
	return *this;
}

/**
 * @brief Overloads the multiplication operator to perform matrix multiplication.
 * 
 * @param matrix2 The matrix to multiply with the current matrix.
 * @return Matrix& A reference to the resulting matrix after multiplication.
 *         The result is dynamically allocated and returned by reference.
 */
Matrix& Matrix::operator * (Matrix& matrix2){
	if (this->n_column != matrix2.n_row) {
        cout << "Matrix multiplication error: number of columns of A must equal number of rows of B\n";
        exit(EXIT_FAILURE);
    }
	
	Matrix *m_aux = new Matrix(this->n_row, matrix2.n_column);	
	
	for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= matrix2.n_column; j++) {
            (*m_aux)(i, j) = 0;
            for (int k = 1; k <= this->n_column; k++) {
                (*m_aux)(i, j) += (*this)(i, k) * matrix2(k, j);
            }
        }
    }
	
	return *m_aux;
}

/**
 * @brief Overloads the division operator to perform matrix division.
 * 
 * This operator divides the current matrix by another matrix by multiplying
 * the current matrix with the inverse of the second matrix. It checks if the
 * dimensions of the matrices are compatible for multiplication.
 * 
 * @param matrix2 The matrix to divide the current matrix by.
 * @return A reference to the resulting matrix after division.
 * @note The function exits the program if the dimensions of the matrices are incompatible.
 */
Matrix& Matrix::operator / (Matrix& matrix2) {
    if (this->n_column != matrix2.n_row) {
        cout << "Matrix division: incompatible matrix dimensions for multiplication\n";
        exit(EXIT_FAILURE);
    }

    Matrix invMatrix2 = inv(matrix2);

    return *this * invMatrix2;  
}

/**
 * @brief Overloads the unary minus operator to negate all elements of the matrix.
 * 
 * This operator multiplies all elements of the matrix by -1.0, effectively
 * negating the matrix.
 * 
 * @return A reference to the negated matrix.
 */
Matrix& Matrix::operator-() {
	return (*this)*(-1.0);
}

/**
 * @brief Overloads the addition operator to add a scalar to each element of the matrix.
 * 
 * @param scalar The scalar value to be added to each element of the matrix.
 * @return Matrix& A reference to a new matrix containing the result of the addition.
 */
Matrix& Matrix::operator + (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + scalar;
        }
    }

    return *m_aux;
}

/**
 * @brief Overloads the subtraction operator to subtract a scalar from each element of the matrix.
 * 
 * @param scalar The scalar value to subtract from each element of the matrix.
 * @return Matrix& A reference to a new matrix containing the result of the subtraction.
 */
Matrix& Matrix::operator - (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - scalar;
        }
    }

    return *m_aux;
}

/**
 * @brief Overloads the multiplication operator to scale the matrix by a scalar value.
 * 
 * @param scalar The scalar value to multiply each element of the matrix by.
 * @return Matrix& A reference to a new matrix containing the scaled values.
 */
Matrix& Matrix::operator * (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) * scalar;
        }
    }

    return *m_aux;
}

/**
 * @brief Overloads the division operator to divide all elements of the matrix by a scalar.
 * 
 * @param scalar The scalar value by which each element of the matrix will be divided.
 *               Must not be zero or close to zero (less than 1e-10 in absolute value).
 * @return Matrix& A reference to a new matrix containing the result of the division.
 *                 Note: The returned matrix is dynamically allocated and must be managed
 *                 (e.g., deleted) by the caller to avoid memory leaks.
 */
Matrix& Matrix::operator / (double scalar) {
    if (abs(scalar) < 1e-10) {
        cout << "Matrix division by zero error\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) / scalar;
        }
    }

    return *m_aux;
}


/**
 * @brief Overloads the output stream operator to print the contents of a Matrix object.
 * 
 * @param o The output stream to which the matrix will be written.
 * @param m The Matrix object to be printed.
 * @return A reference to the output stream after writing the matrix.
 */
ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

/**
 * @brief Creates a matrix filled with zeros.
 * 
 * @param n_row Number of rows in the matrix.
 * @param n_column Number of columns in the matrix.
 * @return Matrix& Reference to the newly created matrix filled with zeros.
 *         Note: The returned reference points to a dynamically allocated object.
 *         The caller is responsible for managing its memory.
 */
Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

/**
 * @brief Creates a matrix initialized with zeros.
 * 
 * @param v_size The size of the matrix to be created.
 * @return A reference to the newly created matrix filled with zeros.
 */
Matrix& zeros(const int v_size) {
	Matrix *m_aux = new Matrix(v_size);
	
	for(int j = 1; j <= v_size; j++) {
		(*m_aux)(j) = 0;
	}
	
	return (*m_aux);
}

/**
 * @brief Creates an identity matrix of the specified size.
 * 
 * @param size The size of the identity matrix (number of rows and columns).
 * @return A reference to the created identity matrix.
 * 
 * @note The returned matrix is dynamically allocated and must be manually deleted 
 *       to avoid memory leaks.
 */
Matrix& eye(int size) {
    Matrix *m_aux = new Matrix(size, size);  

    for (int i = 1; i <= size; i++) {
        for (int j = 1; j <= size; j++) {
            (*m_aux)(i, j) = (i == j) ? 1 : 0;
        }
    }

    return *m_aux;
}

/**
 * @brief Transposes the given matrix.
 * 
 * @param m The matrix to be transposed. It is passed by reference.
 * @return Matrix& A reference to the newly created transposed matrix.
 *         Note: The returned reference points to a dynamically allocated
 *         matrix, and it is the caller's responsibility to manage its memory.
 */
Matrix& transpose(Matrix &m) {
    Matrix *m_aux = new Matrix(m.n_column, m.n_row); 

    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(j, i) = m(i, j);  
        }
    }

    return *m_aux;  
}

/**
 * @brief Computes the inverse of a square matrix.
 * 
 * This function calculates the inverse of the given square matrix using 
 * Gaussian elimination. If the matrix is not square or is singular, 
 * the function will terminate the program with an error message.
 * 
 * @param m The square matrix to be inverted. The matrix must have the 
 *          same number of rows and columns.
 * @return A reference to a new matrix object containing the inverse 
 *         of the input matrix.
 * 
 * @note The function dynamically allocates memory for the inverse matrix 
 *       and returns it. The caller is responsible for managing the 
 *       memory of the returned matrix to avoid memory leaks.
 */
Matrix& inv(Matrix& m) {
    if (m.n_row != m.n_column) {
        std::cout << "Matrix inverse: matrix must be square\n";
        exit(EXIT_FAILURE);
    }

    Matrix *inv = new Matrix(m);   
    Matrix *m_aux = &eye(m.n_column);  

    for (int i = 1; i <= m.n_row; i++) {
        double pivot = (*inv)(i, i);
        if (fabs(pivot) < 1e-10) {
            std::cout << "Matrix inverse: matrix is singular (no inverse)\n";
            exit(EXIT_FAILURE);
        }

        for (int j = 1; j <= m.n_column; j++) {
            (*inv)(i, j) /= pivot;
            (*m_aux)(i, j) /= pivot;
        }

        for (int k = 1; k <= m.n_row; k++) {
            if (k != i) {
                double factor = (*inv)(k, i);
                for (int j = 1; j <= m.n_column; j++) {
                    (*inv)(k, j) -= (*inv)(i, j) * factor;
                    (*m_aux)(k, j) -= (*m_aux)(i, j) * factor;
                }
            }
        }
    }

    return *m_aux;  
}

/**
 * @brief Computes the Euclidean norm of a matrix.
 * 
 * @param a Reference to a Matrix object for which the norm is calculated.
 * @return The Euclidean norm of the matrix as a double.
 */
double norm(Matrix& a){
	 double r=0;
	 for(int i=1;i<=a.n_column;i++){
		 r +=pow(a(i),2);
	 }
	 double solv=sqrt(r);
	 return solv;
}

/**
 * @brief Computes the dot product of two matrices.
 * 
 * @param a Reference to the first matrix.
 * @param b Reference to the second matrix.
 * @return The dot product of the two matrices as a double.
 */
double dot(Matrix& a,Matrix& b){
	 double r=0;
	 for(int i=1;i<=a.n_column;i++)
		 r +=a(i)*b(i);
	 return r;
}

/**
 * @brief Computes the cross product of two 3-dimensional matrices.
 * 
 * @param a The first matrix operand.
 * @param b The second matrix operand.
 * @return Matrix& A reference to a new matrix containing the result of the cross product.
 * 
 * @note The returned matrix is dynamically allocated and must be manually deleted to avoid memory leaks.
 */
Matrix& cross(Matrix& a,Matrix& b){
	Matrix *m_aux = new Matrix(3);
	(*m_aux)(1) = a(2) * b(3) - a(3) * b(2);  
    (*m_aux)(2) = a(3) * b(1) - a(1) * b(3);  
    (*m_aux)(3) = a(1) * b(2) - a(2) * b(1);  
    
    return (*m_aux);
}

/**
 * @brief Extracts a subvector from a matrix.
 * 
 * @param m The matrix from which the subvector will be extracted.
 * @param n1 The starting column index (1-based) of the subvector.
 * @param n2 The ending column index (1-based) of the subvector.
 * @return A reference to a newly allocated matrix containing the extracted subvector.
 * 
 * @note The function will terminate the program if:
 *       - n1 or n2 are less than or equal to 0.
 *       - n1 or n2 exceed the number of columns in the matrix.
 *       - n1 is greater than n2.
 */
Matrix& extract_vector(Matrix& m, const int n1, const int n2) {
    if (n1 <= 0 || n2 <= 0 || n1 > m.n_column || n2 > m.n_column || n1 > n2) {
        cout << "Matrix extract_vector: error in n1/n2\n";
        exit(EXIT_FAILURE);
    }

    int size = n2 - n1 + 1;
    Matrix *m_aux = new Matrix(size);  

    for (int i = 1; i <= size; i++) {
        (*m_aux)(i) = m(n1 + i - 1);
    }

    return *m_aux;
}

/**
 * @brief Combines two vector matrices into a single vector matrix.
 * 
 * This function takes two matrices, which must be vectors (either a single row or a single column),
 * and concatenates them into a new vector matrix.
 * 
 * @param m1 The first vector matrix to be combined.
 * @param m2 The second vector matrix to be combined.
 * @return A reference to the newly created vector matrix containing the elements of both input matrices.
 * 
 * @note The caller is responsible for managing the memory of the returned matrix, as it is dynamically allocated.
 * @warning If either of the input matrices is not a vector, the program will terminate with an error message.
 */
Matrix& union_vector(Matrix& m1, Matrix& m2) {
    if (m1.n_row != 1 && m1.n_column != 1) {
        cout << "Matrix union_vector: first matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }
    if (m2.n_row != 1 && m2.n_column != 1) {
        cout << "Matrix union_vector: second matrix must be a vector\n";
        exit(EXIT_FAILURE);
    }

    int total_size = m1.n_column + m2.n_column;
    Matrix *m_aux = new Matrix(total_size); 

    for (int i = 1; i <= m1.n_column; i++) {
        (*m_aux)(i) = m1(i);
    }

    for (int i = 1; i <= m2.n_column; i++) {
        (*m_aux)(m1.n_column + i) = m2(i);
    }

    return *m_aux;
}

/**
 * @brief Extracts a specific row from a matrix and returns it as a new matrix.
 * 
 * @param m The matrix from which the row will be extracted.
 * @param row The index of the row to extract (1-based index).
 * @return Matrix& A reference to a new matrix containing the extracted row.
 *         Note: The returned reference points to dynamically allocated memory.
 *         The caller is responsible for managing the memory to avoid leaks.
 */
Matrix& extract_row(Matrix& m, const int row){
	if (row <= 0 || row > m.n_row) {
        cout << "Matrix extract_row: invalid row number\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(m.n_column);  

    for (int i = 1; i <= m.n_column; i++) {
        (*m_aux)(i) = m(row, i); 
    }

    return *m_aux;
}

/**
 * @brief Extracts a specific column from a matrix and returns it as a new matrix.
 * 
 * @param m The input matrix from which the column will be extracted.
 * @param column The column number to extract (1-based index).
 * @return A reference to a new matrix containing the extracted column.
 *         Note: The returned matrix is dynamically allocated and must be managed properly to avoid memory leaks.
 */
Matrix& extract_column(Matrix& m, const int column){
	if (column <= 0 || column > m.n_column) {
        cout << "Matrix extract_column: invalid column number\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(m.n_row);  

    for (int i = 1; i <= m.n_row; i++) {
        (*m_aux)(i) = m(i,column); 
    }

    return *m_aux;
}

/**
 * @brief Assigns a column of a matrix with values from a vector.
 * 
 * @param m The matrix to which the column will be assigned.
 * @param v The vector containing the values to assign to the column.
 * @param column The index of the column to be assigned (1-based index).
 * @return Matrix& A reference to the modified matrix.
 */
Matrix& assign_column(Matrix& m, Matrix& v, const int column){

	if (m.n_row != v.n_column) {
		cout << "Matrix assign_column: invalid row number\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux = &m;
	for (int i = 1; i <= m.n_row; i++) {
        (*m_aux)(i, column) = v(i);  
    }

    return *m_aux;
	 
}

/**
 * @brief Assigns a row of a matrix with the values from another matrix (vector).
 * 
 * @param m The matrix whose row will be modified.
 * @param v The matrix (vector) containing the values to assign to the row.
 * @param row The row index in the matrix `m` to be assigned (1-based index).
 * @return Matrix& A reference to the modified matrix `m`.
 */
Matrix& assign_row(Matrix& m, Matrix& v, const int row){
	if (m.n_column!= v.n_column) {
		cout << "Matrix assign_row: invalid row number\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux = &m;
	for (int j = 1; j <= m.n_column; j++) {
        (*m_aux)(row, j) = v(j);  
    }

    return *m_aux;
}
