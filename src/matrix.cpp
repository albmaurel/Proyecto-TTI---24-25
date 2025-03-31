#include "..\include\matrix.h"

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

double& Matrix::operator () (const int row, const int column) {

    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
        cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
    }

    return this->data[row - 1][column - 1];
}


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

Matrix& Matrix::operator = (Matrix& matrix2)
{
	if (this == &matrix2){ 
        return *this;
	}	
	if (this->n_row != matrix2.n_row || this->n_column != matrix2.n_column) {
        cout << "Matrix assign: error in n_row/n_column4\n";
        exit(EXIT_FAILURE);
    }
	
    for (int i = 1; i <= this->n_row; i++)
        for (int j = 1; j <= this->n_column; j++)
            (*this)(i,j)= matrix2(i,j);

    return *this;
}

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

Matrix& Matrix::operator / (Matrix& matrix2) {
    if (this->n_column != matrix2.n_row) {
        cout << "Matrix division: incompatible matrix dimensions for multiplication\n";
        exit(EXIT_FAILURE);
    }

    Matrix invMatrix2 = inverse(matrix2);

    return *this * invMatrix2;  
}

Matrix& Matrix::operator + (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + scalar;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator - (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - scalar;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator * (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) * scalar;
        }
    }

    return *m_aux;
}

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


ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& eye(int size) {
    Matrix *m_aux = new Matrix(size, size);  

    for (int i = 1; i <= size; i++) {
        for (int j = 1; j <= size; j++) {
            (*m_aux)(i, j) = (i == j) ? 1 : 0;
        }
    }

    return *m_aux;
}

Matrix& transpose(Matrix &m) {
    Matrix *m_aux = new Matrix(m.n_column, m.n_row); 

    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(j, i) = m(i, j);  
        }
    }

    return *m_aux;  
}

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