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

double& Matrix::operator () (const int row, const int column) {

    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
        cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
    }

    return this->data[row - 1][column - 1];
}

double& Matrix::operator () (const int v_size) {

    if (v_size <= 0 || v_size > this->n_column*this->n_row) {
        cout << "Matrix get: error in column\n";
        exit(EXIT_FAILURE);
    }

return this->data[(v_size-1)/this->n_column][(v_size-1)%this->n_column];
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

    Matrix invMatrix2 = inv(matrix2);

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

Matrix& zeros(const int v_size) {
	Matrix *m_aux = new Matrix(v_size);
	
	for(int j = 1; j <= v_size; j++) {
		(*m_aux)(j) = 0;
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

double norm(Matrix& a){
	 double r=0;
	 for(int i=1;i<=a.n_column;i++){
		 r +=pow(a(i),2);
	 }
	 double solv=sqrt(r);
	 return solv;
}

double dot(Matrix& a,Matrix& b){
	 double r=0;
	 for(int i=1;i<=a.n_column;i++)
		 r +=a(i)*b(i);
	 return r;
}

Matrix& cross(Matrix& a,Matrix& b){
	Matrix *m_aux = new Matrix(3);
	(*m_aux)(1) = a(2) * b(3) - a(3) * b(2);  
    (*m_aux)(2) = a(3) * b(1) - a(1) * b(3);  
    (*m_aux)(3) = a(1) * b(2) - a(2) * b(1);  
    
    return (*m_aux);
}

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

Matrix& extract_column(Matrix& m, const int column){
	if (column <= 0 || column > m.n_column) {
        cout << "Matrix extract_row: invalid column number\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(m.n_row);  

    for (int i = 1; i <= m.n_row; i++) {
        (*m_aux)(i) = m(i,column); 
    }

    return *m_aux;
}

Matrix& assign_column(Matrix& m, Matrix& v, const int column){
	if (m.n_row != v.n_row) {
		cout << "Matrix assign_column: invalid row number\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux = &m;
	for (int i = 1; i <= m.n_row; i++) {
        (*m_aux)(i, column) = v(i);  
    }

    return *m_aux;
	 
}

Matrix& assign_row(Matrix& m, Matrix& v, const int row){
	if (m.n_column!= v.n_column) {
		cout << "Matrix assign_row: invalid column number\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux = &m;
	for (int j = 1; j <= m.n_column; j++) {
        (*m_aux)(row, j) = v(j);  
    }

    return *m_aux;
}
