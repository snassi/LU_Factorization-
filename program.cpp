#include <iostream>
#include <string>
#include <fstream>
#include <exception>

using namespace std;

class ArrayException : public exception{};
class ArrayMemoryException : public ArrayException {};
class ArrayBoundsException : public ArrayException {};

template<class DataType>
class AbstractArrayClass {
public:
	virtual int size() const = NULL;
	virtual DataType& operator[] (int k) = NULL;
	friend ostream& operator << <DataType>(ostream& s, AbstractArrayClass<DataType>& ac) {
		s << "[";
		for (int i = 0; i < ac.size(); i++) {
			if (i > 0) {
				s << ",";
			}
			s << ac[i];
		}
		s << "]";
		return s;
	}
};


const int ARRAY_CLASS_DEFAULT_SIZE = 1;
template <class DataType>
class ArrayClass : virtual public AbstractArrayClass<DataType>
{

	DataType* paObject;
	int _size;
	void copy(const ArrayClass<DataType>& ac);
public:
	ArrayClass();
	ArrayClass(int n);
	//ArrayClass(int n, DataType* val);
	ArrayClass(const ArrayClass<DataType>& ac);
//	ArrayClass(int n, const DataType* value);
	ArrayClass(int n, const DataType& val);
	virtual ~ArrayClass();
	virtual int size() const;
	virtual DataType& operator[] (int k);
	void operator= (const ArrayClass<DataType>& ac);
	void insert(int n, const DataType& val);



};
template<class DataType>
void ArrayClass<DataType>::insert(int n, const DataType& val)
{
	if (n > _size) throw ArrayMemoryException();
	if (paObject == NULL) throw ArrayMemoryException();


	paObject[n] = val;
}


template<class DataType>
ArrayClass<DataType>::ArrayClass(int n) {
	_size = 0;
	paObject = new DataType[n];
	if (paObject == NULL) throw ArrayMemoryException();
	_size = n;

}
template<class DataType>
ArrayClass<DataType>::~ArrayClass() {

	if (paObject != NULL) delete[] paObject;
	paObject = NULL;
	_size = 0;
}
template<class DataType>
ArrayClass<DataType>::ArrayClass(int n, const DataType& val)
{
	_size = 0;
	paObject = new DataType[n];
	if (paObject == NULL) throw ArrayMemoryException();
	_size = n;
	for (int i = 0; i < n; i++)
		paObject[i] = val;
}
template<class DataType>
ArrayClass<DataType>::ArrayClass() {
	_size = 0;
	paObject = new DataType[ARRAY_CLASS_DEFAULT_SIZE];
	if (paObject == NULL) throw ArrayMemoryException();
	_size = ARRAY_CLASS_DEFAULT_SIZE;
}
template<class DataType>
ArrayClass<DataType>::ArrayClass(const ArrayClass<DataType>& ac) {
	copy(ac);
}
/*
template<class DataType>
ArrayClass<DataType>::ArrayClass(int n, const DataType* value) {
	_size = 0;
	paObject = new DataType[n];
	if (paObject == NULL) throw ArrayMemoryException();
	_size = n;
	for (int i = 0; i < n; i++)
		paObject[i] = value[i];
}
*/
template<class DataType>
void ArrayClass<DataType>::copy(const ArrayClass<DataType>& ac) {
	_size = 0;
	paObject = new DataType[ac._size];
	if (paObject == NULL) throw ArrayMemoryException();
	_size = ac._size;
	for (int i = 0; i < _size; i++) {
		paObject[i] = ac.paObject[i];
	}
}
template <class DataType>
void ArrayClass<DataType>::operator=
(const ArrayClass<DataType>& ac) {
	if (&ac != this) {
		if (paObject != NULL)delete[] paObject;
		copy(ac);
	}
}
template<class DataType>
int ArrayClass<DataType>::size() const
{
	return _size;
}
template<class DataType>
DataType& ArrayClass<DataType>::operator[](int k) {
	if ((k < 0) || (k >= size())) throw ArrayBoundsException();
	return paObject[k];
}

//THE MATRIX CLASS

class MatrixIncompatibleException : public ArrayException {};

template<class DataType>
class Matrix : public AbstractArrayClass<ArrayClass <DataType>>
{
protected:
	ArrayClass<ArrayClass<DataType>*>* theRows;
	void copy(Matrix<DataType>& m);
	void deleteRows();
public:
	Matrix();
	Matrix(int n, int m);
	Matrix(int n, int m, DataType v);
	Matrix(Matrix& m);
	virtual ~Matrix();
	void operator= (Matrix& m);
	void operator= (const DataType* list);
	virtual int size() const;
	int columns();
	int rows();
	virtual ArrayClass<DataType>& operator[] (int index);
};

template<class DataType>
Matrix<DataType>::Matrix()
{
	theRows = new ArrayClass<ArrayClass<DataType>*>(1, NULL);
	if (theRows == NULL) {
		throw ArrayMemoryException();
	}
	(*theRows)[0] = new ArrayClass<DataType>();
	if ((*theRows)[0] == NULL)
	{
		throw ArrayMemoryException();
	}
}

template<class DataType>
Matrix<DataType>::Matrix(int n, int m) {
	theRows = new ArrayClass<ArrayClass<DataType>*>(n, NULL);
	if (theRows == NULL) {
		throw ArrayMemoryException();
	}
	for (int i = 0; i < n; i++) {
		(*theRows)[i] = new ArrayClass<DataType>(m);
		if ((*theRows)[i] == NULL)
		{
			throw ArrayMemoryException();
		}
	}
}

template<class DataType>
Matrix<DataType>::Matrix(int n, int m, DataType v) {
	theRows = new ArrayClass<ArrayClass<DataType>*>(n, NULL);
	if (theRows == NULL) {
		throw ArrayMemoryException();
	}
	for (int i = 0; i < n; i++) {
		(*theRows)[i] = new ArrayClass<DataType>(m, v);
		if ((*theRows)[i] == NULL)
		{
			throw ArrayMemoryException();
		}
	}
}

template<class DataType>
void Matrix<DataType>::deleteRows() {
	if (theRows != NULL) {
		for (int i = 0; i < theRows->size(); i++) {
			if ((*theRows)[i] != NULL) delete (*theRows)[i];
			(*theRows)[i] = NULL;
		}
		delete theRows;
	}
	theRows = NULL;
}

template<class DataType>
Matrix<DataType>::~Matrix() {
	deleteRows();
}

template<class DataType>
void Matrix<DataType>::copy(Matrix<DataType>& m)
{
	deleteRows();
	theRows = new ArrayClass < ArrayClass<DataType>* >
		(m.size(), NULL);
	if (theRows == NULL) throw ArrayMemoryException();
	for (int i = 0; i < m.size(); i++) {
		(*theRows)[i] = new ArrayClass<DataType>(m[i]);
		if ((*theRows)[i] == NULL) throw ArrayMemoryException();
	}
}

template<class DataType>
Matrix<DataType>::Matrix(Matrix<DataType>& m) {
	theRows = NULL;
	copy(m);
}

template<class DataType>
void Matrix<DataType>::operator= (Matrix<DataType>& m) {
	copy(m);
}

template<class DataType>
void Matrix<DataType>::operator= (const DataType* list) {
	copy(list);
}

template<class DataType>
int Matrix<DataType>::rows() {
	return theRows->size();
}

template<class DataType>
int Matrix<DataType>::columns() {
	return (*this)[0].size();
}

template<class DataType>
int Matrix<DataType>::size() const
{
	return theRows->size();
}

template<class DataType>
ArrayClass<DataType>& Matrix<DataType>::operator[](int index) {
	return (*(*theRows)[index]);
}

Matrix<double> operator*(Matrix<double>& MatA, Matrix<double>& MatB) {
	if (MatA.columns() != MatB.rows())
		throw MatrixIncompatibleException();
	Matrix<double> MatC(MatA.rows(), MatB.columns(), 0);
	for (int i = 0; i < MatA.rows(); i++) 
		for (int j = 0; j < MatB.columns(); j++) {
			for (int k = 0; k < MatA.columns(); k++)
				MatC[i][j] = MatC[i][j] + 
				MatA[i][k] * MatB[k][j];
		}
	return MatC;
}
/*
template<class DataType>
class AbstractVector : virtual public AbstractArrayClass<DataType> {
public:
	virtual void insert(const DataType& item, int index) = NULL;
		// insert a new object at position index in the vector
	virtual void remove(int index) = NULL;
		// remove the object at postion index of the vector
	virtual void add(const DataType& item) = NULL;
		// adds item at end of a vector
};

template<class DataType>
class Vector :
	virtual public ArrayClass<DataType>,
	virtual public AbstractVector<DataType>
{
protected:
	int _currSize;
	int _incFactor;
public:
	Vector();
	Vector(int n);
};
*/
Matrix<double> upper(Matrix<double> A, Matrix<double> L, Matrix<double> U, int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];

		}
	}
	return U;
}

Matrix<double> lower(Matrix<double> A, Matrix<double> L, Matrix<double> U, int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];

		}
	}
	return L;
}

/*
void crout(double const **A, double **L, double **U, int n) {
int i, j, k;
double sum = 0;

for (i = 0; i < n; i++) {
U[i][i] = 1;
}

for (j = 0; j < n; j++) {
for (i = j; i < n; i++) {
sum = 0;
for (k = 0; k < j; k++) {
sum = sum + L[i][k] * U[k][j];
}
L[i][j] = A[i][j] - sum;
}

for (i = j; i < n; i++) {
sum = 0;
for(k = 0; k < j; k++) {
sum = sum + L[j][k] * U[k][i];
}
if (L[j][j] == 0) {
printf("det(L) close to 0!\n Can't divide by 0...\n");
exit(EXIT_FAILURE);
}
U[j][i] = (A[j][i] - sum) / L[j][j];
}
}
}
*/



class GraphException : public exception {};
class GraphDuplicateEdge : public GraphException{};
class GraphEdgeOutOfBounds : public GraphException {};
class GraphMemory : public GraphException {};
class GraphNegativeCount : public GraphException {};
class GraphNonExistentEdge : public GraphException {};
class GraphVertexOutOfBounds : public GraphException{};

template<class VertexObject, class EdgeObject>
class AbstractGraph
{
public:
	virtual ~AbstractGraph();
	virtual int vertexCount() = NULL;
		// returns the number of vertices
	virtual int edgeCount() = NULL;
		// returns the number of edges
	virtual bool hasEdge(int start, int end) = NULL;
		// return true if the edges (start,end) exists in the 
		// graph, otherwise false
	virtual VertexObject& vertexInfo(int v) = NULL;
		// return any data stored with vertex v
	virtual double edgeWeight(int start, int end) = NULL;
		// returns the weight of the edge (start,end) if the 
		// edge exists, otherwise returns 0;
		// in an unweughted graph, returns 1 if the edge exists
	virtual EdgeObject& edgeInfo(int start, int end) = NULL;
		// returns any data stored with the edge (start, end)
//	virtual Vector<int> neighbors (int v) = NULL;
		// returns the neighbors of vertex v
	virtual void setVertexInfo
		(int v, VertexObject& info) = NULL;
		// set the data stored with the vertex v
	virtual void setEdgeInfo
		(int v, int end, EdgeObject& info) = NULL;
		// sets the data stored with the edge (start, end)
	virtual void deleteEdge(int start, int end) = NULL;
		// deletes the edge (start,end)
};

int main() {

	cout << "Welcome to the Matrix Decomposer, this Program uses the Crout Matrix Decomposition Algorithm" << endl;
	cout << "Along with that this Program calculates the determinant of a Matriz that you imput" << endl;
	cout << "This gave me a little stress, so be thankful..." << endl;

	while (true) {

		cout << "Give Dimensions for a nXn matrix: " << endl;
		int a;
		int low;
		int up;
		double entry;
		double determinant;
		cin >> a;
		
		Matrix<double> L3(a, a, 0);
		Matrix<double> U3(a, a, 0);

		Matrix<double> A3(a, a, 0);

		for (int i = 0; i < a; i++) {

			for (int j = 0; j < a; j++) {
				cout << "Enter in for row: " << i+1 << " and column: " << j+1 << ": ";
				cin >> entry;
				A3[i][j] = entry;
			}
		}

		L3 = lower(A3, L3, U3, a);
		U3 = upper(A3, L3, U3, a);

		cout << "The Lower Triangular Matrix is: " << L3 << endl;
		cout << "The Upper Triangular Matrix is: " << U3  << endl;
		cout << "The product of the Lower and Upper Triangular Matrices is: " << L3 * U3 << endl;
		cout << "The Original Matrix is: " << A3 << endl;
		determinant = 1;
		for (int i = 0; i < a; i++) {
			determinant = determinant * L3[i][i];
		}
		cout << "The Determinant of the Matrix is: " << determinant << endl;


	}
}
