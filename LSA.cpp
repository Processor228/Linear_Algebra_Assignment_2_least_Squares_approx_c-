#include <bits/stdc++.h>

using namespace std;

struct DimensionalException : public exception{
    const char * what() {
        return "Error: the dimensional problem occurred\n";
    }
};

template <typename T>
class ColumnVector {
protected:
    vector<T> entries;
public:
    explicit ColumnVector (int size) {
        for (int row = 0; row < size; row++) {
            entries.push_back(0);
        }
    }

    friend istream & operator>>(istream& in, ColumnVector & v) {
        for (T & entry : v.entries) {
            in >> entry;
        }
        return in;
    }

    friend ostream & operator<<(ostream& out, ColumnVector v) {
        for (T entry : v.entries) {
            out << entry << '\n';  ///StupidConversion(entry) << '\n';
        }
        return out;
    }

    int dim() {
        return this->entries.size();
    }

    template<typename K>
    ColumnVector<T> operator+(ColumnVector<K> & add) {
        ColumnVector<T> res(add.dim());
        for (int i = 0 ; this->entries.size(); i++) {
            res[i] = add.entries[i] + this->entries[i];
        }
        return res;
    }

    template<typename K>
    ColumnVector<T> operator-(ColumnVector<K> & add) {
        ColumnVector<T> res(add.dim());
        for (int i = 0 ; this->entries.size(); i++) {
            res[i] = this->entries[i] - add.entries[i];
        }
        return res;
    }

    template<typename K>
    ColumnVector<K> operator=(ColumnVector<K> v) {
        for (int i = 0 ; i < v.dim(); i++) {
            this->entries[i] = v.entries[i];
        }
        return *this;
    }

    double mag() {
        double res = 0;
        for (int i = 0; i < this->entries.size(); i++) {
            res += this->entries[i] * this->entries[i];
        }
        return sqrt(res);
    }

    T& operator[](int entry) {
        return entries[entry];
    }

    void discardAlmostZeros() {
        for (int k = 0; k < this->dim(); k++) {
            this->entries[k] = discardSmallerNumber(this->entries[k]);
        }
    }
};

template<typename Type>
class Matrix {
protected:
    vector<vector<Type>> entries;
public:

    friend istream & operator>>(istream& in, Matrix<Type>& m) {
        vector<unsigned int> shapeOfMatrix = m.shape();
        for (int row = 0; row < shapeOfMatrix[0]; row++) {
            for (int col = 0; col < shapeOfMatrix[1]; col++) {
                in >> m[row][col];
            }
        }
        return in;
    }

    friend ostream & operator<<(ostream& out, Matrix<Type>& m) {
        vector<unsigned int> shapeOfMatrix = m.shape();
        for (int row = 0; row < shapeOfMatrix[0]; row++) {
            for (int col = 0; col < shapeOfMatrix[1]; col++) {
                if (col == shapeOfMatrix[1]-1){
                    out << m[row][col];   //StupidConversion(m[row][col]);
                    continue;
                }
                out << m[row][col] << " ";  // StupidConversion(m[row][col]) << " ";
            }
            cout << "\n";
        }
        return out;
    }

    explicit Matrix (vector<vector<Type>>& numbers) {
        entries = numbers;
    }

    ~Matrix() = default;

    Matrix (int rows, int columns) {
        for (int row = 0; row < rows; row++) {
            vector<Type> newRow;
//            newRow.assign(columns, 0);
            for (int i = 0 ; i < columns; i++) {
                newRow.push_back(0);
            }
//            newRow.resize(8);
            entries.push_back(newRow);
        }
    }

    vector<unsigned int> shape() {
        vector<unsigned int> res = {static_cast<unsigned int>(entries.size()), static_cast<unsigned int>(entries[0].size())};
        return res;
    }

    template<typename K>
    ColumnVector<K> operator*(ColumnVector<K> v) {
        ColumnVector<K> res(this->shape()[0]);
        for (int i = 0 ; i < this->shape()[0]; i++) {
            K dotProd = 0;
            for (int j = 0 ; j < v.dim(); j++) {
                dotProd += v[j] * this->entries[i][j];
            }
            res[i] = dotProd;
        }
        return res;
    };

    template<typename T>
    Matrix<T>& operator=(Matrix<T>* a) {
        vector<unsigned int> shapeOfThis = shape();
        if (a->shape()[0] != shapeOfThis[0] || a->shape()[1] != shapeOfThis[1]) {
            throw DimensionalException();
        }

        for (int i = 0; i < shapeOfThis[0]; i++) {
            for (int j = 0; j < shapeOfThis[1]; j++) {
                entries[i][j] = a->entries[i][j];
            }
        }
        return *this;
    }

    template<class T>
    auto operator+(Matrix<T>& a) {   /// tested, works good
        vector<unsigned int> shapeOfAdded = a.shape();
        vector<unsigned int> shapeOfThis = shape();
        for (int i = 0; i < shapeOfAdded.size(); i++) {
            if (shapeOfAdded[i] != shapeOfThis[i]) {
                throw DimensionalException();
            }
        }

        auto res = new Matrix<T>(shapeOfThis[0], shapeOfThis[1]);

        for (int i = 0; i < shapeOfThis[0]; i++) {
            for (int j = 0; j < shapeOfThis[1]; j++) {
                res -> entries[i][j] = entries[i][j] + a.entries[i][j];
            }
        }
        return res;
    }

    auto operator-(Matrix & a) {   /// tested, works fine
        vector<unsigned int> shapeOfAdded = a.shape();
        vector<unsigned int> shapeOfThis = shape();
        for (int i = 0; i < shapeOfAdded.size(); i++) {
            if (shapeOfAdded[i] != shapeOfThis[i]) {
                throw DimensionalException();
            }
        }

        auto* res = new Matrix(shapeOfThis[0], shapeOfThis[1]);

        for (int i = 0; i < shapeOfThis[0]; i++) {
            for (int j = 0; j < shapeOfThis[1]; j++) {
                res->entries[i][j] = entries[i][j] - a.entries[i][j];
            }
        }
        return res;
    }

    template<typename t>
    auto operator*(Matrix<t>& a) {
        vector<unsigned int> shapeOfMultiplier = a.shape();
        vector<unsigned int> shapeOfThis = shape();

        if(shapeOfThis[1] != shapeOfMultiplier[0]) {
            throw DimensionalException();
        }

        auto res = new Matrix<t>(shapeOfThis[0], shapeOfMultiplier[1]);

        for (int colOfMultiplier = 0; colOfMultiplier < shapeOfMultiplier[1]; colOfMultiplier++) {
            for (int rowOfMultiplied = 0; rowOfMultiplied < shapeOfThis[0]; rowOfMultiplied++) {
                t newEl = 0;
                for (int colOfMultiplied = 0; colOfMultiplied < shapeOfMultiplier[0]; colOfMultiplied++) {
                    newEl +=t(entries[rowOfMultiplied][colOfMultiplied]) * t(a[colOfMultiplied][colOfMultiplier]);
                }
                res->operator[](rowOfMultiplied)[colOfMultiplier] = newEl;
            }
        }
        return res;
    }

    template<typename t>
    Matrix<t> operator*(Matrix<t>* a) {
        vector<unsigned int> shapeOfMultiplier = a->shape();
        vector<unsigned int> shapeOfThis = shape();

        if(shapeOfThis[1] != shapeOfMultiplier[0]) {
            throw DimensionalException();
        }

        auto res = new Matrix<t>(shapeOfThis[0], shapeOfMultiplier[1]);

        for (int colOfMultiplier = 0; colOfMultiplier < shapeOfMultiplier[1]; colOfMultiplier++) {
            for (int rowOfMultiplied = 0; rowOfMultiplied < shapeOfThis[0]; rowOfMultiplied++) {
                Type newEl = 0;
                for (int colOfMultiplied = 0; colOfMultiplied < shapeOfMultiplier[0]; colOfMultiplied++) {
                    newEl += t(entries[rowOfMultiplied][colOfMultiplied]) * (a->operator[](colOfMultiplied)[colOfMultiplier]);
                }
                res->operator[](rowOfMultiplied)[colOfMultiplier] = newEl;
            }
        }
        return *res;
    }

    template<typename t>
    auto operator*(Matrix<t>* a) {
        vector<unsigned int> shapeOfMultiplier = a->shape();
        vector<unsigned int> shapeOfThis = shape();

        if(shapeOfThis[1] != shapeOfMultiplier[0]) {
            throw DimensionalException();
        }

        auto res = new Matrix<t>(shapeOfThis[0], shapeOfMultiplier[1]);

        for (int colOfMultiplier = 0; colOfMultiplier < shapeOfMultiplier[1]; colOfMultiplier++) {
            for (int rowOfMultiplied = 0; rowOfMultiplied < shapeOfThis[0]; rowOfMultiplied++) {
                Type newEl = 0;
                for (int colOfMultiplied = 0; colOfMultiplied < shapeOfMultiplier[0]; colOfMultiplied++) {
                    newEl += t(entries[rowOfMultiplied][colOfMultiplied]) * (a->operator[](colOfMultiplied)[colOfMultiplier]);
                }
                res->operator[](rowOfMultiplied)[colOfMultiplier] = newEl;
            }
        }
        return res;
    }

    auto T() {
        auto res = new Matrix(entries[0].size(), entries.size());
        for (int i = 0; i < entries.size(); i++) {
            for (int j = 0; j < entries[0].size(); j++) {
                res -> entries[j][i] = entries[i][j];
            }
        }
        return *res;
    }

    void assign (int row, int col, Type val) {
        entries[row][col] = val;
    }

    Type& getRef(int row, int col) {
        return entries[row][col];
    }

    vector<Type>& operator[](int row) {
        return this->entries[row];
    }

    void discardAlmostZeros() {
        for (int i = 0; i < this->shape()[0]; i++) {
            for (int k = 0; k < this->shape()[1]; k++) {
                this->entries[i][k] = discardSmallerNumber(this->entries[i][k]);
            }
        }
    }

    int pivotPermutationNeeds (int pivot) {
        int max = abs(this->entries[pivot][pivot]);
        int bestPivot = pivot;
        for (int i = pivot; i < this->entries.size(); i++) {
            if (abs(this->entries[i][pivot]) > max) {
                max = abs(this->entries[i][pivot]);
                bestPivot = i;
            }
        }
        return bestPivot;
    }
};

template<typename T>
class SquareMatrix : public Matrix<T> {
public:
    explicit SquareMatrix(int n) : Matrix<T>(n, n) {

    }
    template<typename K>
    explicit SquareMatrix<T>(Matrix<K>* m) : Matrix<T>(m->shape()[0], m->shape()[0]) {
        for (int row = 0 ; row < m->shape()[0]; row++) {
            for (int col = 0 ; col < m->shape()[1]; col ++) {
                this->entries[row][col] = m->getRef(row, col);
            }
        }
    }

    template<typename K>
    explicit SquareMatrix<T>(Matrix<K> m) : Matrix<T>(m.shape()[0], m.shape()[0]) {
        for (int row = 0 ; row < m.shape()[0]; row++) {
            for (int col = 0 ; col < m.shape()[1]; col ++) {
                this->entries[row][col] = m.getRef(row, col);
            }
        }
    }

/**
 * @param pivot pivots already done
 * @return int row to permute current pivot with
 */
    int pivotPermutationNeeds (int pivot) {
        int max = abs(this->entries[pivot][pivot]);
        int bestPivot = pivot;
        for (int i = pivot; i < this->entries.size(); i++) {
            if (abs(this->entries[i][pivot]) > max) {
                max = abs(this->entries[i][pivot]);
                bestPivot = i;
            }
        }
        return bestPivot;
    }
};

template <typename T>
class EliminationMatrix : public SquareMatrix<T> {
public:
    explicit EliminationMatrix(int n) : SquareMatrix<T>(n) {
        for (int i = 0; i < n; i++) {
            this->entries[i][i] = 1;
        }
    }
    template<typename type>  // n is for size
    explicit EliminationMatrix(int n, int row, int col, Matrix<type>& toElim) : SquareMatrix<T>(n) {
        for (int i = 0; i < n; i++) {
            this->entries[i][i] = 1;
        }
        this->entries[row][col] = -(double)toElim[row][col] / toElim[col][col];
    }
};

//template <typename T>
class PermutationMatrix : public SquareMatrix<int> {
public:
    explicit PermutationMatrix(int n, int row1, int row2) : SquareMatrix<int>(n) {
        this->entries[row2][row1] = 1;
        this->entries[row1][row2] = 1;
        for (int i = 0; i < n; i++) {
            if (i != row1 && i != row2) {
                this->entries[i][i] = 1;
            }
        }
    }
};

class IdentityMatrix : public EliminationMatrix<int> {

public:
    explicit IdentityMatrix(int n) : EliminationMatrix<int>(n) {

    }
};


SquareMatrix<double> Inverse (SquareMatrix<double>& A) {
    int size = A.shape()[0];
    int step = 0;

    SquareMatrix<double> Aug(size);
    Aug = SquareMatrix<double>(IdentityMatrix(size));
    for (int pivot = 0; pivot < size; pivot++) {
        int pivotPermutation = A.pivotPermutationNeeds(pivot);
        PermutationMatrix p(size, pivot, pivotPermutation);
        A = SquareMatrix<double>(p * A);
        Aug = SquareMatrix<double>(p * Aug);
        for (int row = pivot + 1; row < size; row++) {
            if (abs(A[row][pivot]) < pow(10, -10)) {
                continue;
            }
            EliminationMatrix<double> elim (size, row, pivot, A);
            A = SquareMatrix<double>(elim * A);
            Aug = SquareMatrix<double>(elim * Aug);
        }
    }
    for (int pivot = A.shape()[0]-1; pivot >= 1; pivot--) {
        for (int row = pivot - 1; row >= 0; row--) {
            EliminationMatrix<double> elim (size, row, pivot, A);
            A = SquareMatrix<double>(elim * A);
            Aug = SquareMatrix<double>(elim * Aug);
        }
    }
    for (int pivot = 0 ; pivot < A.shape()[0]; pivot++) {
        double num = A[pivot][pivot];
        A[pivot][pivot] /= num;
        for (int col = 0; col < A.shape()[1]; col ++) {
            Aug[pivot][col] /= num;
        }
    }
    return Aug;
}

int main () {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(4);
    cout.setf(ios::fixed, ios::floatfield);

    int points;
    int degree;
    cin >> points;
    vector<pair<double, double>> coordinates;
    for (int i = 0; i < points; i++) {
        int x, y;
        cin >> x >> y;
        coordinates.emplace_back(x, y);
    }
    cin >> degree;

    //-----------(Performing actions with input data)------------------------------//
    Matrix<double> A(points, degree + 1);
    ColumnVector<double> b(points);

    for (int row = 0; row < points; row++) {
        for (int col = 0; col < degree + 1; col++) {  // forming the matrix A
            A[row][col] = pow(coordinates[row].first, col);
        }
    }
    for (int i = 0; i < points; i++) {
        b[i] = coordinates[i].second;     // forming the b vector
    }
	
    cout << "A:\n" << A;
    SquareMatrix<double>A_T_A (degree + 1);
    A_T_A = SquareMatrix<double>(A.T() * A);
    cout << "A_T*A:\n" << A_T_A;

    SquareMatrix<double>A_T_A_INV = Inverse(A_T_A);
    cout << "(A_T*A)^-1:\n" << A_T_A_INV;

    ColumnVector<double> A_Tb = A.T() * b;
    cout << "A_T*b:\n" << A_Tb;

    ColumnVector<double> x_hat = A_T_A_INV * A_Tb;
    cout << "x~:\n";
    cout << x_hat;

    return 0;
}

