#include <cmath>
#include <iostream>
#include <vector>
#include <complex> // std::complex を使用する場合

using namespace std;

// 既存のライブラリ (Vector_lib, Matrix_lib) の型定義と関数をここに含める (簡略化のため、必要な部分のみ)

// Vector_lib.h から必要な部分をコピー (簡略化)
class Vector {
public:
    Vector() : size_(3), maxsize_(3), p_(new double[3]), p1_(p_ - 1) { initialize(); } // デフォルトサイズを3に設定
    explicit Vector(int n) : size_(n), maxsize_(n), p_(new double[n]), p1_(p_ - 1) { initialize(); }
    Vector(const Vector& X) : size_(X.size_), maxsize_(X.size_), p_(new double[X.size_]), p1_(p_ - 1) {
        for (int i = 1; i <= size_; ++i) (*this)(i) = X(i);
    }
    ~Vector() { delete[] p_; }

    int size() const { return size_; }
    double& operator()(int i) const { return *(p1_ + i); }
    Vector& operator=(const Vector& X) {
        if (&X == this) return *this;
        resize(X.size_);
        for (int i = 1; i <= size_; ++i) (*this)(i) = X(i);
        return *this;
    }
    Vector operator-(const Vector& X) const {
        Vector result(size_);
        for (int i = 1; i <= size_; ++i) result(i) = (*this)(i) - X(i);
        return result;
    }
    Vector operator+(const Vector& X) const {
        Vector result(size_);
        for (int i = 1; i <= size_; ++i) result(i) = (*this)(i) + X(i);
        return result;
    }
    double operator*(const Vector& X) const {
        double sum = 0.0;
        for (int i = 1; i <= size_; ++i) sum += (*this)(i) * X(i);
        return sum;
    }
    Vector operator*(double c) const {
        Vector result(size_);
        for (int i = 1; i <= size_; ++i) result(i) = (*this)(i) * c;
        return result;
    }
    Vector operator/(double c) const {
        Vector result(size_);
        for (int i = 1; i <= size_; ++i) result(i) = (*this)(i) / c;
        return result;
    }
    friend ostream& operator<<(ostream& os, const Vector& X); // フレンド関数として宣言


    void resize(int n) {
        if (n > maxsize_) {
            delete[] p_;
            maxsize_ = n;
            p_ = new double[maxsize_];
            p1_ = p_ - 1;
        }
        size_ = n;
        initialize();
    }

private:
    double* p_;
    double* p1_;
    int size_;
    int maxsize_ = 0;
    void initialize() { for (int i = 1; i <= size_; ++i) (*this)(i) = 0.0; }
};

ostream& operator<<(ostream& os, const Vector& X) { // 実装を追加
    for (int i = 1; i <= X.size(); ++i) {
        os << X(i) << endl;
    }
    return os;
}


double norm(const Vector& X) {
    double w = 0.0;
    for (int i = 1; i <= X.size(); ++i) w += X(i) * X(i);
    return sqrt(w);
}

// Matrix_lib.h から必要な部分をコピー (簡略化)
class Matrix {
public:
    Matrix() : row_(3), col_(3), maxsize_(9), p_(new double[9]), p1_(p_ - 4) { initialize(); } // デフォルトサイズを3x3に設定
    explicit Matrix(int n) : row_(n), col_(n), maxsize_(n * n), p_(new double[n * n]), p1_(p_ - (n + 1)) { initialize(); }
    Matrix(int m, int n) : row_(m), col_(n), maxsize_(m * n), p_(new double[m * n]), p1_(p_ - (n + 1)) { initialize(); }
    Matrix(const Matrix& A) : row_(A.row_), col_(A.col_), maxsize_(A.maxsize_), p_(new double[A.maxsize_]), p1_(p_ - (A.col_ + 1)) {
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) (*this)(i, j) = A(i, j);
        }
    }
    ~Matrix() { delete[] p_; }

    int row() const { return row_; }
    int col() const { return col_; }
    double* operator[](int i) const { return p_ + i * col_; }
    double& operator()(int i, int j) const { return *(p1_ + i * col_ + j); }
    Matrix& operator=(const Matrix& A) {
        if (&A == this) return *this;
        resize(A.row_, A.col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) (*this)(i, j) = A(i, j);
        }
        return *this;
    }
    Matrix operator*(const Matrix& B) const {
        Matrix result(row_, B.col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= B.col(); ++j) {
                double sum = 0.0;
                for (int k = 1; k <= col_; ++k) sum += (*this)(i, k) * B(k, j);
                result(i, j) = sum;
            }
        }
        return result;
    }
    Matrix operator*(double c) const {
        Matrix result(row_, col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) result(i, j) = (*this)(i, j) * c;
        }
        return result;
    }
    Matrix operator+(const Matrix& B) const {
        Matrix result(row_, col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) result(i, j) = (*this)(i, j) + B(i, j);
        }
        return result;
    }
    Matrix operator-(const Matrix& B) const {
        Matrix result(row_, col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) result(i, j) = (*this)(i, j) - B(i, j);
        }
        return result;
    }
    Matrix operator-() const {
        Matrix result(row_, col_);
        for (int i = 1; i <= row_; ++i) {
            for (int j = 1; j <= col_; ++j) result(i, j) = -(*this)(i, j);
        }
        return result;
    }

    static Matrix identity(int n) {
        Matrix I(n);
        for (int i = 1; i <= n; ++i) I(i, i) = 1.0;
        return I;
    }

    friend ostream& operator<<(ostream& os, const Matrix& A); // フレンド関数として宣言


    void resize(int m, int n) {
        if (m * n > maxsize_) {
            delete[] p_;
            maxsize_ = m * n;
            p_ = new double[maxsize_];
            p1_ = p_ - (n + 1);
        }
        row_ = m;
        col_ = n;
        initialize();
    }


private:
    int row_, col_, maxsize_;
    double* p_;
    double* p1_;
    void initialize() { for (int i = 1; i <= row_; ++i) for (int j = 1; j <= col_; ++j) (*this)(i, j) = 0.0; }
};

ostream& operator<<(ostream& os, const Matrix& A) { // 実装を追加
    for (int i = 1; i <= A.row(); ++i) {
        for (int j = 1; j <= A.col(); ++j) {
            os << A(i, j) << " ";
        }
        os << endl;
    }
    return os;
}


Matrix trans(const Matrix& A) {
    Matrix w(A.col(), A.row());
    for (int i = 1; i <= A.row(); ++i) {
        for (int j = 1; j <= A.col(); ++j) w(j, i) = A(i, j);
    }
    return w;
}

Matrix tensor2(const Vector& x, const Vector& y) {
    Matrix w(x.size(), y.size());
    for (int i = 1; i <= x.size(); ++i)
        for (int j = 1; j <= y.size(); ++j) w(i, j) = x(i) * y(j);
    return w;
}


// Householder transformation (実装は以前と同じ)
void householder(Matrix& A) {
    int n = A.row();
    Vector v(n);
    Vector x(n);
    Vector u(n);
    Matrix P;
    Matrix H = Matrix(n);
    Matrix I = Matrix::identity(n);

    for (int k = 1; k <= n - 2; ++k) {
        for (int i = 1; i <= n; ++i) x(i) = 0.0;
        for (int i = k + 1; i <= n; ++i) x(i) = A(i, k);

        double norm_x = norm(x);
        if (norm_x < 1e-15) continue;

        if (x(k + 1) > 0) {
            norm_x = -norm_x;
        }
        v = x;
        v(k + 1) -= norm_x;

        double norm_v_sq = v * v;
        if (norm_v_sq < 1e-15) continue;

        u = v / sqrt(norm_v_sq);
        P = I - 2.0 * tensor2(u, u);

        H = trans(P) * A * P;
        for (int i = k; i <= n; ++i) {
            for (int j = k; j <= n; ++j) {
                A(i, j) = H(i, j);
            }
        }
        for (int i = k + 1; i <= n; ++i) A(k, i) = A(i, k) = (i == k + 1) ? H(k, i) : 0.0;
    }
}

// QR decomposition (実装は以前と同じ)
void qr_decomposition(Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.row();
    Q = Matrix(n);
    R = Matrix(n);
    Matrix tempA = A;
    Matrix tempQ = Matrix::identity(n);

    for (int j = 1; j <= n - 1; ++j) {
        for (int i = j + 1; i <= n; ++i) {
            if (abs(tempA(i, j)) > 1e-15) {
                double r = sqrt(tempA(j, j) * tempA(j, j) + tempA(i, j) * tempA(i, j));
                double c = tempA(j, j) / r;
                double s = -tempA(i, j) / r;

                Matrix G = Matrix::identity(n);
                G(j, j) = c;
                G(i, i) = c;
                G(j, i) = s;
                G(i, j) = -s;

                tempA = G * tempA;
                tempQ = tempQ * trans(G);
            }
        }
    }
    R = tempA;
    Q = tempQ;
}

// eigenvalues_qr 関数 (実装は以前と同じ)
Vector eigenvalues_qr(Matrix A, int max_iterations = 1000, double tolerance = 1e-10) {
    int n = A.row();
    Vector eigenvalues(n);
    Matrix Q, R;

    for (int iter = 0; iter < max_iterations; ++iter) {
        qr_decomposition(A, Q, R);
        A = R * Q;

        double off_diagonal_norm_sq = 0;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                if (i != j) {
                    off_diagonal_norm_sq += A(i, j) * A(i, j);
                }
            }
        }

        if (sqrt(off_diagonal_norm_sq) < tolerance) {
            break;
        }
    }

    for (int i = 1; i <= n; ++i) {
        eigenvalues(i) = A(i, i);
    }
    return eigenvalues;
}


// ダブルQR法による固有値計算関数 (簡略版 - シフトなし)
Vector eigenvalues_double_qr(Matrix A, int max_iterations = 1000, double tolerance = 1e-10) {
    int n = A.row();
    Vector eigenvalues(n);
    Matrix H = A;

    householder(H); // ヘッセンベルグ行列への変換 (対称行列でも適用可能)

    for (int iter = 0; iter < max_iterations; ++iter) {
        Matrix Q, R;
        qr_decomposition(H, Q, R); // QR分解
        H = R * Q; // RQで更新 (ダブルシフトは省略)


        double off_diagonal_norm_sq = 0;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                if (i != j) {
                    off_diagonal_norm_sq += H(i, j) * H(i, j);
                }
            }
        }


        if (sqrt(off_diagonal_norm_sq) < tolerance) {
            break; // 収束判定 (オフダイアゴナル要素が十分小さい)
        }
    }

    // 固有値を抽出 (対角成分から - 複素固有値ペアの処理は省略)
    for (int i = 1; i <= n; ++i) {
        eigenvalues(i) = H(i, i); // 実部のみ (複素固有値は近似的な実部となる可能性あり)
    }

    return eigenvalues;
}


int main() {
    // stdsize(3); // stdsize は削除

    Matrix A(3);
    A = 1, 2, 3,
        4, 5, 6,
        7, 8, 9; // 非対称行列の例

    cout << "元の行列 A:" << endl;
    cout << A << endl;

    Vector eigenvalues = eigenvalues_double_qr(A); // ダブルQR法を使用 (簡易版)
    cout << "固有値 (ダブルQR法 - 簡易版):" << endl;
    cout << eigenvalues << endl;

    // stop(); // stop() はコメントアウトまたは削除。exit(0) や return 0 で代替可能
    return 0;
}