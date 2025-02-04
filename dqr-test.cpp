#include "../pch.h"
#include <complex>
#include <vector>

using namespace std;

// 単位行列の生成
matrix_ create_identity(int n) {
    matrix_ I(n);
    for (int i = 1; i <= n; ++i) {
        I(i, i) = 1.0;
    }
    return I;
}

// 行列のノルムを計算
double matrix_norm(const matrix_& A) {
    double sum = 0.0;
    for (int i = 1; i <= A.row(); ++i) {
        for (int j = 1; j <= A.col(); ++j) {
            sum += A(i, j) * A(i, j);
        }
    }
    return sqrt(sum);
}

// Wilkinsonシフトの計算
double wilkinson_shift(const matrix_& H, int n) {
    if (abs(H(n, n-1)) < 1e-14) return H(n, n);
    
    double a = H(n-1, n-1);
    double b = H(n-1, n);
    double c = H(n, n-1);
    double d = H(n, n);
    
    double tr = a + d;
    double det = a * d - b * c;
    double disc = tr * tr - 4.0 * det;
    
    if (disc < 0) return tr / 2.0;
    
    double sqrt_disc = sqrt(disc);
    double lambda1 = (tr + sqrt_disc) / 2.0;
    double lambda2 = (tr - sqrt_disc) / 2.0;
    
    return abs(lambda1 - H(n,n)) < abs(lambda2 - H(n,n)) ? lambda1 : lambda2;
}

// 2x2ブロックの固有値計算
pair<complex<double>, complex<double> > eigenvalues_2x2(const matrix_& H, int i) {
    double a = H(i, i);
    double b = H(i, i+1);
    double c = H(i+1, i);
    double d = H(i+1, i+1);
    
    double tr = a + d;
    double det = a * d - b * c;
    double disc = tr * tr - 4.0 * det;
    
    if (disc >= 0) {
        double sqrt_disc = sqrt(disc);
        double lambda1 = (tr + sqrt_disc) / 2.0;
        double lambda2 = (tr - sqrt_disc) / 2.0;
        return make_pair(complex<double>(lambda1, 0), complex<double>(lambda2, 0));
    } else {
        double real = tr / 2.0;
        double imag = sqrt(-disc) / 2.0;
        return make_pair(complex<double>(real, imag), complex<double>(real, -imag));
    }
}

// QR分解
void qr_decomposition(matrix_& A, matrix_& Q, matrix_& R) {
    int n = A.row();
    Q = create_identity(n);
    R = A;
    
    for (int j = 1; j <= n - 1; ++j) {
        for (int i = j + 1; i <= n; ++i) {
            if (abs(R(i, j)) < 1e-14) continue;
            
            double a = R(j, j);
            double b = R(i, j);
            double r = hypot(a, b);
            if (r < 1e-14) continue;
            
            double c = a / r;
            double s = -b / r;
            
            for (int k = j; k <= n; ++k) {
                double temp = R(j, k);
                R(j, k) = c * temp - s * R(i, k);
                R(i, k) = s * temp + c * R(i, k);
            }
            
            for (int k = 1; k <= n; ++k) {
                double temp = Q(k, j);
                Q(k, j) = c * temp - s * Q(k, i);
                Q(k, i) = s * temp + c * Q(k, i);
            }
        }
    }
    
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j < i; ++j) {
            if (abs(R(i, j)) < 1e-13) R(i, j) = 0.0;
        }
    }
    
    Q = trans(Q);
}

// ダブルQR法による固有値計算
vector<complex<double> > eigenvalues_double_qr(matrix_ A, int max_iterations = 200, double tolerance = 1e-12) {
    int n = A.row();
    vector<complex<double> > eigenvalues;
    matrix_ H = A;
    
    hess(H);
    
    double initial_norm = matrix_norm(H);
    if (initial_norm < 1e-14) {
        for (int i = 1; i <= n; ++i) {
            eigenvalues.push_back(complex<double>(0.0, 0.0));
        }
        return eigenvalues;
    }
    
    tolerance *= initial_norm;
    
    int current_size = n;
    int iteration_count = 0;
    
    while (current_size > 1 && iteration_count < max_iterations) {
        if (current_size == 2) {
            pair<complex<double>, complex<double> > vals = eigenvalues_2x2(H, 1);
            eigenvalues.push_back(vals.first);
            eigenvalues.push_back(vals.second);
            break;
        }
        
        if (abs(H(current_size, current_size-1)) < tolerance) {
            eigenvalues.push_back(complex<double>(H(current_size, current_size), 0));
            current_size--;
            continue;
        }
        
        double shift = wilkinson_shift(H, current_size);
        
        matrix_ Q(current_size), R(current_size);
        matrix_ H_shifted = H;
        for (int i = 1; i <= current_size; ++i) {
            H_shifted(i, i) -= shift;
        }
        
        qr_decomposition(H_shifted, Q, R);
        
        matrix_ H_new = R * Q;
        for (int i = 1; i <= current_size; ++i) {
            H_new(i, i) += shift;
        }
        
        for (int i = 1; i <= current_size; ++i) {
            for (int j = 1; j <= current_size; ++j) {
                H(i, j) = H_new(i, j);
                if (abs(H(i, j)) < tolerance) H(i, j) = 0.0;
            }
        }
        
        iteration_count++;
    }
    
    if (current_size == 1) {
        eigenvalues.push_back(complex<double>(H(1, 1), 0));
    }
    
    return eigenvalues;
}

int main() {
    cout << "テストケース1: 実数固有値を持つ対称行列" << endl;
    matrix_ A(3);
    A = 4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0;

    cout << "行列 A:" << endl << A << endl;
    vector<complex<double> > eigenvalues1 = eigenvalues_double_qr(A);
    cout << "固有値:" << endl;
    for (size_t i = 0; i < eigenvalues1.size(); ++i) {
        const complex<double>& lambda = eigenvalues1[i];
        if (abs(lambda.imag()) < 1e-10) {
            cout << lambda.real() << endl;
        } else {
            cout << lambda.real() << " + " << lambda.imag() << "i" << endl;
        }
    }

    cout << "\nテストケース2: 複素固有値を持つ非対称行列" << endl;
    matrix_ B(3);
    B = 1.0, -1.0, 2.0,
        2.0, 3.0, -1.0,
        1.0, 1.0, 2.0;

    cout << "行列 B:" << endl << B << endl;
    vector<complex<double> > eigenvalues2 = eigenvalues_double_qr(B);
    cout << "固有値:" << endl;
    for (size_t i = 0; i < eigenvalues2.size(); ++i) {
        const complex<double>& lambda = eigenvalues2[i];
        if (abs(lambda.imag()) < 1e-10) {
            cout << lambda.real() << endl;
        } else {
            cout << lambda.real() << " + " << lambda.imag() << "i" << endl;
        }
    }

    cout << "\nテストケース3: 純虚数固有値を持つ行列" << endl;
    matrix_ C(2);
    C = 0.0, -1.0,
        1.0, 0.0;

    cout << "行列 C:" << endl << C << endl;
    vector<complex<double> > eigenvalues3 = eigenvalues_double_qr(C);
    cout << "固有値:" << endl;
    for (size_t i = 0; i < eigenvalues3.size(); ++i) {
        const complex<double>& lambda = eigenvalues3[i];
        if (abs(lambda.imag()) < 1e-10) {
            cout << lambda.real() << endl;
        } else {
            cout << lambda.real() << " + " << lambda.imag() << "i" << endl;
        }
    }

    return 0;
}