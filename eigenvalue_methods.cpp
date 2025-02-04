#include "eigenvalue_methods.h"
using namespace std;

// 計算パラメータ
namespace params {
    const double eps = 1e-10;     // 収束判定値
    const int max_iter = 100;     // 最大反復回数を100回に減らす
}

// 単位行列の生成
Matrix create_identity(int n) {
    Matrix I(n);
    for (int i = 1; i <= n; ++i) {
        I(i, i) = 1.0;
    }
    return I;
}

// 行列のノルムを計算
double matrix_norm(const Matrix& A) {
    double sum = 0.0;
    for (int i = 1; i <= A.row(); ++i) {
        for (int j = 1; j <= A.col(); ++j) {
            sum += A(i, j) * A(i, j);
        }
    }
    return sqrt(sum);
}

// べき乗法による最大固有値計算
void power_method(const Matrix& A, double& eigenval, Vector& eigenvec) {
    int n = A.row();
    Vector x(n), x_new(n);
    
    // 初期ベクトルの設定(全ての要素を1に)
    for(int i = 1; i <= n; i++) x(i) = 1.0;
    normalize(x);
    
    cout << "反復計算開始" << endl;
    for(int iter = 0; iter < params::max_iter; iter++) {
        x_new = A * x;
        double lambda = norm(x_new);
        x_new = x_new / lambda;
        
        // 収束判定と途中経過の出力
        if(iter % 10 == 0) {
            cout << iter << "回目: 固有値 = " << lambda << endl;
        }
        
        if(norm(x_new - x) < params::eps) {
            eigenval = lambda;
            eigenvec = x_new;
            cout << "収束しました(" << iter + 1 << "回の反復)" << endl;
            return;
        }
        x = x_new;
    }
    cout << "警告: 最大反復回数に達しました" << endl;
    // 最後の推定値を返す
    eigenval = norm(A * x);
    eigenvec = x;
}

// 逆べき乗法による特定の固有値計算
void inverse_power_method(const Matrix& A, double shift, double& eigenval, Vector& eigenvec) {
    int n = A.row();
    Vector x(n), x_new(n);
    Matrix A_shifted = A;
    
    // シフト行列の作成 (A - σI)
    for(int i = 1; i <= n; i++) {
        A_shifted(i,i) -= shift;
    }
    
    // 初期ベクトルの設定(全ての要素を1に)
    for(int i = 1; i <= n; i++) x(i) = 1.0;
    normalize(x);
    
    // LU分解の準備
    int* p = new int[n+1];
    LUdcp(A_shifted, p);
    
    double prev_rayleigh = 0.0;
    cout << "反復計算開始" << endl;
    for(int iter = 0; iter < params::max_iter; iter++) {
        // (A - σI)y = x を解く
        x_new = x;
        LUslv(A_shifted, x_new, p);
        
        // 新しいベクトルを正規化
        double lambda = norm(x_new);
        x_new = x_new / lambda;
        
        // レイリー商を計算して固有値を推定
        Vector Ax = A * x_new;
        double rayleigh = x_new * Ax;  // 内積演算子*を使用
        
        // 収束判定と途中経過の出力
        if(iter % 10 == 0) {
            cout << iter << "回目: 固有値 ≈ " << rayleigh << endl;
        }
        
        // レイリー商の変化が小さければ収束とみなす
        if(abs(rayleigh - prev_rayleigh) < params::eps) {
            eigenval = rayleigh;
            eigenvec = x_new;
            cout << "収束しました(" << iter + 1 << "回の反復)" << endl;
            delete[] p;
            return;
        }
        
        prev_rayleigh = rayleigh;
        x = x_new;
    }
    
    // 最大反復回数に達した場合は最後の推定値を返す
    Vector Ax = A * x;
    eigenval = x * Ax;
    eigenvec = x;
    delete[] p;
    cout << "警告: 最大反復回数に達しました" << endl;
}

// Wilkinsonシフトの計算
double wilkinson_shift(const Matrix& H, int n) {
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
pair<complex<double>, complex<double>> eigenvalues_2x2(const Matrix& H, int i) {
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
void qr_decomposition(Matrix& A, Matrix& Q, Matrix& R) {
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
vector<complex<double>> eigenvalues_double_qr(Matrix A, int max_iterations, double tolerance) {
    int n = A.row();
    vector<complex<double>> eigenvalues;
    Matrix H = A;
    
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
            pair<complex<double>, complex<double>> vals = eigenvalues_2x2(H, 1);
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
        
        Matrix Q(current_size), R(current_size);
        Matrix H_shifted = H;
        for (int i = 1; i <= current_size; ++i) {
            H_shifted(i, i) -= shift;
        }
        
        qr_decomposition(H_shifted, Q, R);
        
        Matrix H_new = R * Q;
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

// 統合インターフェース
vector<complex<double>> compute_eigenvalues(const Matrix& A, const string& method, double shift) {
    vector<complex<double>> eigenvalues;
    
    if (method == "qr") {
        // ダブルQR法(全ての固有値を計算)
        return eigenvalues_double_qr(A);
    }
    else if (method == "power") {
        // べき乗法(最大固有値のみ)
        double eigenval;
        Vector eigenvec;
        power_method(A, eigenval, eigenvec);
        eigenvalues.push_back(complex<double>(eigenval, 0.0));
    }
    else if (method == "inverse") {
        // 逆べき乗法(シフト値に最も近い固有値)
        double eigenval;
        Vector eigenvec;
        inverse_power_method(A, shift, eigenval, eigenvec);
        eigenvalues.push_back(complex<double>(eigenval, 0.0));
    }
    else {
        cout << "エラー: 未知の計算方法です" << endl;
    }
    
    return eigenvalues;
}