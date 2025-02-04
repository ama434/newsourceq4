#include "pch.h"
using namespace std;

// 計算パラメータ
namespace params {
    const double eps = 1e-10;     // 収束判定値
    const int max_iter = 1000;    // 最大反復回数
}

void power_method(const Matrix& A, double& eigenval, Vector& eigenvec) {
    int n = A.row();
    Vector x(n), x_new(n);
    
    // 初期ベクトルの設定（全ての要素を1に）
    for(int i = 1; i <= n; i++) x(i) = 1.0;
    normalize(x);
    
    cout << "反復計算開始" << endl;
    for(int iter = 0; iter < params::max_iter; iter++) {
        x_new = A * x;
        double lambda = norm(x_new);
        x_new = x_new / lambda;
        
        // 収束判定と途中経過の出力
        if(iter % 10 == 0) {  // 10回ごとに経過を表示
            cout << iter << "回目: 固有値 = " << lambda << endl;
        }
        
        if(norm(x_new - x) < params::eps) {
            eigenval = lambda;
            eigenvec = x_new;
            cout << "収束しました（" << iter + 1 << "回の反復）" << endl;
            return;
        }
        x = x_new;
    }
    cout << "警告: 最大反復回数に達しました" << endl;
}

int main() {
    // テスト行列の定義
    Matrix A(3,3);
    A = 2, -1, 0,
       -1,  2, -1,
        0, -1,  2;
    
    // 最大固有値と対応する固有ベクトルを格納する変数
    double max_eigenval;
    Vector max_eigenvec;
    
    cout << "べき乗法による最大固有値・対応する固有ベクトルの計算：" << endl;
    cout << "入力行列 A:" << endl << A << endl;
    
    // べき乗法による計算
    power_method(A, max_eigenval, max_eigenvec);
    
    cout << endl << "計算結果：" << endl;
    cout << "最大固有値 = " << max_eigenval << endl;
    cout << "対応する固有ベクトル：" << endl << max_eigenvec << endl;
    
    // 結果の検証
    cout << endl << "結果の検証：" << endl;
    Vector Av = A * max_eigenvec;
    Vector lambda_v = max_eigenvec * max_eigenval;
    cout << "A * v =" << endl << Av << endl;
    cout << "λ * v =" << endl << lambda_v << endl;
    
    return 0;
}