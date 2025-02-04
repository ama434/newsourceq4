#include "pch.h"
using namespace std;

// 計算パラメータ
namespace params {
    const double eps = 1e-10;  // 収束判定値
}

void jacobi_method(Matrix& A, Vector& eigenvals, Matrix& eigenvecs) {
    int n = A.row();
    eigenvecs.resize(n, n);
    eigenvals.resize(n);
    
    // 初期化：固有ベクトル行列を単位行列に
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            eigenvecs(i,j) = (i == j) ? 1.0 : 0.0;
        }
    }
    
    while(true) {
        // 最大非対角要素の探索
        double max_val = 0.0;
        int p = 1, q = 1;
        for(int i = 1; i <= n; i++) {
            for(int j = i+1; j <= n; j++) {
                if(abs(A(i,j)) > max_val) {
                    max_val = abs(A(i,j));
                    p = i;
                    q = j;
                }
            }
        }
        
        if(max_val < params::eps) break;  // 収束判定
        
        // 回転角の計算
        double theta;
        if(abs(A(p,p) - A(q,q)) < params::eps) {
            theta = pi/4.0;
        } else {
            theta = 0.5 * atan(2.0 * A(p,q) / (A(p,p) - A(q,q)));
        }
        
        // 回転行列による相似変換
        Matrix R(n,n);
        for(int i = 1; i <= n; i++) R(i,i) = 1.0;
        R(p,p) = cos(theta);
        R(p,q) = -sin(theta);
        R(q,p) = sin(theta);
        R(q,q) = cos(theta);
        
        A = trans(R) * A * R;
        eigenvecs = eigenvecs * R;
    }
    
    // 固有値を取り出し
    for(int i = 1; i <= n; i++) {
        eigenvals(i) = A(i,i);
    }
}

int main() {
    // テスト行列の定義
    Matrix A(7, 7);
    A = -8, 2, 32, 3, -9, -2, 23,
      0, 3, 52, -43, 57, -63, -2,
      0, 0, -9, -5, 33, 45, -3,
      0, 0, 0, -4, -3, 22, 6,
      0, 0, 0, 0, -6, 8, -123, 
      0, 0, 0, 0, 0, 6, 8, 
      0, 0, 0, 0, 0, 0, -13;


    // 固有値・固有ベクトルを格納する変数
    Vector eigenvals;
    Matrix eigenvecs;
    
    cout << "ヤコビ法による固有値・固有ベクトルの計算：" << endl;
    cout << "入力行列 A:" << endl << A << endl;
    
    // ヤコビ法による計算
    jacobi_method(A, eigenvals, eigenvecs);
    
    cout << "固有値：" << endl << eigenvals << endl;
    cout << "固有ベクトル（列ベクトルとして格納）：" << endl << eigenvecs << endl;
    
    // 結果の検証
    cout << "検証：A * v_i = λ_i * v_i" << endl;
    for(int i = 1; i <= A.row(); i++) {
        Vector v(A.row());
        for(int j = 1; j <= A.row(); j++) {
            v(j) = eigenvecs(j,i);
        }
        Vector Av = A * v;
        Vector lambdav = v * eigenvals(i);
        cout << i << "番目の固有ベクトルの検証：" << endl;
        cout << "A * v =" << endl << Av << endl;
        cout << "λ * v =" << endl << lambdav << endl << endl;
    }
    
    return 0;
}