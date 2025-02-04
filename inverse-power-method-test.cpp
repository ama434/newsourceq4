#include "../pch.h"
#include "eigenvalue_methods.h"
using namespace std;

int main() {
    cout << "テストケース1: 対称行列の固有値計算" << endl;
    Matrix A(3,3);
    A = 4.0, -2.0, 0.0,
       -2.0,  4.0, -2.0,
        0.0, -2.0,  4.0;
    
    cout << "行列 A:" << endl << A << endl;
    
    // 各手法での固有値計算
    cout << "\n1. ダブルQR法による全固有値:" << endl;
    vector<complex<double>> qr_eigenvals = compute_eigenvalues(A, "qr");
    for(const auto& val : qr_eigenvals) {
        if(abs(val.imag()) < 1e-10) {
            cout << val.real() << endl;
        } else {
            cout << val.real() << " + " << val.imag() << "i" << endl;
        }
    }
    
    cout << "\n2. べき乗法による最大固有値:" << endl;
    vector<complex<double>> power_eigenvals = compute_eigenvalues(A, "power");
    cout << power_eigenvals[0].real() << endl;
    
    cout << "\n3. 逆べき乗法による固有値 (シフト値 = 2.0):" << endl;
    vector<complex<double>> inverse_eigenvals = compute_eigenvalues(A, "inverse", 2.0);
    cout << inverse_eigenvals[0].real() << endl;
    
    cout << "\nテストケース2: 複素固有値を持つ非対称行列" << endl;
    Matrix B(3,3);
    B = 1.0, -1.0, 2.0,
        2.0,  3.0, -1.0,
        1.0,  1.0, 2.0;
    
    cout << "行列 B:" << endl << B << endl;
    vector<complex<double>> eigenvals_B = compute_eigenvalues(B, "qr");
    cout << "固有値:" << endl;
    for(const auto& val : eigenvals_B) {
        if(abs(val.imag()) < 1e-10) {
            cout << val.real() << endl;
        } else {
            cout << val.real() << " + " << val.imag() << "i" << endl;
        }
    }
    
    cout << "\nテストケース3: 純虚数固有値を持つ行列" << endl;
    Matrix C(2,2);
    C = 0.0, -1.0,
        1.0,  0.0;
    
    cout << "行列 C:" << endl << C << endl;
    vector<complex<double>> eigenvals_C = compute_eigenvalues(C, "qr");
    cout << "固有値:" << endl;
    for(const auto& val : eigenvals_C) {
        if(abs(val.imag()) < 1e-10) {
            cout << val.real() << endl;
        } else {
            cout << val.real() << " + " << val.imag() << "i" << endl;
        }
    }
    
    return 0;
}