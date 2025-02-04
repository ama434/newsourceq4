#ifndef _eigenvalue_methods_h
#define _eigenvalue_methods_h

#include "../pch.h"
#include <complex>
#include <vector>

// べき乗法による最大固有値計算
void power_method(const Matrix& A, double& eigenval, Vector& eigenvec);

// 逆べき乗法による特定の固有値計算
void inverse_power_method(const Matrix& A, double shift, double& eigenval, Vector& eigenvec);

// QR分解
void qr_decomposition(Matrix& A, Matrix& Q, Matrix& R);

// Wilkinsonシフトの計算
double wilkinson_shift(const Matrix& H, int n);

// 2x2ブロックの固有値計算
std::pair<std::complex<double>, std::complex<double>> eigenvalues_2x2(const Matrix& H, int i);

// ダブルQR法による固有値計算
std::vector<std::complex<double>> eigenvalues_double_qr(Matrix A, int max_iterations = 200, double tolerance = 1e-12);

// 統合インターフェース
std::vector<std::complex<double>> compute_eigenvalues(const Matrix& A, const std::string& method = "qr", double shift = 0.0);

#endif // _eigenvalue_methods_h