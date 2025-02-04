#include "pch.h"
using namespace std;

// パラメータ設定用の名前空間
namespace params {
    // 2元連立1階常微分方程式の係数
    const double a11 = 1.0;     // dx1/dt = a11*x1 + a12*x2
    const double a12 = 2.0;     
    const double a21 = 2.0;  // dx2/dt = a21*x1 + a22*x2
    const double a22 = 1.0;
    
    // 計算条件
    const double t_end = 5.0;   // 終了時刻
    const double dt = 0.01;     // 時間刻み幅
    
    // 初期条件
    const double x1_0 = 1.0;    // x1の初期値
    const double x2_0 = 0.0;    // x2の初期値
}

// 微分方程式の右辺を定義
Vector func(const Vector& x, double t) {
    Vector f(2);
    // 2元連立1階常微分方程式
    f(1) = params::a11 * x(1) + params::a12 * x(2);
    f(2) = params::a21 * x(1) + params::a22 * x(2);
    return f;
}

int main() {
    // 初期条件の設定
    Vector x(2);
    x(1) = params::x1_0;
    x(2) = params::x2_0;
    
    // 時間設定
    double t = 0.0;
    
    // 結果をファイルに出力
    ofstream data_file("rk_data.txt");
    ofstream plot_file("rk_plot.gp");
    
    digits(data_file, 8);  // ファイル出力の精度設定
    
    // データファイルにヘッダーとして計算条件を出力
    data_file << "# 2元連立1階常微分方程式の数値解" << endl;
    data_file << "# dx1/dt = " << params::a11 << "*x1 + " << params::a12 << "*x2" << endl;
    data_file << "# dx2/dt = " << params::a21 << "*x1 + " << params::a22 << "*x2" << endl;
    data_file << "# 時刻 t       x1          x2" << endl;
    
    // 初期値を出力
    data_file << t << " " << x(1) << " " << x(2) << endl;
    
    // メインの計算ループ
    while(t < params::t_end) {
        rk(x, func, t, params::dt);
        data_file << t << " " << x(1) << " " << x(2) << endl;
    }
    data_file.close();
    
    // gnuplotスクリプトの生成
    plot_file << "set grid" << endl;
    plot_file << "set title '2元連立1階常微分方程式の数値解'" << endl;
    plot_file << "set xlabel '時刻 t'" << endl;
    plot_file << "set ylabel '解 x1, x2'" << endl;
    plot_file << "plot 'rk_data.txt' using 1:2 title 'x1(t)' with lines, \\" << endl;
    plot_file << "     'rk_data.txt' using 1:3 title 'x2(t)' with lines" << endl;
    plot_file << "pause -1 'Press any key to exit'" << endl;
    plot_file.close();
    
    // gnuplotを実行
    cout << "数値計算が完了しました。" << endl;
    cout << "結果は rk_data.txt に保存されました。" << endl;
    cout << "gnuplotでグラフを表示します..." << endl;
    system("gnuplot rk_plot.gp");
    
    return 0;
}