#include "pch.h"
#include <iomanip>
using namespace std;

class MatrixInputHelper {
private:
    Matrix& A;
    int currentRow;
    int currentCol;
    
    void clearInputBuffer() {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    
    bool getValidSize(int& n) {
        cout << "行列のサイズを入力してください (n > 0): ";
        if (!(cin >> n) || n <= 0) {
            cout << "無効なサイズです。正の整数を入力してください。" << endl;
            clearInputBuffer();
            return false;
        }
        return true;
    }
    
    bool getValidElement(double& element) {
        if (!(cin >> element)) {
            cout << "無効な入力です。数値を入力してください。" << endl;
            clearInputBuffer();
            return false;
        }
        return true;
    }

    string formatNumber(double num) {
        ostringstream oss;
        if (abs(num - round(num)) < 1e-10) {
            // 整数の場合
            oss << setw(8) << fixed << setprecision(0) << num;
        } else {
            // 小数の場合は3桁まで表示
            oss << setw(8) << fixed << setprecision(3) << num;
        }
        return oss.str();
    }

public:
    MatrixInputHelper(Matrix& matrix) : A(matrix), currentRow(1), currentCol(1) {}
    
    bool inputSize() {
        int n;
        if (!getValidSize(n)) return false;
        
        stdsize(n);
        A.resize(n);
        currentRow = 1;
        currentCol = 1;
        return true;
    }
    
    bool inputElement() {
        cout << "A(" << currentRow << "," << currentCol << ") = ";
        double element;
        if (!getValidElement(element)) return false;
        
        A(currentRow, currentCol) = element;
        
        currentCol++;
        if (currentCol > A.col()) {
            currentCol = 1;
            currentRow++;
        }
        return true;
    }
    
    bool isComplete() {
        return currentRow > A.row();
    }
    
    void displayProgress() {
        #ifdef _WIN32
            system("color"); // Windowsでの色表示有効化
        #endif
        
        cout << "\n現在の行列:\n" << endl;
        for (int i = 1; i <= A.row(); ++i) {
            cout << "    ";
            for (int j = 1; j <= A.col(); ++j) {
                if (i == currentRow && j == currentCol) {
                    // 現在の入力位置をハイライト表示（背景色：青、文字色：白）
                    cout << "\033[44;37m" << setw(8) << "?" << "\033[0m";
                } else if (i < currentRow || (i == currentRow && j < currentCol)) {
                    // 入力済みの要素
                    cout << formatNumber(A(i,j));
                } else {
                    // 未入力の要素
                    cout << "\033[90m" << setw(8) << "?" << "\033[0m";
                }
            }
            cout << endl;
        }
    }
    
    void goBack() {
        currentCol--;
        if (currentCol < 1) {
            currentRow--;
            if (currentRow < 1) {
                currentRow = 1;
                currentCol = 1;
            } else {
                currentCol = A.col();
            }
        }
    }
};

void displayEigenvalues(const Vector& zr, const Vector& zi) {
    cout << "\n固有値:" << endl;
    cout << setprecision(14) << fixed;
    
    for(int i = 1; i <= zr.size(); ++i) {
        if(abs(zi(i)) < 1e-10) {
            cout << "λ" << i << " = " << setw(20) << zr(i) << endl;
        } else {
            cout << "λ" << i << " = " << setw(20) << zr(i);
            if(zi(i) >= 0) {
                cout << " + " << setw(20) << zi(i) << "i" << endl;
            } else {
                cout << " - " << setw(20) << abs(zi(i)) << "i" << endl;
            }
        }
    }
}

int main() {
    Matrix A;
    Vector zr, zi;
    MatrixInputHelper inputHelper(A);
    
    cout << "固有値計算プログラム" << endl;
    cout << "===================" << endl;
    
    while (!inputHelper.inputSize()) {
        cout << "再試行しますか？ (y/n): ";
        char response;
        cin >> response;
        if (response != 'y' && response != 'Y') return 0;
        cin.ignore();
    }
    
    while (!inputHelper.isComplete()) {
        inputHelper.displayProgress();
        cout << "\n[Enter]：次の要素を入力" << endl;
        cout << "[b]：一つ前の要素に戻る" << endl;
        cout << "[q]：プログラムを終了" << endl;
        cout << "選択してください: ";
        
        string input;
        getline(cin, input);
        
        if (input == "b" || input == "B") {
            inputHelper.goBack();
        }
        else if (input == "q" || input == "Q") {
            return 0;
        }
        else {
            if (!inputHelper.inputElement()) {
                cout << "入力をやり直してください。" << endl;
            }
        }
    }
    
    zr.resize(A.row());
    zi.resize(A.row());
    eig(A, zr, zi);
    displayEigenvalues(zr, zi);
    
    cout << "\nEnterキーを押して終了...";
    cin.ignore();
    cin.get();
    
    return 0;
}