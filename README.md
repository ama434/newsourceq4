# 数値計算法テスト対策用コード

このリポジトリは数値計算法の試験範囲から、特に重要な以下のアルゴリズムの実装を提供します：

1. ルンゲ・クッタ法による微分方程式の解法
2. ヤコビ法による固有値・固有ベクトルの計算
3. べき乗法による最大固有値の計算

> [!WARNING]
> このソースコードは macOS 上で制作されたものです。そのため、**Windows での動作は保証されません**。
> Windows での動作を確認する場合は、適宜修正を行ってください。

## 必要条件

- C++ コンパイラ（g++推奨）
- GNU Make
- gnuplot（グラフ表示用）

### OS別の準備

#### macOS
```bash
# HomebrewでC++コンパイラとgnuplotをインストール
brew install gcc
brew install gnuplot
```

#### Ubuntu/Debian
```bash
# 必要なパッケージをインストール
sudo apt update
sudo apt install build-essential
sudo apt install gnuplot
```

#### Windows
1. MinGW または Cygwin で g++ と make をインストール
2. gnuplot をインストール（[公式サイト](http://www.gnuplot.info/)から）

## セットアップ手順

1. 配布されている `matrix` を任意のディレクトリに配置
```
/path/to/your/workspace/
└── matrix/                  # 元のマトリックスライブラリ
```

2. 本リポジトリを `matrix` ディレクトリ内に配置
```bash
cd matrix
git clone https://github.com/ama434/newsourceq4.git
```

3. ディレクトリ構成の確認
```
/path/to/your/workspace/
└── matrix/                  # 元のマトリックスライブラリ
    ├── AutoDiff_lib.cpp
    ├── AutoDiff_lib.h
    ├── Matrix_lib.cpp
    ├── Matrix_lib.h
    ├── Vector_lib.cpp
    ├── Vector_lib.h
    ├── samples/
    │   ├── eig.cpp
    │   └── ...
    └── newsourceq4/        # このリポジトリ
        ├── Makefile
        ├── rk-solver.cpp   # ルンゲ・クッタ法
        ├── jacobi-test.cpp # ヤコビ法
        └── power-method.cpp # べき乗法
```

## 使用方法

`newsourceq4` ディレクトリで作業します：
```bash
cd newsourceq4
```

### 1. ルンゲ・クッタ法
```bash
# コンパイルと実行
make MAIN_SRC=rk-solver.cpp
./matrix

# グラフの表示（プログラムの指示に従ってください）
gnuplot rk_plot.gp  # Unix系の場合
```

パラメータの変更は `rk-solver.cpp` 内の `params` 名前空間で行えます：
```cpp
namespace params {
    // 微分方程式の係数
    const double a11 = 0.0;
    const double a12 = 1.0;
    const double a21 = -400.0;
    const double a22 = -6.0;
    
    // その他のパラメータ
    ...
}
```

### 2. ヤコビ法
```bash
make MAIN_SRC=jacobi-test.cpp
./matrix
```

### 3. べき乗法
```bash
make MAIN_SRC=power-method.cpp
./matrix
```

### 参照用の main.cpp の使用
親ディレクトリの `main.cpp` を使用する場合：
```bash
make
./matrix
```

## テストでの使用方法

1. **微分方程式を解く場合**：
   - 運動方程式や共振回路の式を手計算で1階の連立微分方程式に変形
   - 係数を `rk-solver.cpp` の `params` に設定
   - 実行して結果を確認

2. **固有値問題を解く場合**：
   - 対称行列の場合は `jacobi-test.cpp` を使用
   - 最大固有値のみ必要な場合は `power-method.cpp` を使用

## トラブルシューティング

### Makefile のパス関連のエラー
- `newsourceq4` ディレクトリが `matrix` ディレクトリの直下にあることを確認
- 相対パス `..` が親ディレクトリの `matrix` を正しく指していることを確認

### gnuplot の実行エラー
- gnuplot が正しくインストールされていることを確認
- Windows の場合、gnuplot のパスが正しく設定されていることを確認

## ライセンス
このコードは教育目的で提供されています。