
使用する行列は PPE3594

ここでは右辺ベクトルを変更した場合の性能評価をおこなう
以下の3種類の場合をファイルに書き込む
* all-one vector (引数が 0)
* srand(1)
* srand(2)

Fault-case/
    では 不適切な seed を選択した場合の結果を示している
    Reference-BLASでは NaN が
    OpenBLASでは NaN が出ずに収束しなかった