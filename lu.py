import multiprocessing
import numpy as np
from functools import partial

# 定义常量
T = 5
maxsize = 3000
maxrow = 3000
numBasis = 90000

# 初始化全局矩阵
gRows = np.zeros((maxrow, maxsize), dtype=int)
gBasis = np.zeros((numBasis, maxsize), dtype=int)

# 从文件中读取消元子和被消元行
def readBasis():
    iToBasis = {}
    with open("消元子.txt", "r") as file:
        for line in file:
            if not line.strip():
                continue
            indices = list(map(int, line.split()))
            row = indices[0]
            for pos in indices[1:]:
                index, offset = divmod(pos, 32)
                gBasis[row, index] |= (1 << offset)
            iToBasis[row] = gBasis[row]
    return iToBasis

def readRowsFrom(pos):
    iToFirst = {}
    with open("被消元行.txt", "r") as file:
        lines = file.readlines()[pos:pos+maxsize]
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        indices = list(map(int, line.split()))
        for tmp in indices:
            index, offset = divmod(tmp, 32)
            gRows[i, index] |= (1 << offset)
        iToFirst[i] = indices[0]
    return iToFirst

# 更新首项位置
def update(row):
    flag = False
    for i in range(maxsize-1, -1, -1):
        if gRows[row, i] != 0:
            flag = True
            pos = i * 32
            for k in range(31, -1, -1):
                if gRows[row, i] & (1 << k):
                    offset = k
                    break
            newfirst = pos + offset
            return newfirst if flag else -1
    return -1

# 核函数，执行并行消去操作
def kernel_task(i, gRows, gBasis, iToFirst, iToBasis):
    while iToFirst.get(i, -1) != -1:
        firstPos = iToFirst[i]
        if firstPos in iToBasis:
            for j in range(maxsize):
                gRows[i, j] ^= iToBasis[firstPos][j]
            iToFirst[i] = update(i)
        else:
            for j in range(maxsize):
                gBasis[firstPos, j] = gRows[i, j]
            iToBasis[firstPos] = gBasis[firstPos]
            del iToFirst[i]
    return gRows[i], iToFirst

def parallel_gaussian_elimination():
    currentBatchStart = 0
    while True:
        iToFirst = readRowsFrom(currentBatchStart)
        if not iToFirst:
            break

        with multiprocessing.Pool(processes=4) as pool:
            results = pool.map(partial(kernel_task, gRows=gRows, gBasis=gBasis, iToFirst=iToFirst, iToBasis=iToBasis), range(len(iToFirst)))

        for i, (row, first) in enumerate(results):
            gRows[i] = row
            if first != -1:
                iToFirst[i] = first

        currentBatchStart += maxsize
        if len(iToFirst) < maxsize:
            break

if __name__ == "__main__":
    iToBasis = readBasis()
    parallel_gaussian_elimination()
