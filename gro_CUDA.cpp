#include <cuda.h>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

using namespace std;

const int T = 5;
const int maxsize = 3000;
const int maxrow = 3000;
const int numBasis = 90000;

map<int, int *> iToBasis;
map<int, int> iToFirst;
map<int, int *> ans;

fstream RowFile("被消元行.txt", ios::in | ios::out);
fstream BasisFile("消元子.txt", ios::in | ios::out);

int gRows[maxrow][maxsize];
int gBasis[numBasis][maxsize];

// CUDA设备函数，用于更新首项位置
__device__ void updateFirstPos(int *d_gRows, int *d_iToFirst, int row, int maxsize) {
    int newFirstPos = -1;
    for (int i = maxsize - 1; i >= 0; i--) {
        if (d_gRows[row * maxsize + i] != 0) {
            for (int bit = 31; bit >= 0; bit--) {
                if (d_gRows[row * maxsize + i] & (1 << bit)) {
                    newFirstPos = i * 32 + bit;
                    break;
                }
            }
            if (newFirstPos != -1) break;
        }
    }
    d_iToFirst[row] = newFirstPos;
}

// CUDA核函数，执行并行消去操作
__global__ void kernelGaussianElimination(int *d_gRows, int *d_gBasis, int *d_iToFirst, int numRows, int maxsize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < numRows) {
        while (d_iToFirst[idx] != -1) {
            int firstPos = d_iToFirst[idx];
            if (firstPos != -1) {
                for (int j = 0; j < maxsize; j++) {
                    d_gRows[idx * maxsize + j] ^= d_gBasis[firstPos * maxsize + j];
                }
                // 更新首项位置
                updateFirstPos(d_gRows, d_iToFirst, idx, maxsize);
            }
        }
    }
}


void reset() {
    memset(gRows, 0, sizeof(gRows));
    memset(gBasis, 0, sizeof(gBasis));

    RowFile.close();
    BasisFile.close();

    RowFile.open("被消元行.txt", ios::in | ios::out);
    BasisFile.open("消元子.txt", ios::in | ios::out);

    iToBasis.clear();
    iToFirst.clear();
    ans.clear();
}

void readBasis() {
    for (int i = 0; i < maxrow; i++) {
        if (BasisFile.eof()) {
            return;
        }
        string tmp;
        bool flag = false;
        int row = 0;
        getline(BasisFile, tmp);
        stringstream s(tmp);
        int pos;
        while (s >> pos) {
            if (!flag) {
                row = pos;
                flag = true;
                iToBasis.insert(pair<int, int *>(row, gBasis[row]));
            }
            int index = pos / 32;
            int offset = pos % 32;
            gBasis[row][index] |= (1 << offset);
        }
        flag = false;
        row = 0;
    }
}

int readRowsFrom(int pos) {
    iToFirst.clear();
    if (RowFile.is_open())
        RowFile.close();
    RowFile.open("被消元行.txt", ios::in | ios::out);
    memset(gRows, 0, sizeof(gRows));
    string line;
    for (int i = 0; i < pos; i++) {
        getline(RowFile, line);
    }
    for (int i = pos; i < pos + maxsize; i++) {
        int tmp;
        getline(RowFile, line);
        if (line.empty()) {
            cout << "End of File!" << endl;
            return i;
        }
        bool flag = false;
        stringstream s(line);
        while (s >> tmp) {
            if (!flag) {
                iToFirst.insert(pair<int, int>(i - pos, tmp));
            }
            int index = tmp / 32;
            int offset = tmp % 32;
            gRows[i - pos][index] |= (1 << offset);
            flag = true;
        }
    }
    return -1;
}

void update(int row) {
    bool flag = false;
    for (int i = maxsize - 1; i >= 0; i--) {
        if (gRows[row][i] == 0)
            continue;
        else {
            if (!flag)
                flag = true;
            int pos = i * 32;
            int offset = 0;
            for (int k = 31; k >= 0; k--) {
                if (gRows[row][i] & (1 << k)) {
                    offset = k;
                    break;
                }
            }
            int newfirst = pos + offset;
            iToFirst.erase(row);
            iToFirst.insert(pair<int, int>(row, newfirst));
            break;
        }
    }
    if (!flag)
        iToFirst.erase(row);
    return;
}

void writeResult(ofstream &out) {
    for (auto it = ans.rbegin(); it != ans.rend(); it++) {
        int *result = it->second;
        int max = it->first / 32 + 1;
        for (int i = max; i >= 0; i--) {
            if (result[i] == 0)
                continue;
            int pos = i * 32;
            for (int k = 31; k >= 0; k--) {
                if (result[i] & (1 << k)) {
                    out << k + pos << " ";
                }
            }
        }
        out << endl;
    }
}

void cuda_parallel() {
    int *d_gRows, *d_gBasis, *d_iToFirst;
    size_t rows_size = maxrow * maxsize * sizeof(int);
    size_t basis_size = numBasis * maxsize * sizeof(int);
    size_t iToFirst_size = maxrow * sizeof(int);

    cudaMalloc(&d_gRows, rows_size);
    cudaMalloc(&d_gBasis, basis_size);
    cudaMalloc(&d_iToFirst, iToFirst_size);

    cudaMemcpy(d_gRows, gRows, rows_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gBasis, gBasis, basis_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_iToFirst, &iToFirst[0], iToFirst_size, cudaMemcpyHostToDevice);

    int threadsPerBlock = 1024;
    int numBlocks = (maxrow + threadsPerBlock - 1) / threadsPerBlock;

    kernelGaussianElimination<<<numBlocks, threadsPerBlock>>>(d_gRows, d_gBasis, d_iToFirst, maxrow, maxsize);

    cudaMemcpy(gRows, d_gRows, rows_size, cudaMemcpyDeviceToHost);

    cudaFree(d_gRows);
    cudaFree(d_gBasis);
    cudaFree(d_iToFirst);
}

int main() {

    reset();
    readBasis();
    int currentBatchStart = 0;
    int endFlag;

    while (true) {
        endFlag = readRowsFrom(currentBatchStart);
        if (endFlag == -1) {
            currentBatchStart += maxsize;
        } else {
            break;
        }
        cuda_parallel();
    }

    ofstream outFile("结果.txt");
    writeResult(outFile);
    outFile.close();
    return 0;
}
