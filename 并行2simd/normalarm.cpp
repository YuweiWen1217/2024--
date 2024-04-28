#include <iostream>
#include <chrono>
#include <arm_neon.h>
using namespace std;

const int N = 1024;
const int MAX = 1000;
float m[N][N];


void print()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << m[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

void initialize()
{
    // 遍历矩阵的每个元素，将其初始化为特定的值
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // 这里可以根据需要设置初始化值
            m[i][j] = (i * 2 + j) % 1000;
        }
    }
}

void Normal(float M[][N])
{
    float A[N][N];
    memcpy(A, M, sizeof(float) * N * N);
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
        {
            float u = A[i][k] / A[k][k];
            for (int j = k + 1; j < N; j++)
                A[i][j] -= u * A[k][j];
            A[i][k] = 0.0;
        }
}
void NEON(float M[][N]) {
    float A[N][N]__attribute__((aligned(64)));
    memcpy(A, M, sizeof(float) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            float div = A[j][i] / A[i][i];
            float32x4_t ratio = vmovq_n_f32(div);
            for (int k = i; k < N; k += 4) {
                float32x4_t v3 = vld1q_f32(&A[i][k]);  // 加载A[i][k]到向量
                float32x4_t v4 = vld1q_f32(&A[j][k]);  // 加载A[j][k]到向量
                // 存储结果回数组
                vst1q_f32(&A[j][k], vmlsq_f32(v4, ratio, v3));
            }
        }
    }
}
void NEON_aligned(float M[][N]) {
    float A[N][N] __attribute__((aligned(64)));
    memcpy(A, M, sizeof(float) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            float div = A[j][i] / A[i][i];
            // 计算比例
            float32x4_t ratio = vmovq_n_f32(div);
            for (int k = i / 4 * 4; k < N; k += 4) {
                float32x4_t v3 = vld1q_f32(&A[i][k]);
                float32x4_t v4 = vld1q_f32(&A[j][k]);
                vst1q_f32(&A[j][k], vmlsq_f32(v4, ratio, v3));
            }
        }
    }
}


int main()
{
    cout<<"2";
    auto total_duration = std::chrono::microseconds(0); // 总时间初始化为 0
    for (int i = 0; i < 1; ++i) {
        initialize();
        auto start = std::chrono::high_resolution_clock::now(); // 获取起始时间点
        NEON(m);
        // 在这里执行你要计时的代码
        auto end = std::chrono::high_resolution_clock::now();   // 获取结束时间点
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start); // 计算时间差
        total_duration += duration; // 累加每次执行所花费的时间
    }
    std::cout << "NEON: " << total_duration.count() << " microseconds" << std::endl;

    total_duration = std::chrono::microseconds(0); // 总时间初始化为 0
    for (int i = 0; i < 1; ++i) {
        initialize();
        auto start = std::chrono::high_resolution_clock::now(); // 获取起始时间点
        Normal(m);
        auto end = std::chrono::high_resolution_clock::now();   // 获取结束时间点
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start); // 计算时间差
        total_duration += duration; // 累加每次执行所花费的时间
    }
    std::cout << "Normal: " << total_duration.count() << " microseconds" << std::endl;

    total_duration = std::chrono::microseconds(0); // 总时间初始化为 0
    for (int i = 0; i < 1; ++i) {
        initialize();
        auto start = std::chrono::high_resolution_clock::now(); // 获取起始时间点
        NEON_aligned(m);
        auto end = std::chrono::high_resolution_clock::now();   // 获取结束时间点
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start); // 计算时间差
        total_duration += duration; // 累加每次执行所花费的时间
    }
    std::cout << "NEON_aligned: " << total_duration.count() << " microseconds" << std::endl;
    return 0;
}