#include <iostream>
#include <chrono>
#include <ratio>
#include <ctime>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <cstring>
using namespace std;
const int N = 200;
const int MAX = 1000;
float A[N][N];
void initialize()
{
    srand(static_cast<unsigned>(time(0)));
    // 构造原始上三角矩阵:也即算法处理后期望得到的正确结果
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
            a[i][j] = 0.0;
        a[i][i] = 1.0;
        for (int j = i + 1; j < N; j++)
            a[i][j] = rand() % MAX;
    }

    // 后面每一行全部加上前面的行，可看作算法的逆向操作
    // 这样构造出的矩阵一定可以化为标准的上三角形式
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
            for (int j = 0; j < N; j++)
                a[i][j] = (int)(a[i][j] + a[k][j]) % MAX;
}
void Normal()
{
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
        {
            float u = a[i][k] / a[k][k];
            for (int j = k + 1; j < N; j++)
                // 正确做法应该是0 - n-1，但已知必然化为上三角矩阵，可在后续直接处理
                a[i][j] -= u * a[k][j];
            a[i][k] = 0.0; // 化为上三角矩阵，左下角清0
        }
}

void SSE_ALIGN()
{ // 对齐的SSE算法
    for (int k = 0; k < N; k++)
    {
        __m128 t1 = _mm_set1_ps(A[k][k]);
        int j = k + 1;

        // 对非对齐部分进行串行处理
        while ((long long)(&A[k][j]) % 16)
        {
            A[k][j] = A[k][j] / A[k][k];
            j++;
        }
        // cout << &m[k][j]<<endl;
        for (; j + 4 <= N; j += 4)
        {
            __m128 t2 = _mm_load_ps(&A[k][j]); // 已对齐，用load和store指令
            t2 = _mm_div_ps(t2, t1);
            _mm_store_ps(&A[k][j], t2);
        }
        // 对行末进行串行处理
        for (; j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            __m128 vik = _mm_set1_ps(A[i][k]);
            j = k + 1;
            // 下面同理
            while ((long long)(&A[k][j]) % 16)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
                j++;
            }
            for (; j + 4 <= N; j += 4)
            {
                __m128 vkj = _mm_load_ps(&A[k][j]);
                __m128 vij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_mul_ps(vik, vkj);
                vij = _mm_sub_ps(vij, vx);
                _mm_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void SSE()
{
    for (int k = 0; k < N; k++)
    {
        __m128 t1 = _mm_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 4 <= N; j += 4)
        {
            __m128 t2 = _mm_loadu_ps(&A[k][j]);
            t2 = _mm_div_ps(t2, t1);
            _mm_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            __m128 vik = _mm_set1_ps(A[i][k]);
            for (j = k + 1; j + 4 <= N; j += 4)
            {
                __m128 vkj = _mm_loadu_ps(&A[k][j]);
                __m128 vij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_mul_ps(vik, vkj);
                vij = _mm_sub_ps(vij, vx);
                _mm_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void AVX()
{
    for (int k = 0; k < N; k++)
    {
        __m256 t1 = _mm256_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 8 <= N; j += 8)
        {
            __m256 t2 = _mm256_loadu_ps(&A[k][j]);
            t2 = _mm256_div_ps(t2, t1);
            _mm256_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8)
            {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void AVX_ALIGN()
{
    for (int k = 0; k < N; k++)
    {
        __m256 t1 = _mm256_set1_ps(A[k][k]);
        int j = k + 1;
        while ((long long)(&A[k][j]) % 32)
        {
            A[k][j] = A[k][j] / A[k][k];
            j++;
        }
        for (; j + 8 <= N; j += 8)
        {
            __m256 t2 = _mm256_load_ps(&A[k][j]);
            t2 = _mm256_div_ps(t2, t1);
            _mm256_store_ps(&A[k][j], t2);
        }
        for (; j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            j = k + 1;
            while ((long long)(&A[k][j]) % 32)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
                j++;
            }
            for (; j + 8 <= N; j += 8)
            {
                __m256 vkj = _mm256_load_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}
int main()
{
    long long head, tail, freq;
    double total_time = 0;
    for (int i = 0; i < 100; i++)
    {
        initialize();
        QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        AVX_ALIGN()
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        total_time += (tail - head) * 1000.0 / freq;

    }
    std::cout << "time: " << total_time << " ms" << std::endl;
    return 0;
}
