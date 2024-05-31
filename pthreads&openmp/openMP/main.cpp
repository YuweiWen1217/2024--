#include<iostream>
#include<windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
using namespace std;
#define N 200
#define NUM_THREADS 7
float** A = NULL;

long long head, tail, freq;


void A_init() {
    A = new float* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
    }
    for (int i = 0; i < N; i++) {
        A[i][i] = 1.0;
        for (int j = i + 1; j < N; j++) {
            A[i][j] = rand() % 1000;
        }

    }
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] += A[k][j];
                A[i][j] = (int)A[i][j] % 1000;
            }
        }
    }
}


void deleteA() {
    for (int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete A;
}


void print() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void LU() {
    int i, j, k;
    for (k = 0; k < N; k++) {
        for (j = k + 1; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for ( i = k + 1; i < N; i++) {
            for (j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void LU_omp() {
    int i =0 , j = 0, k = 0;
    float tmp = 0;
    #pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k,tmp)
    for (k = 0; k < N; k++) {
        #pragma omp single
        { tmp = A[k][k];
        for (j = k + 1; j < N; j++) {
            A[k][j] = A[k][j] / tmp;
        }
        A[k][k] = 1.0;
        }

        #pragma omp for
        for (i = k + 1; i < N; i++) {
            for (j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void LU_omp_guided() {
    int i = 0, j = 0, k = 0;
    float tmp = 0;
#pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k,tmp)
    for (k = 0; k < N; k++) {
#pragma omp single
        { tmp = A[k][k];
        for (j = k + 1; j < N; j++) {
            A[k][j] = A[k][j] / tmp;
        }
        }
        A[k][k] = 1.0;
#pragma omp for schedule(guided)
        for (i = k + 1; i < N; i++) {
            for (j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void LU_omp_dynamic() {
    int i = 0, j = 0, k = 0;
    float tmp = 0;
#pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k,tmp)
    for (k = 0; k < N; k++) {
#pragma omp single
        { tmp = A[k][k];
        for (j = k + 1; j < N; j++) {
            A[k][j] = A[k][j] / tmp;
        }
        }
        A[k][k] = 1.0;
#pragma omp for schedule(dynamic)
        for (i = k + 1; i < N; i++) {
            for (j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}



void avx_omp_static() {
    int i = 0, j = 0, k = 0;

#pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k)
    for (k = 0; k < N; k++) {
#pragma omp single
        {
            __m256 t1 = _mm256_set1_ps(A[k][k]);
            int j = 0;
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 t2 = _mm256_loadu_ps(&A[k][j]);
                t2 = _mm256_div_ps(t2, t1);
                _mm256_storeu_ps(&A[k][j], t2);
            }
            for (; j < N; j++) {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
        }
#pragma omp for schedule(static)
        for (i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void avx_omp_dynamic() {
    int i = 0, j = 0, k = 0;

#pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k)
    for (k = 0; k < N; k++) {
#pragma omp single
        {
            __m256 t1 = _mm256_set1_ps(A[k][k]);
            int j = 0;
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 t2 = _mm256_loadu_ps(&A[k][j]);
                t2 = _mm256_div_ps(t2, t1);
                _mm256_storeu_ps(&A[k][j], t2);
            }
            for (; j < N; j++) {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
        }
#pragma omp for schedule(dynamic)
        for (i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void avx_optimized() {
    for (int k = 0; k < N; k++) {
        __m256 t1 = _mm256_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 8 <= N; j += 8) {
            __m256 t2 = _mm256_loadu_ps(&A[k][j]);
            t2 = _mm256_div_ps(t2, t1);
            _mm256_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void avx_omp_guided() {
    int i = 0, j = 0, k = 0;

#pragma omp parallel  num_threads(NUM_THREADS),private(i,j,k)
    for (k = 0; k < N; k++) {
#pragma omp single
        {
            __m256 t1 = _mm256_set1_ps(A[k][k]);
            int j = 0;
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 t2 = _mm256_loadu_ps(&A[k][j]);
                t2 = _mm256_div_ps(t2, t1);
                _mm256_storeu_ps(&A[k][j], t2);
            }
            for (; j < N; j++) {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
        }
#pragma omp for schedule(guided)
        for (i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void cal(void(*func)()) {
    A_init();
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    func();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);

}


int main() {
    double sum1=0;
    double sum2=0;
    double sum3=0;
    double sum4=0;
    double sum5=0;
    double sum6=0;
    double sum7=0;
    double sum8=0;
    for(int i=0;i<10;i++){
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    cal(LU);
    sum1+= (tail - head) * 1000 / freq;
    deleteA();

    cal(LU_omp);
    sum2+=(tail - head) * 1000 / freq;
    deleteA();

    cal(LU_omp_guided);
    sum3+=(tail - head) * 1000 / freq;
    deleteA();

    cal(LU_omp_dynamic);
    sum4+=(tail - head) * 1000 / freq;
    deleteA();


    cal(avx_optimized);
    sum5+=(tail - head) * 1000 / freq;
    deleteA();

    cal(avx_omp_static);
    sum6+=(tail - head) * 1000 / freq;
    deleteA();

    cal(avx_omp_dynamic);
    sum7+=(tail - head) * 1000 / freq;
    deleteA();

    cal(avx_omp_guided);
    sum8+=(tail - head) * 1000 / freq;
    deleteA();
    }
    cout << "ordinary:" << sum1/10<< "ms" << endl;
    cout << "OpenMp_static:" <<sum2/10 << "ms" << endl;
    cout << "OpenMp_guided:" << sum3/10 << "ms" << endl;
    cout << "OpenMp_dynamic:" << sum4/10 << "ms" << endl;
    cout << "avx_ordinary:" << sum5/10 << "ms" << endl;
    cout << "avx_OpenMP_static:" << sum6/10 << "ms" << endl;
    cout << "avx_OpenMP_dynamic:" << sum7/10 << "ms" << endl;
    cout << "avx_OpenMP_guided:" << sum8/10 << "ms" << endl;
}
