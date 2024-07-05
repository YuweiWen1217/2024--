#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include <iomanip>
#include"cuda_runtime.h"
#include"device_launch_parameters.h"
#include"device_functions.h"
#include<windows.h>

using namespace std;
const int N = 1024;
const int BLOCK_SIZE = 1024;
float **A;
float elm[N][N] = { 0 };

float* result = new float[N * N];
float* temp = new float[N * N];

long long head, tail, freq;

void A_init() {     //����ĳ�ʼ��
    A = new float* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
    }
    for (int i = 0; i < N; i++) {
        A[i][i] = 1.0;
        for (int j = i + 1; j < N; j++) {
            A[i][j] = rand() % 5000;
        }

    }
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] += A[k][j];
                A[i][j] = (int)A[i][j] % 5000;
            }
        }
    }
}

void copy_matrix() {
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			temp[i * N + j] = A[i][j];
		}
	}
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            elm[i][j] = A[i][j];
        }
    }
}

void transfer_result() {
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			elm[i][j] = result[i * N + j];
		}
	}
}

void print_result() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << elm[i][j] << " ";
        }
        cout << endl;
    }
}

void LU() {    //��ͨ��Ԫ�㷨
    for (int k = 0; k < N; k++) {
        for (int j = k + 1; j < N; j++) {
            elm[k][j] = elm[k][j] / elm[k][k];
        }
        elm[k][k] = 1.0;

        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1; j < N; j++) {
                elm[i][j] = elm[i][j] - elm[i][k] * elm[k][j];
            }
            elm[i][k] = 0;
        }
    }
}


__global__ void division_kernel(float* data, int k, int N) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;//�����߳�����
    int element = data[k * N + k];
    int temp = data[k * N + tid];
    data[k * N + tid] = (float)temp / element;
    return;
}

__global__ void eliminate_kernel(float* data, int k, int N) {
    int tx = blockDim.x * blockIdx.x + threadIdx.x;
    if (tx == 0)
        data[k * N + k] = 1.0;//�Խ���Ԫ����Ϊ 1
    int row = k + 1 + blockIdx.x;//ÿ���鸺��һ��
    while (row < N) {
        int tid = threadIdx.x;
        while (k + 1 + tid < N) {
            int col = k + 1 + tid;
            float temp_1 = data[(row * N) + col];
            float temp_2 = data[(row * N) + k];
            float temp_3 = data[k * N + col];
            data[(row * N) + col] = temp_1 - temp_2 * temp_3;
            tid = tid + blockDim.x;
        }
        __syncthreads();//����ͬ��
        if (threadIdx.x == 0) {
            data[row * N + k] = 0;
        }
        row += gridDim.x;
    }
    return;
}

int main() {
    A_init();
    copy_matrix();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    LU();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "CPU_LU:" << (tail - head) * 1000 / freq << "ms" << endl;
    print_result();       //ƽ����ȥ��ִ�к�ʱ�����

	cudaError_t ret;//���ڴ����飬�� CUDA �ӿڵ��óɹ��᷵�� cudaSucess
	float* gpudata;
	int size = N * N * sizeof(float);
	ret = cudaMalloc(&gpudata, size);//�����Դ�ռ�
	if (ret != cudaSuccess) {
		printf("cudaMalloc gpudata failed!\n");

	}
	ret = cudaMemcpy(gpudata, temp, size, cudaMemcpyHostToDevice);//�����ݴ����� GPU ��
	if (ret != cudaSuccess) {
		printf("cudaMemcpyHostToDevice failed!\n");
	}

	dim3 dimBlock(BLOCK_SIZE, 1);//�߳̿�
	dim3 dimGrid(1, 1);//�߳�����

	cudaEvent_t start, stop;//��ʱ��
	float elapsedTime = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);//��ʼ��ʱ

	for (int k = 0; k < N; k++) {
		division_kernel << <dimGrid, dimBlock >> > (gpudata, k, N);//�����������ĺ˺���
		cudaDeviceSynchronize();//CPU �� GPU ֮���ͬ������
		ret = cudaGetLastError();
		if (ret != cudaSuccess) {
			printf("division_kernel failed, %s\n", cudaGetErrorString(ret));
		}
		eliminate_kernel << <dimGrid, dimBlock >> > (gpudata, k,N);//������ȥ����ĺ˺���
		cudaDeviceSynchronize();
		ret = cudaGetLastError();
		if (ret != cudaSuccess) {
			printf("eliminate_kernel failed, %s\n", cudaGetErrorString(ret));
		}
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);//ֹͣ��ʱ
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("GPU_LU:%f ms\n", elapsedTime);
	cudaError_t cudaStatus2 = cudaGetLastError();
	if (cudaStatus2 != cudaSuccess) {
		fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus2));
	}
	ret = cudaMemcpy(result, gpudata, size, cudaMemcpyDeviceToHost);//�����ݴ��� CPU ��
	if (ret != cudaSuccess) {
		printf("cudaMemcpyDeviceToHost failed!\n");
	}
	transfer_result();
    print_result();
	cudaFree(gpudata);//�ͷ��Դ�ռ䣬�� CUDA �ӿڷ���Ŀռ������ cudaFree �ͷ�
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
}

