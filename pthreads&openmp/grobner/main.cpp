#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<iomanip>
#include<windows.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
#include<pthread.h>
#include<omp.h>
using namespace std;

#define NUM_THREADS 7


struct threadParam_t {
	int t_id;
	int num;
};

const int maxsize = 3000;
const int maxrow = 60000;
const int numBasis = 100000;
pthread_mutex_t lock;
long long head, tail, freq;
map<int, int*>ans;

fstream RowFile("被消元行.txt", ios::in | ios::out);
fstream BasisFile("消元子.txt", ios::in | ios::out);


int gRows[maxrow][maxsize];
int gBasis[numBasis][maxsize];

int ifBasis[numBasis] = { 0 };

void reset() {
	memset(gRows, 0, sizeof(gRows));
	memset(gBasis, 0, sizeof(gBasis));
	memset(ifBasis, 0, sizeof(ifBasis));
	RowFile.close();
	BasisFile.close();
	RowFile.open("被消元行.txt", ios::in | ios::out);
	BasisFile.open("消元子.txt", ios::in | ios::out);
	ans.clear();
}

int readBasis() {
	for (int i = 0; i < numBasis; i++) {
		if (BasisFile.eof()) {
			cout << "读取消元子" << i - 1 << "行" << endl;
			return i - 1;
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
				ifBasis[row] = 1;
			}
			int index = pos / 32;
			int offset = pos % 32;
			gBasis[row][index] = gBasis[row][index] | (1 << offset);
		}
		flag = false;
		row = 0;
	}
}

int readRowsFrom(int pos) {
	if (RowFile.is_open())
		RowFile.close();
	RowFile.open("被消元行.txt", ios::in | ios::out);
	memset(gRows, 0, sizeof(gRows));
	string line;
	for (int i = 0; i < pos; i++) {
		getline(RowFile, line);
	}
	for (int i = pos; i < pos + maxrow; i++) {
		int tmp;
		getline(RowFile, line);
		if (line.empty()) {
			cout << "读取被消元行 " << i << " 行" << endl;
			return i;
		}
		bool flag = false;
		stringstream s(line);
		while (s >> tmp) {
			int index = tmp / 32;
			int offset = tmp % 32;
			gRows[i - pos][index] = gRows[i - pos][index] | (1 << offset);
			flag = true;
		}
	}
	cout << "read max rows" << endl;
	return -1;

}

int findfirst(int row) {
	int first;
	for (int i = maxsize - 1; i >= 0; i--) {
		if (gRows[row][i] == 0)
			continue;
		else {
			int pos = i * 32;
			int offset = 0;
			for (int k = 31; k >= 0; k--) {
				if (gRows[row][i] & (1 << k))
				{
					offset = k;
					break;
				}
			}
			first = pos + offset;
			return first;
		}
	}
	return -1;
}



void writeResult(ofstream& out) {
	for (auto it = ans.rbegin(); it != ans.rend(); it++) {
		int* result = it->second;
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

void GE() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);

	int num = (flag == -1) ? maxrow : flag;
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	for (int i = 0; i < num; i++) {
		while (findfirst(i)!= -1) {
			int first =findfirst(i);
			if (ifBasis[first]==1) {
				for (int j = 0; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];

				}
			}
			else {
				for (int j = 0; j < maxsize; j++) {
					gBasis[first][j] = gRows[i][j];
				}
				ifBasis[first] = 1;
				ans.insert(pair<int, int*>(first, gBasis[first]));
				break;
			}
		}
	}
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "Ordinary time:" <<setprecision(5) <<(tail - head) * 1000.000 / freq << "ms" << endl;
}

void GE_omp() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);
	//int i = 0, j = 0;
	int t_id = omp_get_thread_num();
	int num = (flag == -1) ? maxrow : flag;
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(guided)
	for (int i = 0; i < num; i++) {
		while (findfirst(i) != -1) {
			int first = findfirst(i);
			if (ifBasis[first] == 1) {
				for (int j = 0; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];

				}
			}
			else {
#pragma omp critical
				if (ifBasis[first] == 0) {
					for (int j = 0; j < maxsize; j++) {
						gBasis[first][j] = gRows[i][j];
					}
					ifBasis[first] = 1;
					ans.insert(pair<int, int*>(first, gBasis[first]));
				}
			}

		}
	}
	}
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "Omp time:" <<setprecision(8) << (tail - head) * 1000.000 / freq << "ms" << endl;
}


void AVX_GE() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);
	int num = (flag == -1) ? maxrow : flag;
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	for (int i = 0; i < num; i++) {
		while (findfirst(i) != -1) {
			int first = findfirst(i);
			if (ifBasis[first]==1) {
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					__m256i vj = _mm256_loadu_si256((__m256i*) & gBasis[first][j]);
					__m256i vx = _mm256_xor_si256(vij, vj);
					_mm256_storeu_si256((__m256i*) & gRows[i][j], vx);
				}
				for (; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];
				}
			}
			else {
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					_mm256_storeu_si256((__m256i*) & gBasis[first][j], vij);
				}
				for (; j < maxsize; j++) {
					gBasis[first][j] = gRows[i][j];
				}
				ifBasis[first] = 1;
				ans.insert(pair<int, int*>(first, gBasis[first]));
				break;

			}
		}
	}
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "AVX time:" <<setprecision(8) << (tail - head) * 1000.000 / freq << "ms" << endl;
}

void AVX_GE_omp() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);
	int num = (flag == -1) ? maxrow : flag;
	int i = 0, j = 0;
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
#pragma omp parallel  num_threads(NUM_THREADS),private(i,j)
#pragma omp for schedule(guided)
	for (i = 0; i < num; i++) {
		while (findfirst(i) != -1) {
			int first = findfirst(i);
			if (ifBasis[first]==1) {
				j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					__m256i vj = _mm256_loadu_si256((__m256i*) & gBasis[first][j]);
					__m256i vx = _mm256_xor_si256(vij, vj);
					_mm256_storeu_si256((__m256i*) & gRows[i][j], vx);
				}
				for (; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];
				}
			}
			else {
#pragma omp critical
				{j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					_mm256_storeu_si256((__m256i*) & gBasis[first][j], vij);
				}
				for (; j < maxsize; j++) {
					gBasis[first][j] = gRows[i][j];
				}
				ifBasis[first] = 1;
				ans.insert(pair<int, int*>(first, gBasis[first]));
				}
			}
		}
	}
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "AVX_omp time:" <<setprecision(8) << (tail - head) * 1000.000 / freq << "ms" << endl;
}

void* GE_lock_thread(void* param) {

	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id;
	int num = p->num;
	for (int i = t_id; i < num; i += NUM_THREADS) {
		while (findfirst(i) != -1) {
			int first = findfirst(i);
			if (ifBasis[first]==1) {
				for (int j = 0; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];

				}
			}
			else {
				pthread_mutex_lock(&lock);
				if (ifBasis[first]==1)
				{
					pthread_mutex_unlock(&lock);
					continue;
				}

				for (int j = 0; j < maxsize; j++) {
					gBasis[first][j] = gRows[i][j];
				}
				ifBasis[first] = 1;
				ans.insert(pair<int, int*>(first, gBasis[first]));
				pthread_mutex_unlock(&lock);
				break;
			}

		}
	}
	pthread_exit(NULL);
	return NULL;
}

void GE_pthread() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);

	int num = (flag == -1) ? maxrow : flag;

	pthread_mutex_init(&lock, NULL);

	pthread_t* handle = (pthread_t*)malloc(NUM_THREADS * sizeof(pthread_t));
	threadParam_t* param = (threadParam_t*)malloc(NUM_THREADS * sizeof(threadParam_t));

	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
		param[t_id].t_id = t_id;
		param[t_id].num = num;
		pthread_create(&handle[t_id], NULL, GE_lock_thread, &param[t_id]);
	}

	for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
		pthread_join(handle[t_id], NULL);
	}

	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GE_pthread time:" <<setprecision(8)<< (tail - head) * 1000.000 / freq << "ms" << endl;
	free(handle);
	free(param);
	pthread_mutex_destroy(&lock);
}

void* AVX_lock_thread(void* param) {

	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id;
	int num = p->num;

	for (int i = t_id; i  < num; i += NUM_THREADS) {
		while (findfirst(i) != -1) {
			int first = findfirst(i);
			if (ifBasis[first]==1) {
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					__m256i vj = _mm256_loadu_si256((__m256i*) & gBasis[first][j]);
					__m256i vx = _mm256_xor_si256(vij, vj);
					_mm256_storeu_si256((__m256i*) & gRows[i][j], vx);
				}
				for (; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];
				}
			}
			else {
				pthread_mutex_lock(&lock);
				if (ifBasis[first]==1)
				{
					pthread_mutex_unlock(&lock);
					continue;
				}
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					_mm256_storeu_si256((__m256i*) & gBasis[first][j], vij);
				}
				for (; j < maxsize; j++) {
					gBasis[first][j] = gRows[i][j];
				}
				ifBasis[first] = 1;
				ans.insert(pair<int, int*>(first, gBasis[first]));
				pthread_mutex_unlock(&lock);
				break;
			}
		}
	}
	pthread_exit(NULL);
	return NULL;
}

void AVX_pthread() {
	int begin = 0;
	int flag;
	flag = readRowsFrom(begin);

	int num = (flag == -1) ? maxrow : flag;

	pthread_mutex_init(&lock, NULL);

	pthread_t* handle = (pthread_t*)malloc(NUM_THREADS * sizeof(pthread_t));
	threadParam_t* param = (threadParam_t*)malloc(NUM_THREADS * sizeof(threadParam_t));

	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
		param[t_id].t_id = t_id;
		param[t_id].num = num;
		pthread_create(&handle[t_id], NULL, AVX_lock_thread, &param[t_id]);
	}

	for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
		pthread_join(handle[t_id], NULL);
	}

	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "AVX_pthread time:" <<setprecision(8) << (tail - head) * 1000.000 / freq << "ms" << endl;
	free(handle);
	free(param);
	pthread_mutex_destroy(&lock);
}

int main() {

		ofstream out("消元结果.txt");
		ofstream out1("消元结果(AVX).txt");
		ofstream out2("消元结果(GE_lock).txt");
		ofstream out3("消元结果(AVX_lock).txt");
		ofstream out4("消元结果(GE_omp).txt");
		ofstream out5("消元结果(AVX_omp).txt");
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

		// 串行
		readBasis();
		GE();
		writeResult(out);

		reset();

		// OpenMP
		readBasis();
		GE_omp();
		writeResult(out4);

		reset();

		// AVX+OPenMP
		readBasis();
		AVX_GE_omp();
		writeResult(out5);

		reset();

		// AVX
		readBasis();
		AVX_GE();
		writeResult(out1);

		reset();

		// Pthreads
		readBasis();
		GE_pthread();
		writeResult(out2);

		reset();

		//AVX + Pthreads
		readBasis();
		AVX_pthread();
		writeResult(out3);

		reset();
		out.close();
		out1.close();
		out2.close();
		out3.close();
}
