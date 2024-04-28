#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <time.h>
#include <windows.h>
#include <tmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <immintrin.h>
using namespace std;

// 定义常量 T 表示批次数，maxsize 表示矩阵的最大列数，maxrow 表示矩阵的最大行数，numBasis 表示消元子的最大数量
const int T = 5;
const int maxsize = 3000;
const int maxrow = 3000;
const int numBasis = 90000;

// 用 map 存储消元子和被消元行的信息
// iToBasis 存储行号到消元子的映射
// iToFirst 存储消元子的首项位置
// ans 存储计算结果
map<int, int *> iToBasis;
map<int, int> iToFirst;
map<int, int *> ans;

// 打开被消元行和消元子文件
fstream RowFile("被消元行.txt", ios::in | ios::out);
fstream BasisFile("消元子.txt", ios::in | ios::out);

// 声明二维数组存储被消元行和消元子的矩阵
int gRows[maxrow][maxsize];	   // 被消元行矩阵
int gBasis[numBasis][maxsize]; // 消元子矩阵

// 重置函数，清空矩阵和文件，并重新打开文件
void reset()
{
	// 清空被消元行矩阵和消元子矩阵
	memset(gRows, 0, sizeof(gRows));
	memset(gBasis, 0, sizeof(gBasis));

	// 关闭文件流
	RowFile.close();
	BasisFile.close();

	// 重新打开文件流
	RowFile.open("被消元行.txt", ios::in | ios::out);
	BasisFile.open("消元子.txt", ios::in | ios::out);

	// 清空映射
	iToBasis.clear();
	iToFirst.clear();
	ans.clear();
}
// 读取消元子文件，并将其存储到 gBasis 数组和 iToBasis 映射中
void readBasis()
{
	// 逐行读取消元子文件
	for (int i = 0; i < maxrow; i++)
	{
		if (BasisFile.eof())
		{
			return;
		}
		string tmp;
		bool flag = false;
		int row = 0;
		getline(BasisFile, tmp);
		stringstream s(tmp);
		int pos;
		// 将每行中的元素按照空格分割并读取
		while (s >> pos)
		{
			if (!flag)
			{
				row = pos;
				flag = true;
				// 将行号和对应的消元子存入映射 iToBasis 中
				iToBasis.insert(pair<int, int *>(row, gBasis[row]));
			}
			// 将位置转换为数组索引和位偏移
			int index = pos / 32;
			int offset = pos % 32;
			// 将对应位置的位设置为 1
			gBasis[row][index] = gBasis[row][index] | (1 << offset);
		}
		flag = false;
		row = 0;
	}
}

// 从指定位置读取被消元行数据，并存储到 gRows 数组和 iToFirst 映射中
int readRowsFrom(int pos)
{
	// 清空 iToFirst 映射
	iToFirst.clear();
	// 如果文件已经打开，则关闭文件
	if (RowFile.is_open())
		RowFile.close();
	// 重新打开文件
	RowFile.open("被消元行.txt", ios::in | ios::out);
	// 初始化 gRows 数组为零
	memset(gRows, 0, sizeof(gRows));
	string line;
	// 定位到指定位置
	for (int i = 0; i < pos; i++)
	{
		getline(RowFile, line);
	}
	// 逐行读取数据，并存储到 gRows 数组和 iToFirst 映射中
	for (int i = pos; i < pos + maxsize; i++)
	{
		int tmp;
		getline(RowFile, line);
		if (line.empty())
		{
			cout << "End of File!" << endl;
			return i;
		}
		bool flag = false;
		stringstream s(line);
		// 逐个读取元素并存储到数组中
		while (s >> tmp)
		{
			if (!flag)
			{
				// 将当前行的首个元素和对应的行号存入映射 iToFirst 中
				iToFirst.insert(pair<int, int>(i - pos, tmp));
			}
			// 将位置转换为数组索引和位偏移
			int index = tmp / 32;
			int offset = tmp % 32;
			// 将对应位置的位设置为 1
			gRows[i - pos][index] = gRows[i - pos][index] | (1 << offset);
			flag = true;
		}
	}
	return -1;
}

// 更新指定行的首项位置
void update(int row)
{
	bool flag = false;
	// 从数组最后一个元素开始向前搜索
	for (int i = maxsize - 1; i >= 0; i--)
	{
		// 如果当前位置元素为零，则继续搜索下一个位置
		if (gRows[row][i] == 0)
			continue;
		else
		{
			// 找到第一个非零元素，更新首项位置
			if (!flag)
				flag = true;
			// 计算非零元素的位置
			int pos = i * 32;
			int offset = 0;
			for (int k = 31; k >= 0; k--)
			{
				if (gRows[row][i] & (1 << k))
				{
					offset = k;
					break;
				}
			}
			// 计算新的首项位置并更新映射 iToFirst
			int newfirst = pos + offset;
			iToFirst.erase(row);
			iToFirst.insert(pair<int, int>(row, newfirst));
			break;
		}
	}
	// 如果行中所有元素都为零，则从映射 iToFirst 中移除该行
	if (!flag)
		iToFirst.erase(row);
	return;
}

// 将结果写入输出流
void writeResult(ofstream &out)
{
	// 遍历结果映射 ans，从最后一个结果开始写入
	for (auto it = ans.rbegin(); it != ans.rend(); it++)
	{
		int *result = it->second;
		// 计算当前结果的数组大小
		int max = it->first / 32 + 1;
		// 从数组最后一个元素开始向前搜索
		for (int i = max; i >= 0; i--)
		{
			// 如果当前位置元素为零，则继续搜索下一个位置
			if (result[i] == 0)
				continue;
			int pos = i * 32;
			// 从当前位置向前搜索，找到非零元素的位置并写入输出流
			for (int k = 31; k >= 0; k--)
			{
				if (result[i] & (1 << k))
				{
					out << k + pos << " ";
				}
			}
		}
		// 写入换行符
		out << endl;
	}
}

void serial()
{
	// 当前批次的起始行索引
	int currentBatchStart = 0;
	// 标志，用于指示是否到达文件末尾
	int endFlag;
	// 循环直到处理完所有批次
	while (true)
	{
		// 从文件中读取被消元行数据，并返回标志，指示是否到达文件末尾
		endFlag = readRowsFrom(currentBatchStart);
		// 计算当前批次的行数
		int numRows = (endFlag == -1) ? maxsize : endFlag;
		// 遍历当前批次的每一行
		for (int i = 0; i < numRows; i++)
		{
			// 当当前行存在待消元元素时，执行消元操作
			while (iToFirst.find(i) != iToFirst.end())
			{
				// 获取当前行的首个待消元元素的列索引
				int firstPos = iToFirst.find(i)->second;
				// 如果该列索引对应的消元子存在，则执行消元操作
				if (iToBasis.find(firstPos) != iToBasis.end())
				{
					// 获取对应的消元子数组
					int *basis = iToBasis.find(firstPos)->second;
					// 将当前行与消元子进行异或操作，执行消元
					for (int j = 0; j < maxsize; j++)
					{
						gRows[i][j] = gRows[i][j] ^ basis[j];
					}
					// 更新当前行的首个待消元元素的位置
					update(i);
				}
				// 如果该列索引对应的消元子不存在，则将当前行作为新的消元子
				else
				{
					// 将当前行作为新的消元子，添加到消元子列表中
					for (int j = 0; j < maxsize; j++)
					{
						gBasis[firstPos][j] = gRows[i][j];
					}
					// 更新消元子映射和结果映射
					iToBasis.insert(pair<int, int *>(firstPos, gBasis[firstPos]));
					ans.insert(pair<int, int *>(firstPos, gBasis[firstPos]));
					// 移除当前行的待消元元素
					iToFirst.erase(i);
				}
			}
		}
		// 如果已经处理完所有行，则更新当前批次的起始行索引
		if (endFlag == -1)
			currentBatchStart += maxsize;
		// 如果还有剩余行未处理，则结束循环
		else
			break;
	}
}

void AVX()
{
	int begin = 0; // 起始位置
	int flag;	   // 标志位

	while (true)
	{
		flag = readRowsFrom(begin);				 // 读取被消元行数据
		int num = (flag == -1) ? maxsize : flag; // 计算需要处理的行数
		for (int i = 0; i < num; i++)
		{
			while (iToFirst.find(i) != iToFirst.end())
			{
				int first = iToFirst.find(i)->second; // 获取首项
				if (iToBasis.find(first) != iToBasis.end())
				{
					int *basis = iToBasis.find(first)->second; // 获取对应消元子
					int j = 0;
					for (; j + 8 < maxsize; j += 8)
					{
						// 使用 AVX 指令集加载被消元行和消元子，并执行异或操作
						__m256i vij = _mm256_loadu_si256((__m256i *)&gRows[i][j]);
						__m256i vj = _mm256_loadu_si256((__m256i *)&basis[j]);
						__m256i vx = _mm256_xor_si256(vij, vj);
						_mm256_storeu_si256((__m256i *)&gRows[i][j], vx);
					}
					for (; j < maxsize; j++)
					{
						// 使用普通方式进行异或操作
						gRows[i][j] = gRows[i][j] ^ basis[j];
					}
					update(i); // 更新首项
				}
				else
				{
					int j = 0;
					for (; j + 8 < maxsize; j += 8)
					{
						// 使用 AVX 指令集加载被消元行，并将其存储为消元子
						__m256i vij = _mm256_loadu_si256((__m256i *)&gRows[i][j]);
						_mm256_storeu_si256((__m256i *)&gBasis[first][j], vij);
					}
					for (; j < maxsize; j++)
					{
						// 使用普通方式将被消元行存储为消元子
						gBasis[first][j] = gRows[i][j];
					}
					// 更新相关映射
					iToBasis.insert(pair<int, int *>(first, gBasis[first]));
					ans.insert(pair<int, int *>(first, gBasis[first]));
					iToFirst.erase(i);
				}
			}
		}
		// 判断是否需要继续处理下一批数据
		if (flag == -1)
			begin += maxsize;
		else
			break;
	}
}

void SSE()
{
	// 初始化起始位置和标志位
	int start_position = 0;
	int end_flag;

	// 主循环，处理被消元行的批次
	while (true)
	{
		// 读取当前批次的被消元行
		end_flag = readRowsFrom(start_position);

		// 计算当前批次中的行数
		int num_rows = (end_flag == -1) ? maxsize : end_flag;

		// 对于每一行
		for (int row = 0; row < num_rows; row++)
		{
			// 循环处理每一个非零消元子
			while (iToFirst.find(row) != iToFirst.end())
			{
				// 获取当前行的首项对应的消元子
				int first_elem = iToFirst.find(row)->second;

				// 如果该消元子存在于消元子映射中
				if (iToBasis.find(first_elem) != iToBasis.end())
				{
					// 加载当前行和对应消元子到 SSE 寄存器中
					int *basis = iToBasis.find(first_elem)->second;
					int col = 0;
					for (; col + 4 < maxsize; col += 4)
					{
						__m128i row_vec = _mm_loadu_si128((__m128i *)&gRows[row][col]);
						__m128i basis_vec = _mm_loadu_si128((__m128i *)&basis[col]);
						__m128i result_vec = _mm_xor_si128(row_vec, basis_vec);
						_mm_storeu_si128((__m128i *)&gRows[row][col], result_vec);
					}
					// 处理剩余部分
					for (; col < maxsize; col++)
					{
						gRows[row][col] = gRows[row][col] ^ basis[col];
					}
					// 更新当前行
					update(row);
				}
				// 如果当前行的首项对应的消元子不存在于消元子映射中
				else
				{
					// 存储当前行到 gBasis 数组中
					int col = 0;
					for (; col + 4 < maxsize; col += 4)
					{
						__m128i row_vec = _mm_loadu_si128((__m128i *)&gRows[row][col]);
						_mm_storeu_si128((__m128i *)&gBasis[first_elem][col], row_vec);
					}
					// 处理剩余部分
					for (; col < maxsize; col++)
					{
						gBasis[first_elem][col] = gRows[row][col];
					}
					// 将新的消元子插入到映射中，并从 iToFirst 中删除当前行
					iToBasis.insert(pair<int, int *>(first_elem, gBasis[first_elem]));
					ans.insert(pair<int, int *>(first_elem, gBasis[first_elem]));
					iToFirst.erase(row);
				}
			}
		}
		// 更新起始位置
		if (end_flag == -1)
			start_position += maxsize;
		else
			break;
	}
}
int main()
{
	return 0;
}
