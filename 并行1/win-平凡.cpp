#include <iostream>
#include <vector>
#include <chrono>
#include <windows.h>
using namespace std;
int main()
{
    long long head, tail, freq;
    int n;
    cin >> n;
    int sum[n], b[n][n], a[n];
    auto sumtim = 0.0;
    // cout<<n<<endl;
    for (int f = 0; f < 10; f++)
    {
        for (int i = 0; i < n; i++)
        {
            sum[i] = 0;
            a[i] = 2 * i;
            for (int j = 0; j < 0; j++)
                b[i][j] = i + j;
        }
        QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for (int count = 0; count < 1000; count++)
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    sum[i] += (b[j][i] * a[i]);
        }
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        sumtim += (tail - head) * 1000.0 / freq;
    }
    sumtim = sumtim / 10;
    cout << sumtim << endl;
}