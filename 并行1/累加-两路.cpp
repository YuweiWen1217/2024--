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
    int sum, sum1, sum2, a[n];
    auto sumtim = 0.0;
    // cout<<n<<endl;
    for (int f = 0; f < 10; f++)
    {
        sum = 0;
        sum1 = 0;
        sum2 = 0;
        for (int i = 0; i < n; i++)
        {
            a[i] = i;
        }
        QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for (int count = 0; count < 1000; count++)
        {
            for (int i = 0; i < n; i += 2)
            {
                sum1 += a[i];
                sum2 += a[i + 1];
            }
        }
        sum = sum1 + sum2;
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        sumtim += (tail - head) * 1000.0 / freq;
    }
    sumtim = sumtim / 10;
    cout << sumtim << endl;
}