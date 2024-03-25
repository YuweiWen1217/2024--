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
    int a[n];
    auto sumtim = 0.0;
    // cout<<n<<endl;
    for (int f = 0; f < 10; f++)
    {
        for (int i = 0; i < n; i++)
        {
            a[i] = i;
        }
        QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for (int count = 0; count < 1000; count++)
        {
            for (m = n; m > 1; m /= 2)
                for (i = 0; i > m / 2; i++)
                    a[i] = a[2 * i] + a[2 * i + 1];
        }
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
        sumtim += (tail - head) * 1000.0 / freq;
    }
    sumtim = sumtim / 10;
    cout << sumtim << endl;
}