#include <iostream>
#include <vector>
#include <chrono>
#include <sys/time.h>
using namespace std;
int main()
{
    long long head, tail, freq;
    int n;
    cin >> n;
    int sum[n], b[n][n], a[n];
    auto sumtim = 0.0;
    timeval currentTime;
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
        gettimeofday(&currentTime, NULL);
        long long startTime = currentTime.tv_sec * 1000 + currentTime.tv_usec / 1000;
        for (int count = 0; count < 1000; count++)
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    sum[j] += (b[i][j] * a[j]);
        }
        gettimeofday(&currentTime, NULL);
        long long endTime = currentTime.tv_sec * 1000 + currentTime.tv_usec / 1000;
        sumtim += endTime - startTime;
    }
    sumtim = sumtim / 10;
    cout << sumtim << endl;
}
