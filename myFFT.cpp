#define N 16                 /*FFT点数*/
#define pi 3.141592653589793 /*圆周率*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "myFFT.h"
//数据交换函数
void swap(double &a, double &b)
{
    double t;
    t = a;
    a = b;
    b = t;
}
//位码倒置函数
void bitrp(double xreal[], double ximag[], int n)
{
    // 位反转置换 Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
    {
        p++;
    }
    for (i = 0; i < n; i++)
    {
        a = i;
        b = 0;
        for (j = 0; j < p; j++)
        {
            b = (b << 1) + (a & 1); // b = b * 2 + a % 2;
            a >>= 1;                // a = a / 2;
        }
        if (b > i)
        {
            swap(xreal[i], xreal[b]);
            swap(ximag[i], ximag[b]);
        }
    }
}
void myFFT(double xreal[], double ximag[], int n)
{
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    double treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    bitrp(xreal, ximag, n);
    for (m = 2; m <= n; m *= 2)
    {
        for (k = 0; k < n; k += m)
        {
            for (j = 0; j < m / 2; j++)
            {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;
                treal = cos(2 * pi / n * (t)) * xreal[index2] - sin(2 * pi / n * (t)) * ximag[index2];
                timag = cos(2 * pi / n * (t)) * ximag[index2] + sin(2 * pi / n * (t)) * xreal[index2];
                ureal = xreal[index1];
                uimag = ximag[index1];
                xreal[index1] = ureal + treal;
                ximag[index1] = uimag + timag;
                xreal[index2] = ureal - treal;
                ximag[index2] = uimag - timag;
            }
        }
    }
}
//输出数据进行测试
int main()
{
    double data_real[N], data_imag[N];
    int len = N;
    for (int i = 0; i < N / 2; i++)
    {
        data_real[i] = 1;
        data_imag[i] = 0;
    }
    for (int i = N / 2; i < N; i++)
    {
        data_real[i] = 0;
        data_imag[i] = 0;
    }
    myFFT(data_real, data_imag, len);
    for (int i = 0; i < len; i++)
    {
        if (data_imag[i] >= 0.0)
        {
            printf("%.15f + %.15fi\n", data_real[i], data_imag[i]);
        }
        else
        {
            printf("%.15f - %.15fi\n", data_real[i], fabs(data_imag[i]));
        }
    }
    return 0;
}
