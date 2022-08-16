#include <opencv2/core.hpp>
#include<iostream>
#include<cmath>
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"

using namespace cv;
using namespace std;

double** convolveH(double** a,double coef1,double coef2,int rows,int cols) {
    double** conv = new double* [rows];
    double** down = new double* [rows];
    for (int i = 0; i < rows; i++) {
        conv[i] = new double[cols+1];

    }
    for (int i = 0; i < rows; i++) {
        down[i] = new double[cols/2];

    }
    for (int i = 0; i < rows; i++) {
        conv[i][0] = a[i][0] * coef2;
        for (int j = 1; j < cols; j++) {
            conv[i][j] = a[i][j] * coef2 + a[i][j - 1] * coef1;
        }
        conv[i][cols] = a[i][cols-1] * coef1;
    }
    int k = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 1; j < cols; j += 2) {
            down[i][k] = conv[i][j];
            k++;
        }
        k = 0;
    }
    return down;
}
double** convolveV(double** a, double coef1, double coef2, int rows, int cols) {
    int rows2 = rows + 1;
    double** conv = new double* [rows+1];
    double** down = new double* [rows/2];
    for (int i = 0; i <= rows; i++) {
        conv[i] = new double[cols];

    }
    for (int i = 0; i < rows/2; i++) {
        down[i] = new double[cols];

    }
    for (int i = 0; i < cols; i++) {
        conv[0][i] = a[0][i] * coef2;
        for (int j = 1; j < rows; j++) {
            conv[j][i] = a[j][i] * coef2 + a[j-1][i] * coef1;
        }
        conv[rows][i] = a[rows - 1][i] * coef1;
    }
    int k = 0;
    
        for (int i = 1; i < rows; i+=2) 
        {
            for (int j = 0; j < cols; j++) {
                down[k][j] = conv[i][j];
            }
        k++;
        }
    return down;
}
void main() {
    Mat src = imread("lena.PNG");
    cvtColor(src, src, COLOR_BGR2GRAY, 0);
    Mat real(src.rows/2, src.cols/2, CV_64F);
    Mat real2(src.rows / 4, src.cols / 4, CV_64F);
    double** area = new double* [src.rows];
    double** result = new double* [src.rows];
    for (int i = 0; i < src.rows; i++) {
        area[i] = new double[src.cols];

    }
    for (int i = 0; i < src.rows; i++) {
        result[i] = new double[src.cols];

    }

    MatIterator_<uchar> it, end;
    int j = 0, k = 0;
    for (it = src.begin<uchar>(), end = src.end<uchar>(); it != end; ++it) {
        area[j][k] = (double)(*it);
        k++;
        if (k >= src.cols) {
            k = 0;
            j++;
        }
    }
    double** result1, ** result2, ** result3, ** result4, ** result5, ** result6,** r1, ** r2, ** r3, ** r4, ** r5, ** r6;
    result1 = convolveH(area, 1 / sqrt(2), -1 / sqrt(2), src.rows, src.cols);
    result2 = convolveH(area, 1 / sqrt(2), 1 / sqrt(2), src.rows, src.cols);
    result3 = convolveV(result1, 1 / sqrt(2), -1 / sqrt(2), src.rows, src.cols/2);
    result4 = convolveV(result1, 1 / sqrt(2), 1 / sqrt(2), src.rows, src.cols/2);
    result5 = convolveV(result2, 1 / sqrt(2), -1 / sqrt(2), src.rows, src.cols/2);
    result6 = convolveV(result2, 1 / sqrt(2), 1 / sqrt(2), src.rows, src.cols/2);

    MatIterator_<double> it1, end1;
    Mat m1 = Mat::ones(real.rows, real.cols, real.type());
    k = 0, j = 0;
    for (it1 = real.begin<double>(), end1 = real.end<double>(); it1 != end1; ++it1) {

        *it1 = result6[j][k];
        k++;
        if (k >= src.cols/2) {
            k = 0;
            j++;
        }
    }

    k = 0, j = 0;
    normalize(real, real, 0, 1, NORM_MINMAX);
    imwrite("out1.PNG", real * 255);

    for (it1 = real.begin<double>(), end1 = real.end<double>(); it1 != end1; ++it1) {
        *it1 = result5[j][k];
        k++;
        if (k >= src.cols / 2) {
            k = 0;
            j++;
        }
    }

    k = 0, j = 0;
    normalize(real, real, 0, 1, NORM_MINMAX);
    imwrite("out2.PNG", real * 255);

    for (it1 = real.begin<double>(), end1 = real.end<double>(); it1 != end1; ++it1) {
        *it1 = result4[j][k];
        k++;
        if (k >= src.cols / 2) {
            k = 0;
            j++;
        }
    }

    k = 0, j = 0;
    normalize(real, real, 0, 1, NORM_MINMAX);
    imwrite("out3.PNG", real * 255);
    for (it1 = real.begin<double>(), end1 = real.end<double>(); it1 != end1; ++it1) {
        *it1 = result3[j][k];
        k++;
        if (k >= src.cols / 2) {
            k = 0;
            j++;
        }
    }
    k = 0, j = 0;
    normalize(real, real, 0, 1, NORM_MINMAX);
    imwrite("out4.PNG", real * 255);

    r1 = convolveH(result6, 1 / sqrt(2), -1 / sqrt(2), src.rows/2, src.cols/2);
    r2 = convolveH(result6, 1 / sqrt(2), 1 / sqrt(2), src.rows / 2, src.cols / 2);
    r3 = convolveV(r1, 1 / sqrt(2), -1 / sqrt(2), src.rows / 2, src.cols / 4);
    r4 = convolveV(r1, 1 / sqrt(2), 1 / sqrt(2), src.rows / 2, src.cols / 4);
    r5 = convolveV(r2, 1 / sqrt(2), -1 / sqrt(2), src.rows / 2, src.cols / 4);
    r6 = convolveV(r2, 1 / sqrt(2), 1 / sqrt(2), src.rows / 2, src.cols / 4);

    k = 0, j = 0;
    for (it1 = real2.begin<double>(), end1 = real2.end<double>(); it1 != end1; ++it1) {
        *it1 = r6[j][k];
        k++;
        if (k >= src.cols / 4) {
            k = 0;
            j++;
        }
    }
    k = 0, j = 0;
    normalize(real2, real2, 0, 1, NORM_MINMAX);
    imwrite("out11.PNG", real2 * 255);

    for (it1 = real2.begin<double>(), end1 = real2.end<double>(); it1 != end1; ++it1) {
        *it1 = r5[j][k];
        k++;
        if (k >= src.cols / 4) {
            k = 0;
            j++;
        }
    }

    k = 0, j = 0;
    normalize(real2, real2, 0, 1, NORM_MINMAX);
    imwrite("out22.PNG", real2 * 255);
    for (it1 = real2.begin<double>(), end1 = real2.end<double>(); it1 != end1; ++it1) {
        *it1 = r4[j][k];
        k++;
        if (k >= src.cols / 4) {
            k = 0;
            j++;
        }
    }
    k = 0, j = 0;
    normalize(real2, real2, 0, 1, NORM_MINMAX);
    imwrite("out33.PNG", real2 * 255);

    for (it1 = real2.begin<double>(), end1 = real2.end<double>(); it1 != end1; ++it1) {

        *it1 = r3[j][k];
        k++;
        if (k >= src.cols / 4) {
            k = 0;
            j++;
        }
    }
    k = 0, j = 0;
    normalize(real2, real2, 0, 1, NORM_MINMAX);
    imwrite("out44.PNG", real2 * 255);


}

