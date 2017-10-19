//
//  color_adjust.cpp
//  test
//
//  Created by 郭秀艳 on 2017/6/21.
//  Copyright © 2017年 郭秀艳. All rights reserved.
//

#include "fft2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>
#define intsize sizeof(int)
#define complexsize sizeof(complex)
#define PI 3.1415926

#define IMG_WIDTH 512
#define IMG_HEIGHT 512

#define IMG_WIDTH_512 512
#define IMG_HEIGHT_512 512

#define IMG_WIDTH_64 64
#define IMG_HEIGHT_64 64

#define SIZE 512

complex *complexSrc;
complex *complexDt;
complex *complexDtPct;
double *imgSrc;
double *magPlusOne;
double *normalise;
double *logE;
double *highpass;


int *a,*b;
int nLen, init_nLen ,mLen, init_mLen, N, M;
FILE *dataFile;

complex *A,*A_In,*W;

cv::Mat img;

void readData()
{
    int i,j;
//    int data[16][16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
//    dataFile = fopen("dataIn.txt","r");
//    fscanf(dataFile,"%d %d", &init_mLen, &init_nLen);
    init_nLen = SIZE;
    init_mLen = SIZE;
    
    M = calculate_M(init_mLen);
    N = calculate_M(init_nLen);
    nLen = (int)pow(2,N);
    mLen = (int)pow(2,M);
    A_In = (complex *)malloc(sizeof(complex)*nLen*mLen);
    
    for(i=0; i<init_mLen; i++)
    {
        for(j=0; j<init_nLen; j++)
        {
            A_In[i*nLen+j].real = imgSrc[i*nLen+j];
            A_In[i*nLen+j].image = 0.0;
//            A_In[i*nLen+j].real = complexSrc[i*nLen+j].real;
//            A_In[i*nLen+j].image = complexSrc[i*nLen+j].image;
        }
    }
//    fclose(dataFile);
    
    for(i=0; i<mLen; i++)
    {
        for(j=init_nLen; j<nLen; j++)
        {
            A_In[i*nLen+j].real = 0.0;
            A_In[i*nLen+j].image = 0.0;
        }
    }
    
    for(i=init_mLen; i<mLen; i++)
    {
        for(j=0; j<init_nLen; j++)
        {
            A_In[i*nLen+j].real = 0.0;
            A_In[i*nLen+j].image = 0.0;
        }
    }
    
    return;
    
    printf("Reading initial datas:\n");
    for(i=0; i<init_mLen; i++)
    {
        for(j=0; j<init_nLen; j++)
        {
            if(A_In[i*nLen+j].image < 0)
            {
                printf("%f%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
            else
            {
                printf("%f+%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
        }
        printf("\n");
    }
    
    printf("\n");
    
    printf("Reading formal datas:\n");
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            if(A_In[i*nLen+j].image < 0)
            {
                printf("%f%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
            else
            {
                printf("%f+%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
        }
        printf("\n");
    }
}


void fft(int fft_nLen, int fft_M)
{
    int i;
    int lev,dist,p,t;
    complex B;
    
    W = (complex *)malloc(complexsize*fft_nLen/2);
    
    for(lev=1; lev<=fft_M; lev++)
    {
        dist = (int)pow(2,lev-1);
        for(t=0; t<dist; t++)
        {
            p = t*(int)pow(2,fft_M-lev);
            W[p].real = (double)cos(2*PI*p/fft_nLen);
            W[p].image = (double)(-1*sin(2*PI*p/fft_nLen));
            for(i=t; i<fft_nLen; i=i+(int)pow(2,lev))
            {
                B = Add(A[i],Mul(A[i+dist],W[p]));
                A[i+dist] = Sub(A[i],Mul(A[i+dist],W[p]));
                A[i].real = B.real;
                A[i].image = B.image;
            }
        }
    }
    
    free(W);
}


void printResult_fft()
{
    int i,j;
    
    printf("Output FFT results:\n");
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            if(A_In[i*nLen+j].image < 0)
            {
                printf("%f%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
            else
            {
                printf("%f+%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
        }
        printf("\n");
    }
}

void printResult_Ifft()
{
    int i,j;
    
    printf("Output IFFT results:\n");
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            if(A_In[i*nLen+j].image < 0)
            {
                printf("%f%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
            else
            {
                printf("%f+%fi\n",A_In[i*nLen+j].real,A_In[i*nLen+j].image);
            }
        }
        printf("\n");
    }
    
    free(A_In);
}

int calculate_M(int len)
{
    int i;
    int k;
    
    i = 0;
    k = 1;
    while(k < len)
    {
        k = k*2;
        i++;
    }
    
    return i;
}

void reverse(int len, int M)
{
    int i,j;
    
    a = (int *)malloc(intsize*M);
    b = (int *)malloc(intsize*len);
    
    for(i=0; i<M; i++)
    {
        a[i] = 0;
    }
    
    b[0] = 0;
    for(i=1; i<len; i++)
    {
        j = 0;
        while(a[j] != 0)
        {
            a[j] = 0;
            j++;
        }
        
        a[j] = 1;
        b[i] = 0;
        for(j=0; j<M; j++)
        {
            b[i] = b[i]+a[j]*(int)pow(2,M-1-j);
        }
    }
}

complex Add(complex c1, complex c2)
{
    complex c;
    c.real = c1.real+c2.real;
    c.image = c1.image+c2.image;
    return c;
}

complex Sub(complex c1, complex c2)
{
    complex c;
    c.real = c1.real-c2.real;
    c.image = c1.image-c2.image;
    return c;
}

complex Mul(complex c1, complex c2)
{
    complex c;
    c.real = c1.real*c2.real-c1.image*c2.image;
    c.image = c1.real*c2.image+c2.real*c1.image;
    return c;
}

void Ifft()
{
    int i,j;
    
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            A_In[i*nLen+j].image = -A_In[i*nLen+j].image;
        }
    }
    
    A = (complex *)malloc(complexsize*nLen);
    reverse(nLen,N);
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            A[j].real = A_In[i*nLen+b[j]].real;
            A[j].image = A_In[i*nLen+b[j]].image;
        }
        
        fft(nLen,N);
        for(j=0; j<nLen; j++)
        {
            A_In[i*nLen+j].real = A[j].real/nLen;
            A_In[i*nLen+j].image = A[j].image/nLen;
        }
    }
    free(A);
    free(a);
    free(b);
    
    A = (complex *)malloc(complexsize*mLen);
    reverse(mLen,M);
    for(i=0; i<nLen; i++)
    {
        for(j=0; j<mLen; j++)
        {
            A[j].real = A_In[b[j]*nLen+i].real;
            A[j].image = A_In[b[j]*nLen+i].image;
        }
        
        fft(mLen,M);
        for(j=0; j<mLen; j++)
        {
            A_In[j*nLen+i].real = A[j].real/mLen;
            A_In[j*nLen+i].image = A[j].image/mLen;
        }
    }
    free(A);
    free(a);
    free(b);
}

void shiyangfft2d()
{
    int i,j;
    readData();
#if 1
    A = (complex *)malloc(complexsize*nLen);
    reverse(nLen,N);
    for(i=0; i<mLen; i++)
    {
        for(j=0; j<nLen; j++)
        {
            A[j].real = A_In[i*nLen+b[j]].real;
            A[j].image = A_In[i*nLen+b[j]].image;
        }
        
        fft(nLen,N);
        for(j=0; j<nLen; j++)
        {
            A_In[i*nLen+j].real = A[j].real;
            A_In[i*nLen+j].image = A[j].image;
        }
    }
    
    free(a);
    free(b);
    free(A);
    
    A = (complex *)malloc(complexsize*mLen);
    reverse(mLen,M);
    for(i=0; i<nLen; i++)
    {
        for(j=0; j<mLen; j++)
        {
            A[j].real = A_In[b[j]*nLen+i].real;
            A[j].image = A_In[b[j]*nLen+i].image;
        }
        
        fft(mLen,M);
        for(j=0; j<mLen; j++)
        {
            A_In[j*nLen+i].real = A[j].real;
            A_In[j*nLen+i].image = A[j].image;
        }
    }
    free(A);
//    printResult_fft();
    return;
    if (true) {
        doMagPlusOne();
        doLogE();
        doNormalise();
        doHighPass();
    }
    
    if (true) {
        doHighPass(A_In);
    }
#endif
    Ifft();
    
    char *imgChar = new char[SIZE*SIZE];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            imgChar[idx] = (char) (255*(A_In[idx].real/255.0));
        }
    }
    
    cvNamedWindow("shiyangifft", 1);
    IplImage *ImageIFFT = cvCreateImage(cvSize(SIZE,SIZE), IPL_DEPTH_8U, 1);
    ImageIFFT->imageData = imgChar;
    cvShowImage("shiyangifft", ImageIFFT);
//    cvWaitKey();
    
    printResult_Ifft();
}

void doMagPlusOne()
{
    magPlusOne = new double[SIZE*SIZE];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            double sqroot = sqrt(A_In[idx].real*A_In[idx].real+A_In[idx].image*A_In[idx].image);
            magPlusOne[idx] = 1.0 + sqroot;
        }
    }
}

void doLogE()
{
    logE = new double[SIZE*SIZE];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            logE[idx] = log(magPlusOne[idx]);
        }
    }
}

void doNormalise()
{
    double max = -1.0;
    double min = 10000.0;
    normalise = new double[SIZE*SIZE];
    char *normaliseChar = new char[SIZE*SIZE];
    
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            if (logE[idx] > max) {
                max = logE[idx];
            }
            if (logE[idx] < min) {
                min = logE[idx];
            }
        }
    }
    
    double delta = max - min;
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            normalise[idx] = (logE[idx] - min)/delta;
            if (normalise[idx] < 0.0) {
                printf("shiyang error normalise wrong");
            }
            normaliseChar[idx] = (char) (255.0*normalise[idx]);
        }
    }
    
    cvNamedWindow("normaliseChar", 1);
    IplImage *ImageNMLS = cvCreateImage(cvSize(SIZE,SIZE), IPL_DEPTH_8U, 1);
    ImageNMLS->imageData = normaliseChar;
    cvShowImage("normaliseChar", ImageNMLS);
}

void doHighPass(complex *_in)
{
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            bool corner0 = (i<=init_mLen/2) && (j<=init_nLen/2);
            bool corner1 = (i<=init_mLen/2) && (j>init_nLen/2);
            bool corner2 = (i>init_mLen/2) && (j<=init_nLen/2);
            bool corner3 = (i>init_mLen/2) && (j>init_nLen/2);
            int di=i,dj=j;
            if (corner0) {
                di = i;
                dj = j;
            }
            if (corner1) {
                di = i;
                dj = init_nLen-j;
            }
            if (corner2) {
                di = init_mLen-i;
                dj = j;
            }
            if (corner3) {
                di = init_mLen-i;
                dj = init_nLen-j;
            }
            double dist = sqrt(di*di+dj*dj);
            int idx = i*nLen+j;
//            int idx = i*mLen+j;
//            bool need =dist<SIZE/10;
//            if (need) {
//                _in[idx].real *= 0.0;
//                _in[idx].image *= 0.0;
//            }
            double alpha = dist/SIZE/sqrt(2.0);
            _in[idx].real *= (alpha);
            _in[idx].image *= (alpha);
            
            
        }
    }
}

void doHighPass()
{
    highpass = new double[SIZE*SIZE];
    char *highpassChar = new char[SIZE*SIZE];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            bool corner0 = (i<=init_mLen/2) && (j<=init_nLen/2);
            bool corner1 = (i<=init_mLen/2) && (j>init_nLen/2);
            bool corner2 = (i>init_mLen/2) && (j<=init_nLen/2);
            bool corner3 = (i>init_mLen/2) && (j>init_nLen/2);
            int di=SIZE,dj=SIZE;
            if (corner0) {
                di = i;
                dj = j;
            }
            if (corner1) {
                di = i;
                dj = init_nLen-j;
            }
            if (corner2) {
                di = init_mLen-i;
                dj = j;
            }
            if (corner3) {
                di = init_mLen-i;
                dj = init_nLen-j;
            }
            double dist = sqrt(di*di+dj*dj);
            int idx = i*nLen+j;
            //            int idx = i*mLen+j;
            bool need =dist<SIZE/10;
            if (need) {
                highpass[idx] *= 0.0;
            } else {
                highpass[idx] = normalise[idx];
            }
            highpassChar[idx] = (char) (255.0*highpass[idx]);
        }
    }
    
    cvNamedWindow("highpassChar", 1);
    IplImage *ImageHP = cvCreateImage(cvSize(SIZE,SIZE), IPL_DEPTH_8U, 1);
    ImageHP->imageData = highpassChar;
    cvShowImage("highpassChar", ImageHP);
}

void shiyangcvdft2dimage()
{
    cv::Mat planes[] = {cv::Mat_<float>(img), cv::Mat::zeros(img.size(), CV_32F) };
    cv::Mat complexI;
    merge(planes, 2, complexI);
    dft(complexI, complexI);
    split(complexI, planes);
    

    complexDt = (complex *)malloc(sizeof(complex)*nLen*mLen);
    complexDtPct = (complex *)malloc(sizeof(complex)*nLen*mLen);
    
    printf("\n");
    int errCount = 0;
    for (int j=0; j<SIZE; j++) {
        for (int i=0; i<SIZE; i++) {
            int idx = j*SIZE+i;
            double dtReal = fabs(A_In[idx].real - planes[0].at<float>(j, i));
            double dtImag = fabs(A_In[idx].image - planes[1].at<float>(j, i));
            complexDt[idx].real = dtReal;
            complexDt[idx].image = dtImag;
            
            complexDtPct[idx].real = dtReal/fabs(planes[0].at<float>(j, i));
            complexDtPct[idx].image = dtImag/fabs(planes[1].at<float>(j, i));
            
            if (complexDtPct[idx].real > 0.1 || complexDtPct[idx].image > 0.1) {
                printf("index=%d\n", idx);
                printf("ours: real=%f, imag=%f\n", A_In[idx].real, A_In[idx].image);
                printf("opcv: real=%f, imag=%f\n", planes[0].at<float>(j, i), planes[1].at<float>(j, i));
                errCount++;
            }
        }
    }
    printf("total errCount = %d\n", errCount);
}

void shiyangcvdft2d()
{
    cv::Mat planes0[] = {cv::Mat::zeros(cvSize(16,16), CV_32F), cv::Mat::zeros(cvSize(16,16), CV_32F) };
    cv::Mat planes1[] = {cv::Mat::zeros(cvSize(16,16), CV_32F), cv::Mat::zeros(cvSize(16,16), CV_32F) };
    for (int j=0; j<16; j++) {
        for (int i=0; i<16; i++) {
            planes0[0].at<float>(j, i) = j*16+i;
            planes1[0].at<float>(i, j) = j*16+i;
        }
    }
    cv::Mat complexI0;
    merge(planes0, 2, complexI0);
    
    dft(complexI0, complexI0);
    split(complexI0, planes0);

    cv::Mat complexI1;
    merge(planes1, 2, complexI1);
    
    dft(complexI1, complexI1);
    split(complexI1, planes1);
    
    for (int j=0; j<16; j++) {
        for (int i=0; i<16; i++) {
            int idx=j*16+i;
            printf("0i=%d,real=%f,imag=%f\n",idx, planes0[0].at<float>(j, i), planes0[1].at<float>(j, i));
        }
    }
    for (int j=0; j<16; j++) {
        for (int i=0; i<16; i++) {
            int idx=j*16+i;
            printf("1i=%d,real=%f,imag=%f\n",idx, planes1[0].at<float>(j, i), planes1[1].at<float>(j, i));
        }
    }
}

double areal[16] = {
    -1.0,        0.0,          1.0,   -10.0,
    100.0,     5556.0,      -8888.0,    10.0,
    -99999.0, 765432.876, 1234567890.0,   888.0,
    0.0,  -908765.0,     876578.0,  2348.0
};

double aimage[16] = {
    98765.0,   -12355.0,          0.0,     1.0,
    88888.0,  -777777.0,          9.0, 45444.0,
    66677.0,    111.777,      888.222, 888.444,
    22222.7,    -2348.0,    -384623.0,  2736.0
};

void shiyanggetImgSrc()
{
    imgSrc = new double[SIZE*SIZE];
    bool isReadPNG = true;
    if (isReadPNG) {
        img = cv::imread("/Users/admin/Desktop/beautyTest/dilireba.png", CV_LOAD_IMAGE_GRAYSCALE);
        for (int j=0; j<SIZE; j++) {
            for (int i=0; i<SIZE; i++) {
                int idx = j*SIZE+i;
                imgSrc[idx] = img.at<unsigned char>(j, i);
                if (idx<20) {
                    printf("imgSrc %f\n", imgSrc[idx]);
                }
            }
        }
    } else {
        complexSrc = new complex[16];
        for (int i=0; i<SIZE*SIZE; i++) {
            complexSrc[i].real = areal[i];
            complexSrc[i].image = aimage[i];
            double real = complexSrc[i].real;
            double image = complexSrc[i].image;
            if (image > 0.0) {
                printf("raw %f+%fi\n", real, image);
            } else {
                printf("raw %f%fi\n", real, image);
            }
        }
    }
}

//complex *aaa = new complex[16];
//complex *bbb = new complex[16];
//
//
//
void shiyangtestifft()
{
//    for (int i=0; i<16; i++) {
//        aaa[i].real = areal[i];
//        aaa[i].image = aimage[i];
//    }
//
//    Ifft();
}

void shiyangtestcvidft2d()
{
    cv::Mat complexI;
    cv::Mat planes[] = {cv::Mat::zeros(cvSize(SIZE,SIZE), CV_32F), cv::Mat::zeros(cvSize(SIZE,SIZE), CV_32F) };
    for (int j=0; j<SIZE; j++) {
        for (int i=0; i<SIZE; i++) {
            int idx = j*SIZE+i;
//            int idx = i*SIZE+j;
            planes[0].at<float>(idx) = areal[idx];
            planes[1].at<float>(idx) = aimage[idx];
            float real = planes[0].at<float>(idx);
            float image = planes[1].at<float>(idx);
            if (image > 0.0) {
                printf("cv raw %f+%fi\n", real, image);
            } else {
                printf("cv raw %f%fi\n", real, image);
            }
        }
    }
    
    cv::merge(planes, 2, complexI);
    cv::Mat ifftcv;
    idft(complexI, ifftcv);
    split(ifftcv, planes);
    
    for (int j=0; j<SIZE; j++) {
        for (int i=0; i<SIZE; i++) {
            int idx = j*SIZE+i;
//            int idx = i*SIZE+j;
            float real = planes[0].at<float>(idx);
            float image = planes[1].at<float>(idx);
            if (image > 0.0)
                printf("cv idft %f+%fi\n", real, image);
            else
                printf("cv idft %f%fi\n", real, image);
        }
    }
    
}

#if 0
raw -1.000000+98765.000000i
raw 0.000000-12355.000000i
raw 1.0000000.000000i
raw -10.000000+1.000000i
raw 100.000000+88888.000000i
raw 5556.000000-777777.000000i
raw -8888.000000+9.000000i
raw 10.000000+45444.000000i
raw -99999.000000+66677.000000i
raw 765432.875000+111.777000i
raw 1234567936.000000+888.221985i
raw 888.000000+888.443970i
raw 0.000000+22222.699219i
raw -908765.000000-2348.000000i
raw 876578.000000-384623.000000i
raw 2348.000000+2736.000000i
Output IFFT results:
77200080.000000+53154.492188i
-77168376.000000-32454.142578i
77216896.000000-39757.855469i
-77273560.000000-50080.664062i

-77184552.000000-2778.909668i
77130056.000000-60759.023438i
-77180352.000000+109549.687500i
77243184.000000-54058.757812i

77204208.000000-72526.546875i
-77163464.000000-83684.437500i
76995808.000000+17547.355469i
-77061576.000000+125080.804688i

-77219736.000000+548.214355i
77204888.000000+152203.843750i
-77032352.000000-115118.937500i
77088848.000000-45630.132812i

cv raw -1.000000+98765.000000i
cv raw 0.000000-12355.000000i
cv raw 1.0000000.000000i
cv raw -10.000000+1.000000i
cv raw 100.000000+88888.000000i
cv raw 5556.000000-777777.000000i
cv raw -8888.000000+9.000000i
cv raw 10.000000+45444.000000i
cv raw -99999.000000+66677.000000i
cv raw 765432.875000+111.777000i
cv raw 1234567936.000000+888.221985i
cv raw 888.000000+888.443970i
cv raw 0.000000+22222.699219i
cv raw -908765.000000-2348.000000i
cv raw 876578.000000-384623.000000i
cv raw 2348.000000+2736.000000i
cv idft 1235201280.000000-850471.875000i
cv idft -1234694016.000000+519266.312500i
cv idft 1235470336.000000+636125.687500i
cv idft -1236376960.000000+801290.625000i
cv idft -1234952832.000000+44462.562500i
cv idft 1234080896.000000+972144.375000i
cv idft -1234885632.000000-1752795.000000i
cv idft 1235890944.000000+864940.125000i
cv idft 1235267328.000000+1160424.750000i
cv idft -1234615424.000000+1338951.000000i
cv idft 1231932928.000000-280757.687500i
cv idft -1232985216.000000-2001292.875000i
cv idft -1235515776.000000-8771.437500i
cv idft 1235278208.000000-2435261.500000i
cv idft -1232517632.000000+1841903.000000i
cv idft 1233421568.000000+730082.125000i
#endif
