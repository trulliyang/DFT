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

float *imgSrc;
float *magPlusOne;
float *normalise;
float *logE;
float *highpass;

int *a,*b;
int nLen, init_nLen ,mLen, init_mLen, N, M;
FILE *dataFile;

complex *A,*A_In,*W;
void readData()
{
    int i,j;
//    int data[16][16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
//    dataFile = fopen("dataIn.txt","r");
//    fscanf(dataFile,"%d %d", &init_mLen, &init_nLen);
    init_nLen = 512;
    init_mLen = 512;
    
    
    M = calculate_M(init_mLen);
    N = calculate_M(init_nLen);
    nLen = (int)pow(2,N);
    mLen = (int)pow(2,M);
    A_In = (complex *)malloc(complexsize*nLen*mLen);
    
    for(i=0; i<init_mLen; i++)
    {
        for(j=0; j<init_nLen; j++)
        {
            A_In[i*nLen+j].real = imgSrc[i*nLen+j];
            A_In[i*nLen+j].image = 0.0;
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
            W[p].real = (float)cos(2*PI*p/fft_nLen);
            W[p].image = (float)(-1*sin(2*PI*p/fft_nLen));
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
    
    if (true) {
        doMagPlusOne();
        doLogE();
        doNormalise();
        doHighPass();
    }
    
    if (true) {
        doHighPass(A_In);
    }
    
    Ifft();
    
    char *imgChar = new char[512*512];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            imgChar[idx] = (char) (255*(A_In[idx].real/255.0));
        }
    }
    
    cvNamedWindow("shiyangifft", 1);
    IplImage *ImageIFFT = cvCreateImage(cvSize(512,512), IPL_DEPTH_8U, 1);
    ImageIFFT->imageData = imgChar;
    cvShowImage("shiyangifft", ImageIFFT);
//    cvWaitKey();
    
//    printResult_Ifft();
}

void doMagPlusOne()
{
    magPlusOne = new float[512*512];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            int idx = i*nLen+j;
            float sqroot = sqrt(A_In[idx].real*A_In[idx].real+A_In[idx].image*A_In[idx].image);
            magPlusOne[idx] = 1.0 + sqroot;
        }
    }
}

void doLogE()
{
    logE = new float[512*512];
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
    float max = -1.0;
    float min = 10000.0;
    normalise = new float[512*512];
    char *normaliseChar = new char[512*512];
    
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
    
    float delta = max - min;
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
    IplImage *ImageNMLS = cvCreateImage(cvSize(512,512), IPL_DEPTH_8U, 1);
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
            float dist = sqrt(di*di+dj*dj);
            int idx = i*nLen+j;
//            int idx = i*mLen+j;
//            bool need =dist<512/10;
//            if (need) {
//                _in[idx].real *= 0.0;
//                _in[idx].image *= 0.0;
//            }
            float alpha = dist/512/sqrt(2.0);
            _in[idx].real *= (alpha);
            _in[idx].image *= (alpha);
            
            
        }
    }
}

void doHighPass()
{
    highpass = new float[512*512];
    char *highpassChar = new char[512*512];
    for(int i=0; i<init_mLen; i++)
    {
        for(int j=0; j<init_nLen; j++)
        {
            bool corner0 = (i<=init_mLen/2) && (j<=init_nLen/2);
            bool corner1 = (i<=init_mLen/2) && (j>init_nLen/2);
            bool corner2 = (i>init_mLen/2) && (j<=init_nLen/2);
            bool corner3 = (i>init_mLen/2) && (j>init_nLen/2);
            int di=512,dj=512;
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
            float dist = sqrt(di*di+dj*dj);
            int idx = i*nLen+j;
            //            int idx = i*mLen+j;
            bool need =dist<512/10;
            if (need) {
                highpass[idx] *= 0.0;
            } else {
                highpass[idx] = normalise[idx];
            }
            highpassChar[idx] = (char) (255.0*highpass[idx]);
        }
    }
    
    cvNamedWindow("highpassChar", 1);
    IplImage *ImageHP = cvCreateImage(cvSize(512,512), IPL_DEPTH_8U, 1);
    ImageHP->imageData = highpassChar;
    cvShowImage("highpassChar", ImageHP);
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

void shiyanggetImgSrc()
{
    imgSrc = new float[512*512];
    cv::Mat img = cv::imread("/Users/admin/Desktop/beautyTest/dilireba.png", CV_LOAD_IMAGE_GRAYSCALE);
    for (int j=0; j<512; j++) {
        for (int i=0; i<512; i++) {
            int idx = j*512+i;
            imgSrc[idx] = img.at<unsigned char>(j, i);
            if (idx<20) printf("imgSrc %f\n", imgSrc[idx]);
        }
    }
}
