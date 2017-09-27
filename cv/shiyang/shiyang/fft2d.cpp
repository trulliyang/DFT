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
    init_nLen = 16;
    init_mLen = 16;
    
    
    M = calculate_M(init_mLen);
    N = calculate_M(init_nLen);
    nLen = (int)pow(2,N);
    mLen = (int)pow(2,M);
    A_In = (complex *)malloc(complexsize*nLen*mLen);
    
    for(i=0; i<init_mLen; i++)
    {
        for(j=0; j<init_nLen; j++)
        {
            A_In[i*nLen+j].real = i*nLen+j;
            A_In[i*nLen+j].image = 0.0;
        }
    }
    fclose(dataFile);
    
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
    printResult_fft();
    Ifft();
    printResult_Ifft();
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
