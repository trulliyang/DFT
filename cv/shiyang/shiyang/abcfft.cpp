//
//  abcfft.cpp
//  shiyang
//
//  Created by admin on 2017/10/18.
//  Copyright © 2017年 shiyang. All rights reserved.
//

#include "abcfft.hpp"

#include"math.h"
#include <complex>

#define PI 3.1415926535897932384626

using namespace std;

complex<double> * DataFitFormat(unsigned char *data,int lWidth,int lHeight)//将数组转换为适合FFT处理的数据（数组长度为2的整数次幂）,填充的数据补零操作.当lHeight=1时表示为对一维数组处理.data为对二维数据的一维表示，是按照从左到右，从上到下。
{
    complex<double> *TD;
    int w=1;
    int h=1;
    int wp=0;//存储w的2的幂数
    int hp=0;//存储h的2的幂数
    //////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
    while(w<lWidth)
    {
        w=w*2;
        wp++;
    }
    while(h<lHeight)
    {
        h=h*2;
        hp++;
    }
    TD=new complex<double>[w*h];
    ////////////////////////////////////////////////////////////////////////////////
    for(int i=0;i<h;i++)
    {
        if(i<lHeight)
        {
            for(int j=0;j<w;j++)
            {
                if(j< lWidth)
                {
                    TD[i*w+j]=complex<double>(data[i*lWidth+j],0);//将char数据，准换为实数为data数据，虚数为0的复数
                }
                else
                {
                    TD[i*w+j]=complex<double>(0,0);//对于超出原数据的数据进行补零操作
                }
                
            }
        }
        else
        {
            for(int j=0;j<w;j++)
            {
                TD[i*w+j]=complex<double>(0,0);//对于超出原数据的数据进行补零操作
                
            }
        }
    }
    
    return TD;
    
}

void InitTDAndFD(complex<double> *&TD,complex<double> *&FD,unsigned char *data,int lWidth,int lHeight)//初始化TD和FD(FFT操作)
{
    int w=1;
    int h=1;
    int wp=0;//存储w的2的幂数
    int hp=0;//存储h的2的幂数
    complex<double> *TmpFD;
    //////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
    while(w<lWidth)
    {
        w=w*2;
        wp++;
    }
    while(h<lHeight)
    {
        h=h*2;
        hp++;
    }
    
    TmpFD=new complex<double>[w*h];
    
    for(int i=0;i<w*h;i++)
    {
        TmpFD[i]=complex<double>(0.0,0.0);//FD初值设为0
    }
    
    
    TD=DataFitFormat(data,lWidth,lHeight);//调用已有函数DataFitFormat初始化TD
    FD=TmpFD;
    
}

void FFT_1D(complex<double> *&TD,complex<double>*&FD,int Len)//一维FFT运算，len为一维数组的真实长度。而TD和FD数组的长度都是经过 InitTDAndFD得到的适合FFT处理的数组长度为2的整数次幂的数组 。
{
    //long i,j,k;
    int l=1;
    int lp=0;
    int p=0;
    double angle=0;//中间变量及角度
    complex<double> *W,*X1,*X2,*X;
    
    while(l<Len)
    {
        l=l*2;
        lp++;
    }
    int r=lp;
    
    long N=1<<r;//快速傅里叶变换点数 2的r次幂;
    W=new complex<double>[N/2];//存放旋转系数
    X1=new complex<double>[N];//
    X2=new complex<double>[N];//分配运算的存储器
    for(long i=0;i<N/2;i++)
    {
        angle=-i*PI*2/N;
        W[i]=complex<double>(cos(angle),sin(angle));
    }
    
    memcpy(X1,TD,sizeof(complex<double>)*N);//将TD所在的内存数据拷贝到X1中
    
    ///////////////////////////核心算法：蝶形运算/////////////
    for(long k=0;k<r;k++)
    {
        for(long j=0;j<(1<<k);k++)
        {
            for(long i=0;i<(1<<(r-k-1));i++)
            {
                p=j*(1<<(r-k));
                X2[i+p]=X1[i+p]+X1[i+p+(int)(1<<(r-k-1))];
                X2[i+p+(int)(1<<(r-k-1))]=(X1[i+p]-X1[i+p+(int)(1<<(r-k-1))]) *W[i*(1<<k)];
            }
        }
        X=X1;
        X1=X2;
        X2=X;
    }
    
    /////////////////重新排序，将反序->正序//////////////
    for(int j=0;j<N;j++)
    {
        p=0;
        for(int i=0;i<r;i++)
        {
            if(j&(1<<i))
            {
                p+=1<<(r-i-1);
            }
        }
        FD[j]=X1[p];
    }
    delete W;
    delete X1;
    delete X2;
    
    
}

void FFT_2D(complex<double>*&TD,complex<double> *&FD,int lWidth,int lHeight)//由一维FFT推算出二维FFT。lWidth，lHeight分别为要处理数据的宽和高。TD数组的长度为2的整数次幂，由专门的函数对原始数据进行处理
{
    
    int w=1;
    int h=1;
    int wp=0;//存储w的2的幂数
    int hp=0;//存储h的2的幂数
    complex<double> *TmpTD,*TmpFD;//存放临时的一列的数据
    
    //////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
    while(w<lWidth)
    {
        w=w*2;
        wp++;
    }
    while(h<lHeight)
    {
        h=h*2;
        hp++;
    }
    ////////////////////////先按y方向进行快速的一维FFT运算 //////////////////////////////////
    TmpTD=new complex<double>[h];
    TmpFD=new complex<double>[h];
    
    for(int i=0;i<w;i++)
    {
        //先按y方向进行快速的一维FFT运算
        for(int j=0;j<h;j++)
        {
            TmpTD[j]=TD[j*w+i];
        }
        FFT_1D(TmpTD,TmpFD,lHeight);
        //保存结果
        for(int j=0;j<h;j++)
        {
            TD[j*w+i]=TmpFD[j];
        }
    }
    delete TmpTD;
    delete TmpFD;
    ///////////////////////再按x方向进行快速的一维FFT运算///////////////////////////
    TmpTD=new complex<double>[h];
    TmpFD=new complex<double>[h];
    
    for(int i=0;i<h;i++)
    {
        ////再按x方向进行快速的一维FFT运算
        for(int j=0;j<w;j++)
        {
            TmpTD[j]=TD[i*w+j];
        }
        FFT_1D(TmpTD,TmpFD,lWidth);
        //保存结果
        for(int j=0;j<w;j++)
        {
            FD[i*w+j]=TmpFD[j];
        }
    }
    delete TmpTD;
    delete TmpFD;
    
}

void abcStart() {
//    InitTDAndFD(complex<double> *&TD,complex<double> *&FD,unsigned char *data,int lWidth,int lHeight)
//    void FFT_2D(complex<double>*&TD,complex<double> *&FD,int lWidth,int lHeight)
    int w=16;
    int h=16;
    complex<double> *td = new complex<double>[w*h];
    complex<double> *fd = new complex<double>[w*h];
    for (int i=0; i<h; i++) {
        for (int j=0; j<w; j++) {
            int index = i*w+j;
            td[index] = complex<double>(0,0);
            
            fd[index] = complex<double>(0,0);
        }
    }
    FFT_2D(td, fd, w, h);
    delete td;
    delete fd;
}

