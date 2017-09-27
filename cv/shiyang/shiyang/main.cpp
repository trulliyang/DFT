#if 0

//#include <stdio.h>
//#include <string.h>
//#include <memory.h>
//#include <opencv2/opencv.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc.hpp>
//
//#include "color_adjust.hpp"
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#define eps 1e-6
#define PI 3.14159265354
typedef std::complex<double> complex_t;
using namespace std;

#define NUM 64

//旋转因子的计算
complex_t W(int k,int n,int N){
    return complex_t(cos(2*PI*k*n/N),-sin(2*PI*k*n/N));
}

//格式化 零
complex_t format(complex_t &c){

    if(fabs(c.real())<eps) c.real(0);
    if(fabs(c.imag())<eps) c.imag(0);
    return c;
}

double format(double &c){
    if(fabs(c)<eps) c=0;
    return c;
}

//离散序列的DFT计算,只针对实数序列 ,完全按照DFT的公式来计算,O(n^2)的复杂度
void DFT(vector<double> &x_n,vector<complex_t> &X_k){
    X_k.clear();
    int N=x_n.size();
    for(int k=0;k<N;++k){
        complex_t t(0,0);
        for(int n=0;n<N;++n){
            t+=x_n[n]*W(k,n,N);
        }
        X_k.push_back(format(t));
    }
    
    for(int i=0;i<N;++i){
        cout<<format(X_k[i])<<endl;
    }
    
}

//IDFT的计算,只针对实数序列
void IDFT(vector<complex_t> &X_k,vector<double> &x_n){
    x_n.clear();
    int N=X_k.size();
    for(int n=0;n<N;++n){
        complex_t t(0,0);
        for(int k=0;k<N;++k){
            t+=X_k[k]*W(k,-n,N);
        }
        x_n.push_back(t.real()/N);//运算结果只剩下实部
        //cout<<(t/(double)N)<<endl;
    }
    
    for(int i=0;i<N;++i){
        cout<<format(x_n[i])<<endl;
    }
    
}

void DFT_test(){
    int N=NUM;
    vector<double> x_n(N);
    vector<complex_t> X_k(N);
    for(int i=0;i<N;++i){
        x_n[i]=i;
    }
    DFT(x_n,X_k);
    IDFT(X_k,x_n);
}


//保证N是2的n次幂
int bitlen(int N){
    int n=0;
    while((N&1)==0){
        n++;
        N>>=1;
    }
    return n;
}


int reverse_bit(int n,int len){//bit反转
    int tmp=0;
    while(len--){
        tmp+=((n&1)<<len);
        n>>=1;
    }
    return tmp;
    
}

//序数重排
void resort(vector<complex_t> &x_n,int N){
    vector<complex_t> v(x_n);
    int len=bitlen(N);
    for(int i=0;i<N;++i){
        x_n[i]=v[reverse_bit(i,len)];
    }
}


//基2,FFT算法实现,O(nlogn)的复杂度
void FFT(vector<complex_t> &x_n){
    int N=x_n.size();
    int r=bitlen(N);
    
    vector<complex_t> W(N);
    
    //预先计算旋转因子
    for(int i=0;i<N;++i){
        double angle=-i*2*PI/N;
        W[i]=complex_t(cos(angle),sin(angle));
    }
    
    
    for(int k=0;k<r;++k){//迭代次数
        for(int j=0;j<(1<<k);++j){
            int butterfly=1<<(r-k);
            int p=j*butterfly;
            int s=p+butterfly/2;
            for(int i=0;i<butterfly/2;++i){
                complex_t c=x_n[i+p]+x_n[i+s];
                x_n[i+s]=(x_n[i+p]-x_n[i+s])*W[i*(1<<k)];
                x_n[i+p]=c;
            }
        }
    }
    
    //次序重排
    resort(x_n,N);
    for(int i=0;i<N;++i){
        cout<<format(x_n[i])<<endl;
    }
    
}

//IFFT,与FFT基本一致
void IFFT(vector<complex_t> &x_n){
    int N=x_n.size();
    int r=bitlen(N);
    
    vector<complex_t> W(N);
    
    //预先计算旋转因子
    for(int i=0;i<N;++i){
        double angle=i*2*PI/N;//IFFT的旋转因子与FFT的旋转因子差一个负号
        W[i]=complex_t(cos(angle),sin(angle));
    }
    
    
    for(int k=0;k<r;++k){//迭代次数
        for(int j=0;j<(1<<k);++j){
            int butterfly=1<<(r-k);
            int p=j*butterfly;
            int s=p+butterfly/2;
            for(int i=0;i<butterfly/2;++i){
                complex_t c=x_n[i+p]+x_n[i+s];
                x_n[i+s]=(x_n[i+p]-x_n[i+s])*W[i*(1<<k)];
                x_n[i+p]=c;
            }
        }
    }
    
    //次序重排
    resort(x_n,N);
    for(int i=0;i<N;++i){
        x_n[i]/=N;//IFFT与FFT还差一个系数
        cout<<format(x_n[i])<<endl;
    }
}

void FFT_test(){
    int N=NUM;
    vector<complex_t> x_n;
    complex_t c(0,0);
    for(int i=0;i<N;++i){
        c.real(i);
        x_n.push_back(c);
    }
    
    FFT(x_n);
    IFFT(x_n);
}


int main(){
    DFT_test();
    cout<<endl;
//    FFT_test();
    return 0;
}
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    int ret = 0;
//    int nWidth, nHeight, Stride;
//    const char* path = "/Users/admin/Desktop/beautyTest/";
//    char Image_Path[256], Image_outPath[256];// = (const char*)strcpy("%s1.jpg", path);
//    //char *Image_outPath = strcpy("%s1_adjustColorOut.jpg", path);
//    sprintf(Image_Path, "%s7.png", path);
//    sprintf(Image_outPath, "%s7_adjustColorOut.png", path);
//    IplImage *Image = cvLoadImage(Image_Path, CV_LOAD_IMAGE_ANYCOLOR);
//    if(NULL == Image->imageData)
//    {
//        printf("%s is not avalibal!", Image_Path);
//    }
//    unsigned char *pSrcImage = (unsigned char *)(Image->imageData);
//    
//    nWidth = Image->width;
//    Stride = Image->widthStep;
//    nHeight = Image->height;
//    printf("w=%d,h=%d,std=%d\n", nWidth, nHeight, Stride);
//    unsigned char *pResultImage = (unsigned char*)malloc(nHeight * Stride);
//    if(NULL == pResultImage)
//    {
//        printf("malloc error!");
//    }
//    memset(pResultImage, 0, nHeight * Stride);
//    int average_brightness[3] = {190, 145, 125};
//    ret = color_adjust_RGB(pSrcImage, pResultImage, nWidth, nHeight, Stride, average_brightness);
//    if(0 != ret)
//    {
//        printf("color_adjust is error, error code is %d!", ret); //error code is %d”，ret)；
//    }
//    //memcpy(Image->imageData, pResultImage, nHeight * Stride);
//    cvNamedWindow("colorAdjust_resualt", 0);
//    cvShowImage("colorAdjust_resualt", Image);
//    cvSaveImage(Image_outPath, Image);
//    cvWaitKey();
//    cvDestroyWindow("colorAdjust_resualt");
//    cvReleaseImage(&Image);
//    return ret;
//}
#elif 0
#elif 0

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define SIZE 1000
#define VALUE_MAX 2000
struct Complex_{
    double real;
    double imagin;
};
typedef struct Complex_ Complex;

void setInput(double * data,int n){
    printf("Setinput signal:\n");
    srand((int)time(0));
    for(int i=0;i<SIZE;i++){
        data[i]=rand()%VALUE_MAX;
        printf("%lf\n",data[i]);
    }
    
}
void DFT(double * src,Complex * dst,int size){
    clock_t start,end;
    start=clock();
    
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n;
            real+=src[n]*cos(x/size);
            imagin+=src[n]*(-sin(x/size));
            
        }
        dst[m].imagin=imagin;
        dst[m].real=real;
        if(imagin>=0.0)
            printf("%lf+%lfj\n",real,imagin);
        else
            printf("%lf%lfj\n",real,imagin);
    }
    end=clock();
    printf("DFT use time :%lf for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size);
    
}

void IDFT(Complex *src,Complex *dst,int size){
    //Complex temp[SIZE];
    clock_t start,end;
    start=clock();
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n/size;
            real+=src[n].real*cos(x)-src[n].imagin*sin(x);
            imagin+=src[n].real*sin(x)+src[n].imagin*cos(x);
            
        }
        real/=SIZE;
        imagin/=SIZE;
        if(dst!=NULL){
            dst[m].real=real;
            dst[m].imagin=imagin;
        }
        if(imagin>=0.0)
            printf("%lf+%lfj\n",real,imagin);
        else
            printf("%lf%lfj\n",real,imagin);
    }
    end=clock();
    printf("IDFT use time :%lfs for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size);
    
    
    
}
int main(int argc, const char * argv[]) {
    double input[SIZE];
    Complex dst[SIZE];
    setInput(input,SIZE);
    printf("\n\n");
    DFT(input, dst, SIZE);
    printf("\n\n");
    IDFT(dst, NULL, SIZE);
    
}
#elif 0
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define VALUE_MAX 256
#define WIDTH 128
#define HEIGHT 128

struct Complex_{
    double real;
    double imagin;
};
typedef struct Complex_ Complex;


int Initdata(double (*src)[WIDTH],int size_w,int size_h){
    srand((int)time(0));
    for(int i=0;i<size_w;i++){
        for(int j=0;j<size_h;j++){
            src[i][j]=rand()%VALUE_MAX;
            printf("%lf ",src[i][j]);
        }
        printf(";\n");
    }
    return 0;
}
int DFT2D(double (*src)[WIDTH],Complex (*dst)[WIDTH],int size_w,int size_h){
    for(int u=0;u<size_w;u++){
        for(int v=0;v<size_h;v++){
            double real=0.0;
            double imagin=0.0;
            for(int i=0;i<size_w;i++){
                for(int j=0;j<size_h;j++){
                    double I=src[i][j];
                    double x=M_PI*2*((double)i*u/(double)size_w+(double)j*v/(double)size_h);
                    real+=cos(x)*I;
                    imagin+=-sin(x)*I;
                    
                }
            }
            dst[u][v].real=real;
            dst[u][v].imagin=imagin;
            if(imagin>=0)
                printf("%lf+%lfj ",real,imagin);
            else
                printf("%lf%lfj ",real,imagin);
        }
        printf(";\n");
    }
    return 0;
}
int IDFT2D(Complex (*src)[WIDTH],Complex (*dst)[WIDTH],int size_w,int size_h){
    for(int i=0;i<size_w;i++){
        for(int j=0;j<size_h;j++){
            double real=0.0;
            double imagin=0.0;
            for(int u=0;u<size_w;u++){
                for(int v=0;v<size_h;v++){
                    double R=src[u][v].real;
                    double I=src[u][v].imagin;
                    double x=M_PI*2*((double)i*u/(double)size_w+(double)j*v/(double)size_h);
                    real+=R*cos(x)-I*sin(x);
                    imagin+=I*cos(x)+R*sin(x);
                    
                }
            }
            dst[i][j].real=(1./(size_w*size_h))*real;
            dst[i][j].imagin=(1./(size_w*size_h))*imagin;
            if(imagin>=0)
                printf("%lf+%lfj ",dst[i][j].real,dst[i][j].imagin);
            else
                printf("%lf%lfj ",dst[i][j].real,dst[i][j].imagin);
        }
        printf(";\n");
    }
    return 0;
}


int main() {
    double src[WIDTH][HEIGHT];
    Complex dst[WIDTH][HEIGHT];
    Complex dst_[WIDTH][HEIGHT];
    Initdata(src, WIDTH, HEIGHT);
    printf("\n\n");
    DFT2D(src,dst,WIDTH,HEIGHT);
    printf("\n\n");
    IDFT2D(dst,dst_,WIDTH,HEIGHT);
}

#elif 0
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//mac下M_PI在math.h中有宏定义，所以这里我们选择行的宏定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIZE 1024*2048
#define VALUE_MAX 1000

////////////////////////////////////////////////////////////////////
//定义一个复数结构体
///////////////////////////////////////////////////////////////////
struct Complex_{
    double real;
    double imagin;
};
typedef struct Complex_ Complex;
////////////////////////////////////////////////////////////////////
//定义一个复数计算，包括乘法，加法，减法
///////////////////////////////////////////////////////////////////
void Add_Complex(Complex * src1,Complex *src2,Complex *dst){
    dst->imagin=src1->imagin+src2->imagin;
    dst->real=src1->real+src2->real;
}
void Sub_Complex(Complex * src1,Complex *src2,Complex *dst){
    dst->imagin=src1->imagin-src2->imagin;
    dst->real=src1->real-src2->real;
}
void Multy_Complex(Complex * src1,Complex *src2,Complex *dst){
    double r1=0.0,r2=0.0;
    double i1=0.0,i2=0.0;
    r1=src1->real;
    r2=src2->real;
    i1=src1->imagin;
    i2=src2->imagin;
    dst->imagin=r1*i2+r2*i1;
    dst->real=r1*r2-i1*i2;
}
////////////////////////////////////////////////////////////////////
//在FFT中有一个WN的n次方项，在迭代中会不断用到，具体见算法说明
///////////////////////////////////////////////////////////////////
void getWN(double n,double size_n,Complex * dst){
    double x=2.0*M_PI*n/size_n;
    dst->imagin=-sin(x);
    dst->real=cos(x);
}
////////////////////////////////////////////////////////////////////
//随机生成一个输入，显示数据部分已经注释掉了
//注释掉的显示部分为数据显示，可以观察结果
///////////////////////////////////////////////////////////////////
void setInput(double * data,int  n){
    //printf("Setinput signal:\n");
    srand((int)time(0));
    for(int i=0;i<SIZE;i++){
        data[i]=rand()%VALUE_MAX;
        //printf("%lf\n",data[i]);
    }
    
}
////////////////////////////////////////////////////////////////////
//定义DFT函数，其原理为简单的DFT定义，时间复杂度O（n^2）,
//下面函数中有两层循环，每层循环的step为1，size为n，故为O（n*n）,
//注释掉的显示部分为数据显示，可以观察结果
///////////////////////////////////////////////////////////////////
void DFT(double * src,Complex * dst,int size){
    clock_t start,end;
    start=clock();
    
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n;
            real+=src[n]*cos(x/size);
            imagin+=src[n]*(-sin(x/size));
            
        }
        dst[m].imagin=imagin;
        dst[m].real=real;
        /* if(imagin>=0.0)
         printf("%lf+%lfj\n",real,imagin);
         else
         printf("%lf%lfj\n",real,imagin);*/
    }
    end=clock();
    printf("DFT use time :%lf for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size);
    
}
////////////////////////////////////////////////////////////////////
//定义IDFT函数，其原理为简单的IDFT定义，时间复杂度O（n^2）,
//下面函数中有两层循环，每层循环的step为1，size为n，故为O（n*n）,
///////////////////////////////////////////////////////////////////
void IDFT(Complex *src,Complex *dst,int size){
    clock_t start,end;
    start=clock();
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n/size;
            real+=src[n].real*cos(x)-src[n].imagin*sin(x);
            imagin+=src[n].real*sin(x)+src[n].imagin*cos(x);
            
        }
        real/=SIZE;
        imagin/=SIZE;
        if(dst!=NULL){
            dst[m].real=real;
            dst[m].imagin=imagin;
        }
        if(imagin>=0.0)
            printf("%lf+%lfj\n",real,imagin);
        else
            printf("%lf%lfj\n",real,imagin);
    }
    end=clock();
    printf("IDFT use time :%lfs for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size);
    
    
}
////////////////////////////////////////////////////////////////////
//定义FFT的初始化数据，因为FFT的数据经过重新映射，递归结构
///////////////////////////////////////////////////////////////////
int FFT_remap(double * src,int size_n){
    
    if(size_n==1)
        return 0;
    double * temp=(double *)malloc(sizeof(double)*size_n);
    for(int i=0;i<size_n;i++)
        if(i%2==0)
            temp[i/2]=src[i];
        else
            temp[(size_n+i)/2]=src[i];
    for(int i=0;i<size_n;i++)
        src[i]=temp[i];
    free(temp);
    FFT_remap(src, size_n/2);
    FFT_remap(src+size_n/2, size_n/2);
    return 1;
    
    
}
////////////////////////////////////////////////////////////////////
//定义FFT，具体见算法说明，注释掉的显示部分为数据显示，可以观察结果
///////////////////////////////////////////////////////////////////
void FFT(double * src,Complex * dst,int size_n){
    
    FFT_remap(src, size_n);
    // for(int i=0;i<size_n;i++)
    //    printf("%lf\n",src[i]);
    clock_t start,end;
    start=clock();
    int k=size_n;
    int z=0;
    while (k/=2) {
        z++;
    }
    k=z;
    if(size_n!=(1<<k))
        exit(0);
    Complex * src_com=(Complex*)malloc(sizeof(Complex)*size_n);
    if(src_com==NULL)
        exit(0);
    for(int i=0;i<size_n;i++){
        src_com[i].real=src[i];
        src_com[i].imagin=0;
    }
    for(int i=0;i<k;i++){
        z=0;
        for(int j=0;j<size_n;j++){
            if((j/(1<<i))%2==1){
                Complex wn;
                getWN(z, size_n, &wn);
                Multy_Complex(&src_com[j], &wn,&src_com[j]);
                z+=1<<(k-i-1);
                Complex temp;
                int neighbour=j-(1<<(i));
                temp.real=src_com[neighbour].real;
                temp.imagin=src_com[neighbour].imagin;
                Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
                Sub_Complex(&temp, &src_com[j], &src_com[j]);
            }
            else
                z=0;
        }
        
    }
    
    /* for(int i=0;i<size_n;i++)
     if(src_com[i].imagin>=0.0){
     printf("%lf+%lfj\n",src_com[i].real,src_com[i].imagin);
     }
     else
     printf("%lf%lfj\n",src_com[i].real,src_com[i].imagin);*/
    for(int i=0;i<size_n;i++){
        dst[i].imagin=src_com[i].imagin;
        dst[i].real=src_com[i].real;
    }
    end=clock();
    printf("FFT use time :%lfs for Datasize of:%d\n",(double)(end-start)/CLOCKS_PER_SEC,size_n);
    
}

////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
    double *input = new double[SIZE];
    Complex *dst = new Complex[SIZE];
    setInput(input, SIZE);
    printf("\n\n");
//    DFT(input, dst, SIZE);
    printf("\n\n");
    FFT(input, dst, SIZE);
    //IDFT(dst, NULL, SIZE);
    getchar();
}
#elif 1
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//mac下M_PI在math.h中有宏定义，所以这里我们选择行的宏定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define WIDTH 512
#define HEIGHT 512
#define SIZE (WIDTH*HEIGHT)
#define VALUE_MAX 256

////////////////////////////////////////////////////////////////////
//定义一个复数结构体
///////////////////////////////////////////////////////////////////
struct Complex_{
    double real;
    double imagin;
};
typedef struct Complex_ Complex;

int isBase2(int size_n){
    int k=size_n;
    int z=0;
    while (k/=2) {
        z++;
    }
    k=z;
    if(size_n!=(1<<k))
        return -1;
    else
        return k;
}
////////////////////////////////////////////////////////////////////
//复数基本运算
///////////////////////////////////////////////////////////////////
void Add_Complex(Complex * src1,Complex *src2,Complex *dst){
    dst->imagin=src1->imagin+src2->imagin;
    dst->real=src1->real+src2->real;
}
void Sub_Complex(Complex * src1,Complex *src2,Complex *dst){
    dst->imagin=src1->imagin-src2->imagin;
    dst->real=src1->real-src2->real;
}
void Multy_Complex(Complex * src1,Complex *src2,Complex *dst){
    double r1=0.0,r2=0.0;
    double i1=0.0,i2=0.0;
    r1=src1->real;
    r2=src2->real;
    i1=src1->imagin;
    i2=src2->imagin;
    dst->imagin=r1*i2+r2*i1;
    dst->real=r1*r2-i1*i2;
}
void Copy_Complex(Complex * src,Complex *dst){
    dst->real=src->real;
    dst->imagin=src->imagin;
}
void Show_Complex(Complex * src,int size_n){
    if(size_n==1){
        if(src->imagin>=0.0)
            printf("%lf+%lfj  ",src->real,src->imagin);
        else
            printf("%lf%lfj  ",src->real,src->imagin);
        
    }
    else if(size_n>1){
        for(int i=0;i<size_n;i++)
            if(src[i].imagin>=0.0){
                printf("%lf+%lfj  ",src[i].real,src[i].imagin);
            }
            else
                printf("%lf%lfj  ",src[i].real,src[i].imagin);
        
        
        
    }
    
    
}
////////////////////////////////////////////////////////////////////
//计算WN
///////////////////////////////////////////////////////////////////
void getWN(double n,double size_n,Complex * dst){
    double x=2.0*M_PI*n/size_n;
    dst->imagin=-sin(x);
    dst->real=cos(x);
}
////////////////////////////////////////////////////////////////////
//随机初始化输入
///////////////////////////////////////////////////////////////////
void setInput(double * data,int  n){
    printf("Setinput signal:\n");
    srand((int)time(0));
    for(int i=0;i<n;i++){
//        data[i]=rand()%VALUE_MAX;
        if (i<n/3) {
            data[i]=255;
        } else if (i<n*2/3) {
            data[i]=127;
        } else {
            data[i]=0;
        }
        
    }
    
}
////////////////////////////////////////////////////////////////////
//标准DFT
///////////////////////////////////////////////////////////////////
void DFT(double * src,Complex * dst,int size){
    
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n;
            real+=src[n]*cos(x/size);
            imagin+=src[n]*(-sin(x/size));
            
        }
        dst[m].imagin=imagin;
        dst[m].real=real;
        
    }
}
////////////////////////////////////////////////////////////////////
//IDT，复原傅里叶
///////////////////////////////////////////////////////////////////
void IDFT(Complex *src,Complex *dst,int size){
    for(int m=0;m<size;m++){
        double real=0.0;
        double imagin=0.0;
        for(int n=0;n<size;n++){
            double x=M_PI*2*m*n/size;
            real+=src[n].real*cos(x)-src[n].imagin*sin(x);
            imagin+=src[n].real*sin(x)+src[n].imagin*cos(x);
            
        }
        real/=SIZE;
        imagin/=SIZE;
        if(dst!=NULL){
            dst[m].real=real;
            dst[m].imagin=imagin;
        }
    }
    
    
}
////////////////////////////////////////////////////////////////////
//FFT前，对输入数据重新排序
///////////////////////////////////////////////////////////////////
int FFTReal_remap(double * src,int size_n){
    
    if(size_n==1)
        return 0;
    double * temp=(double *)malloc(sizeof(double)*size_n);
    for(int i=0;i<size_n;i++)
        if(i%2==0)
            temp[i/2]=src[i];
        else
            temp[(size_n+i)/2]=src[i];
    for(int i=0;i<size_n;i++)
        src[i]=temp[i];
    free(temp);
    FFTReal_remap(src, size_n/2);
    FFTReal_remap(src+size_n/2, size_n/2);
    return 1;
    
    
}
int FFTComplex_remap(Complex * src,int size_n){
    
    if(size_n==1)
        return 0;
    Complex * temp=(Complex *)malloc(sizeof(Complex)*size_n);
    for(int i=0;i<size_n;i++)
        if(i%2==0)
            Copy_Complex(&src[i],&(temp[i/2]));
        else
            Copy_Complex(&(src[i]),&(temp[(size_n+i)/2]));
    for(int i=0;i<size_n;i++)
        Copy_Complex(&(temp[i]),&(src[i]));
    free(temp);
    FFTComplex_remap(src, size_n/2);
    FFTComplex_remap(src+size_n/2, size_n/2);
    return 1;
    
    
}
////////////////////////////////////////////////////////////////////
//FFT公式
///////////////////////////////////////////////////////////////////
void FFT(Complex * src,Complex * dst,int size_n){
    
    int k=size_n;
    int z=0;
    while (k/=2) {
        z++;
    }
    k=z;
    if(size_n!=(1<<k))
        exit(0);
    Complex * src_com=(Complex*)malloc(sizeof(Complex)*size_n);
    if(src_com==NULL)
        exit(0);
    for(int i=0;i<size_n;i++){
        Copy_Complex(&src[i], &src_com[i]);
    }
    FFTComplex_remap(src_com, size_n);
    for(int i=0;i<k;i++){
        z=0;
        for(int j=0;j<size_n;j++){
            if((j/(1<<i))%2==1){
                Complex wn;
                getWN(z, size_n, &wn);
                Multy_Complex(&src_com[j], &wn,&src_com[j]);
                z+=1<<(k-i-1);
                Complex temp;
                int neighbour=j-(1<<(i));
                temp.real=src_com[neighbour].real;
                temp.imagin=src_com[neighbour].imagin;
                Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
                Sub_Complex(&temp, &src_com[j], &src_com[j]);
            }
            else
                z=0;
        }
        
    }
    
    
    for(int i=0;i<size_n;i++){
        Copy_Complex(&src_com[i], &dst[i]);
    }
    free(src_com);
    
    
}
void RealFFT(double * src,Complex * dst,int size_n){
    
    
    int k=size_n;
    int z=0;
    while (k/=2) {
        z++;
    }
    k=z;
    if(size_n!=(1<<k))
        exit(0);
    Complex * src_com=(Complex*)malloc(sizeof(Complex)*size_n);
    if(src_com==NULL)
        exit(0);
    for(int i=0;i<size_n;i++){
        src_com[i].real=src[i];
        src_com[i].imagin=0;
    }
    FFTComplex_remap(src_com, size_n);
    for(int i=0;i<k;i++){
        z=0;
        for(int j=0;j<size_n;j++){
            if((j/(1<<i))%2==1){
                Complex wn;
                getWN(z, size_n, &wn);
                Multy_Complex(&src_com[j], &wn,&src_com[j]);
                z+=1<<(k-i-1);
                Complex temp;
                int neighbour=j-(1<<(i));
                temp.real=src_com[neighbour].real;
                temp.imagin=src_com[neighbour].imagin;
                Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
                Sub_Complex(&temp, &src_com[j], &src_com[j]);
            }
            else
                z=0;
        }
        
    }
    
    
    for(int i=0;i<size_n;i++){
        Copy_Complex(&src_com[i], &dst[i]);
    }
    free(src_com);
    
}

////////////////////////////////////////////////////////////////////
//IFFT实现
////////////////////////////////////////////////////////////////////
void IFFT(Complex * src,Complex * dst,int size_n) {
    for(int i=0;i<size_n;i++)
        src[i].imagin=-src[i].imagin;
    FFTComplex_remap(src, size_n);
    int z,k;
    if((z=isBase2(size_n))!=-1)
        k=isBase2(size_n);
    else
        exit(0);
    for(int i=0;i<k;i++){
        z=0;
        for(int j=0;j<size_n;j++){
            if((j/(1<<i))%2==1){
                Complex wn;
                getWN(z, size_n, &wn);
                Multy_Complex(&src[j], &wn, &src[j]);
                z+=1<<(k-i-1);
                Complex temp;
                int neighbour=j-(1<<(i));
                Copy_Complex(&src[neighbour], &temp);
                Add_Complex(&temp, &src[j], &src[neighbour]);
                Sub_Complex(&temp, &src[j], &src[j]);
            }
            else
                z=0;
        }
        
    }
    for(int i=0;i<size_n;i++){
        dst[i].imagin=(1./size_n)*src[i].imagin;
        dst[i].real=(1./size_n)*src[i].real;
    }
    
    
    
    
}

int DFT2D(double *src,Complex *dst,int size_w,int size_h){
    for(int u=0;u<size_w;u++){
        for(int v=0;v<size_h;v++){
            double real=0.0;
            double imagin=0.0;
            for(int i=0;i<size_w;i++){
                for(int j=0;j<size_h;j++){
                    double I=src[i*size_w+j];
                    double x=M_PI*2*((double)i*u/(double)size_w+(double)j*v/(double)size_h);
                    real+=cos(x)*I;
                    imagin+=-sin(x)*I;
                    
                }
            }
            dst[u*size_w+v].real=real;
            dst[u*size_w+v].imagin=imagin;
            
        }
        
    }
    return 0;
}
/*
 
 */
int IDFT2D(Complex *src,Complex *dst,int size_w,int size_h){
    for(int i=0;i<size_w;i++){
        for(int j=0;j<size_h;j++){
            double real=0.0;
            double imagin=0.0;
            for(int u=0;u<size_w;u++){
                for(int v=0;v<size_h;v++){
                    double R=src[u*size_w+v].real;
                    double I=src[u*size_w+v].imagin;
                    double x=M_PI*2*((double)i*u/(double)size_w+(double)j*v/(double)size_h);
                    real+=R*cos(x)-I*sin(x);
                    imagin+=I*cos(x)+R*sin(x);
                    
                }
            }
            dst[i*size_w+j].real=(1./(size_w*size_h))*real;
            dst[i*size_w+j].imagin=(1./(size_w*size_h))*imagin;
            
        }
    }
    return 0;
}
/*
 
 
 
 */
void ColumnVector(Complex * src,Complex * dst,int size_w,int size_h){
    for(int i=0;i<size_h;i++)
        Copy_Complex(&src[size_w*i], &dst[i]);
    
}
/*
 
 */
void IColumnVector(Complex * src,Complex * dst,int size_w,int size_h){
    for(int i=0;i<size_h;i++)
        Copy_Complex(&src[i],&dst[size_w*i]);
    
}
/*
 */
int FFT2D(double *src,Complex *dst,int size_w,int size_h){
    if(isBase2(size_w)==-1||isBase2(size_h)==-1)
        exit(0);
    Complex *temp=(Complex *)malloc(sizeof(Complex)*size_h*size_w);
    if(temp==NULL)
        return -1;
    
    for(int i=0;i<size_h;i++){
        RealFFT(&src[size_w*i], &temp[size_w*i], size_w);
    }
    
    Complex *Column=(Complex *)malloc(sizeof(Complex)*size_h);
    if(Column==NULL)
        return -1;
    for(int i=0;i<size_w;i++){
        ColumnVector(&temp[i], Column, size_w, size_h);
        FFT(Column, Column, size_h);
        IColumnVector(Column, &temp[i], size_w, size_h);
        
    }
    
    
    
    for(int i=0;i<size_h*size_w;i++)
        Copy_Complex(&temp[i], &dst[i]);
//    Show_Complex(dst, size_h*size_w);
    free(temp);
    free(Column);
    return 0;
}
/*
 */
int IFFT2D(Complex *src,Complex *dst,int size_w,int size_h){
    
    if(isBase2(size_w)==-1||isBase2(size_h)==-1)
        exit(0);
    
    Complex *temp=(Complex *)malloc(sizeof(Complex)*size_h*size_w);
    if(temp==NULL)
        return -1;
    Complex *Column=(Complex *)malloc(sizeof(Complex)*size_h);
    if(Column==NULL)
        return -1;
    
    for(int i=0;i<size_w;i++){
        ColumnVector(&src[i], Column, size_w, size_h);
        IFFT(Column, Column, size_h);
        IColumnVector(Column, &src[i], size_w, size_h);
        
    }
    for(int i=0;i<size_w*size_h;i++)
        src[i].imagin=-src[i].imagin;
    for(int i=0;i<size_h;i++){
        IFFT(&src[size_w*i], &temp[size_w*i], size_w);
        
    }
    
    
    for(int i=0;i<size_h*size_w;i++)
        Copy_Complex(&temp[i], &dst[i]);
//    Show_Complex(dst, size_h*size_w);
    free(temp);
    free(Column);
    return 0;
    
}

void generateImg(Complex *src, int size_w, int size_h, char *dst)
{
    for (int i=0; i<size_w*size_h; i++) {
        dst[i] = src[i].real;
    }
}
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

#include "fft2d.h"

int main(int argc, const char * argv[]) {
    getImageSrc();
    shiyangfft2d();
    shiyangcvdft2d();
    return 1;
    
    IplImage *ImagePic = cvLoadImage("/Users/admin/Desktop/beautyTest/dilireba.png",
                                     CV_LOAD_IMAGE_GRAYSCALE);
    cvNamedWindow("pic", 1);
    cvShowImage("pic", ImagePic);
    if (ImagePic->width == WIDTH && ImagePic->height == HEIGHT) {
        printf("size is matched\n");
    } else {
        printf("size is not matched\n");
    }
    
    double *input = new double[SIZE];
    Complex *dst = new Complex[SIZE];
//    memset(input, 255.0, SIZE*sizeof(double));
//    memset(input+SIZE/4*0, 0, SIZE/4);
//    memset(input+SIZE/4*1, 63, SIZE/4);
//    memset(input+SIZE/4*2, 127, SIZE/4);
//    memset(input+SIZE/4*3, 255, SIZE/4);
//    setInput(input, SIZE);
    for (int i=0; i<SIZE; i++) {
        input[i] = ImagePic->imageData[i];
        if (i<20) {
            printf("input[%d]=%u\n", i, (unsigned char)ImagePic->imageData[i]);
        }
    }
    
    
    cvNamedWindow("ori", 0);
    IplImage *ImageOri = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    for (int i=0; i<SIZE; i++) {
        ImageOri->imageData[i] = (char)input[i];
    }
    cvShowImage("ori", ImageOri);
    
    
    
    
    for (int i=0; i<SIZE; i++) {
//        printf("%lf            ", input[i]);
    }
    printf("\n\n");
    printf("\n\n");
    //    DFT(input, dst, SIZE);
    printf("\n\n");
    FFT2D(input, dst, WIDTH, HEIGHT);
    printf("\n\n");
    
    cvNamedWindow("fft", 0);
    IplImage *Image = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    for (int i=0; i<SIZE; i++) {
        Image->imageData[i] = (char)dst[i].real;
        if (i<256) {
            printf("fft2d =%d, real=%lf,imag=%lf\n", i, dst[i].real, dst[i].imagin);
        }
    }
    cvShowImage("fft", Image);
//    cv::imshow("lalala0",Image->imageData);
    
    double *mag = new double[SIZE];
    for (int i=0; i<SIZE; i++) {
        mag[i] = sqrt(dst[i].real*dst[i].real + dst[i].imagin*dst[i].imagin);
//        mag[i] = sqrt(dst[i].real*dst[i].real + dst[i].real*dst[i].real);
        if (i<20) {
            printf("mag[%d]=%lf\n", i, mag[i]);
        }
    }
    
    double *mag2 = new double[SIZE];
    for (int i=0; i<SIZE; i++) {
        mag2[i] = mag[i]+1.0;
        if (i<20) {
            printf("mag2[%d]=%lf\n", i, mag2[i]);
        }
    }

    double logMin = 100.0;
    double logMax = 0.0;
    double *logE = new double[SIZE];
    for (int i=0; i<SIZE; i++) {
        logE[i] = log(mag2[i]);
        if (logE[i] > logMax) {
            logMax = logE[i];
        }
        if (logE[i] < logMin) {
            logMin = logE[i];
        }
        if (i<20) {
            printf("log[%d]=%lf\n", i, logE[i]);
        }
    }
    double logDet = logMax - logMin;
    printf("logMax=%lf,logMin=%f,logDet=%lf\n", logMax, logMin, logDet);
    
    
    double *normaliseDouble = new double[SIZE];
    char *normaliseChar = new char[SIZE];
    for (int i=0; i<SIZE; i++) {
        normaliseDouble[i] = (logE[i]-logMin)/logDet;
        normaliseChar[i] = (char)(255.0*normaliseDouble[i]);
        if (i<20) {
            printf("normaliseDouble[%d]=%lf\n", i, normaliseDouble[i]);
            printf("normaliseChar[%d]=%d\n", i, normaliseChar[i]);
        }
    }
    
    
    cvNamedWindow("normaliseChar", 1);
    IplImage *ImageNMLS = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    ImageNMLS->imageData = normaliseChar;
    cvShowImage("normaliseChar", ImageNMLS);
    
    Complex *highpassComplex = new Complex[SIZE];
    double *highpass = new double[SIZE];
    char *highpasschar = new char[SIZE];
    for (int j=0; j<HEIGHT; j++) {
        for (int i=0; i<WIDTH; i++) {
            int idx = j*HEIGHT+i;
//            bool leftTop = (i<WIDTH/10) && (j<HEIGHT/10);
//            bool rightTop = (i>WIDTH*9/10) && (j<HEIGHT/10);
//            bool leftBottom = (i<WIDTH/10) && (j>HEIGHT*9/10);
//            bool RightBottom = (i>WIDTH*9/10) && (j>HEIGHT*9/10);
//            bool need = leftTop || rightTop || leftBottom || RightBottom;
            float dx = 0.0;
            float dy = 0.0;
            bool isLeft = i<=WIDTH/2;
            bool isUp = j<=HEIGHT/2;
            if (isLeft && isUp) {
                dx = i;
                dy = j;
            } else if (!isLeft && isUp) {
                dx = WIDTH-i;
                dy = j;
            } else if (!isLeft && !isUp) {
                dx = WIDTH-i;
                dy = HEIGHT-j;
            } else if (isLeft && !isUp) {
                dx = i;
                dy = HEIGHT-j;
            }
            float dist = sqrt(dx*dx+dy*dy);
            float need = dist<WIDTH/10;
            
            need = false;
            
            if (!need) {
                highpassComplex[idx].real = 1.0;//normaliseDouble[idx];
                highpassComplex[idx].imagin = 1.0;//normaliseDouble[idx];
                highpass[idx] = normaliseDouble[idx];
            } else {
                highpassComplex[idx].real = 0.0;
                highpassComplex[idx].imagin = 0.0;
                highpass[idx] = 0.0;
            }
            highpasschar[idx] = (char)(255.0*highpass[idx]);
            
        }
    }
    
    cvNamedWindow("highpass", 1);
    IplImage *Imagehighpass = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    Imagehighpass->imageData = highpasschar;
    cvShowImage("highpass", Imagehighpass);
    
    Complex *result = new Complex[SIZE];
    for (int i=0; i<SIZE; i++) {
        result[i].real = highpassComplex[i].real * dst[i].real;
        result[i].imagin = highpassComplex[i].imagin * dst[i].imagin;
    }
    
    Complex *result2= new Complex[SIZE];
    IFFT2D(result, result2, WIDTH, HEIGHT);
    
    char *imgIFFT2 = new char[SIZE];
    for (int i=0; i<SIZE; i++) {
        if (i<40) {
            printf("result2.real=%lf\n", result2[i].real);
        }
        imgIFFT2[i] = result2[i].real;
    }
    cvNamedWindow("ifft2", 1);
    IplImage *ImageIFFT2 = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    ImageIFFT2->imageData = imgIFFT2;
    cvShowImage("ifft2", ImageIFFT2);
    
    
    
    
    Complex *dst2= new Complex[SIZE];
    IFFT2D(dst, dst2, WIDTH, HEIGHT);
    
    char *imgIFFT = new char[SIZE];
    for (int i=0; i<SIZE; i++) {
        imgIFFT[i] = dst2[i].real;
    }
    cvNamedWindow("ifft", 1);
    IplImage *ImageIFFT = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    ImageIFFT->imageData = imgIFFT;
    cvShowImage("ifft", ImageIFFT);
    
//    for (int i=0; i<SIZE; i++) {
//        double dlt1 = (input[i] - dst2[i].real);
//        double dlt2 = (ImageOri->imageData[i] - dst2[i].real);
//        if (dlt1 > 0.0001) {
//            printf("aaaa double i %d dlt1 %lf\n", i, dlt1);
//        } else if (dlt2 > 0.0001) {
//            printf("aaaa char i %d, dlt2 %lf\n", i, dlt2);
//        } else {
//            printf("aaaa ok\n");
//        }
//    }
    

    
    
    cv::Mat image = cv::imread("/Users/admin/Desktop/beautyTest/dilireba.png", CV_LOAD_IMAGE_GRAYSCALE);
    imshow("oriCV", image);
    
    for (int i=0; i<SIZE; i++) {
        if (i<20) {
            printf("image[%d]=%u\n", i, image.at<unsigned char>(i));
        }
    }
    
    for (int i=0; i<256; i++) {
        image.at<unsigned char>(i) = i;
    }
    
    cv::Mat padded;
    int m = cv::getOptimalDFTSize(image.rows);  // Return size of 2^x that suite for FFT
    int n = cv::getOptimalDFTSize(image.cols);
    // Padding 0, result is @padded
    copyMakeBorder(image, padded, 0, m-image.rows, 0, n-image.cols, 0, cv::Scalar::all(0));
    cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F) };

    
    
    
    
    
    cv::Mat complexI;
    merge(planes, 2, complexI);
    
    dft(complexI, complexI);
    
//    cv::Mat inverseTransform;
//    complexI.copyTo(inverseTransform);
//    cv::Mat iPartDft[] = {cv::Mat::zeros(padded.size(),CV_32F),cv::Mat::zeros(padded.size(),CV_32F)};
//    cv::dft(complexI, inverseTransform, cv::DFT_INVERSE);
    
    split(complexI, planes);
    for (int j=0; j<HEIGHT; j++) {
        for (int i=0; i<WIDTH; i++) {
//            bool leftTop = (i<WIDTH/10) && (j<HEIGHT/10);
//            bool rightTop = (i>WIDTH*9/10) && (j<HEIGHT/10);
//            bool leftBottom = (i<WIDTH/10) && (j>HEIGHT*9/10);
//            bool RightBottom = (i>WIDTH*9/10) && (j>HEIGHT*9/10);
//            bool need = leftTop || rightTop || leftBottom || RightBottom;
            
            int idx = j*HEIGHT+i;
            if (idx < 256) {
//                printf("i=%d,j=%d, real=%f,imag=%f\n", i, j, planes[0].at<float>(i,j), planes[1].at<float>(i,j));
//                printf("i=%d,j=%d, real=%f,imag=%f\n", i, j, planes[0].at<float>(j,i), planes[1].at<float>(j,i));
                printf("idx=%d,real=%f,imag=%f\n", idx, planes[0].at<float>(idx), planes[1].at<float>(idx));
            }
            
            
            float dx = 0.0;
            float dy = 0.0;
            bool isLeft = i<=WIDTH/2;
            bool isUp = j<=HEIGHT/2;
            if (isLeft && isUp) {
                dx = i;
                dy = j;
            } else if (!isLeft && isUp) {
                dx = WIDTH-i;
                dy = j;
            } else if (!isLeft && !isUp) {
                dx = WIDTH-i;
                dy = HEIGHT-j;
            } else if (isLeft && !isUp) {
                dx = i;
                dy = HEIGHT-j;
            }
            float dist = sqrt(dx*dx+dy*dy);
            float need = dist<WIDTH/10;
            
            if (!need) {
                
            } else {
                planes[0].at<float>(i,j) = 0.0;
                planes[1].at<float>(i,j) = 0.0;
            }
        }
    }
    merge(planes, 2, complexI);
    
    
    cv::Mat ifftcv;
    idft(complexI, ifftcv, cv::DFT_REAL_OUTPUT);
    normalize(ifftcv,ifftcv,0,1,CV_MINMAX);
    imshow("idftcv",ifftcv);
    
//    split(inverseTransform, iPartDft);
//    imshow("ifftCVMat0", iPartDft[0]);
//    imshow("ifftCVMat1", iPartDft[1]);
    // Back to 8-bits
//    cv::Mat finalImage;
//    inverseTransform.convertTo(finalImage, CV_8U);
//    cvNamedWindow("ifftCV", 1);
//    IplImage *ImageIFFTCV = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
//    ImageIFFTCV->imageData = (char *)finalImage.data;
//    cvShowImage("ifftCV", ImageIFFTCV);
    
//    for (int i=0; i<SIZE; i++) {
//        if (i<10) {
//            printf("cvdft i=%d, complexI=%lf\n", i, complexI.at<float>(i));
//        }
//    }
    
    split(complexI, planes);
//    printf("SIZE=%d\n", SIZE);
//    int all = 0;
//    for (int j=0; j<HEIGHT; j++) {
//        for (int i=0; i<WIDTH; i++) {
//            if (j<1&&i<20) {
//                printf("cvdft i=%d, planes0=%lf, planes1=%lf\n",
//                       i, planes[0].at<float>(i,j), planes[1].at<float>(i,j));
//            }
//        }
//    }
//    printf("all=%d\n", all);
//    imshow("lalala0", planes[0]);
//    imshow("lalala1", planes[1]);
    
    magnitude(planes[0], planes[1], planes[0]);
//    for (int i=0; i<20; i++) {
//        printf("mag planes[0].at(%d)=%f\n", i, planes[0].at<float>(i));
//    }
    
    cv::Mat magI = planes[0];

    // => log(1+sqrt(Re(DFT(I))^2+Im(DFT(I))^2))
    magI += cv::Scalar::all(1);
//    for (int i=0; i<20; i++) {
//        printf("Scalar magI.at(%d)=%f\n", i, magI.at<float>(i));
//    }
    
    
    cv::log(magI, magI);
//    for (int i=0; i<20; i++) {
//        printf("log magI.at(%d)=%f\n", i, magI.at<float>(i));
//    }

    // crop the spectrum
    magI = magI(cv::Rect(0, 0, magI.cols & (-2), magI.rows & (-2)));
    int a = magI.cols;
    int b = magI.cols & (-2);
    int c = magI.rows;
    int d = magI.rows & (-2);
    printf("rect magI.cols=%d, magI.cols & (-2)=%d, magI.rows=%d, magI.rows & (-2)=%d\n", a, b, c, d);
    
    
    
    
    cv::Mat _magI = magI.clone();
    normalize(_magI, _magI, 0, 1, CV_MINMAX);
//    for (int i=0; i<20; i++) {
//        printf("normalize _magI.at(%d)=%f\n", i, _magI.at<float>(i));
//    }
    
    

    // rearrange the quadrants of Fourier image so that the origin is at the image center
    int cx = magI.cols/2;
    int cy = magI.rows/2;
    
    cv::Mat q0(magI, cv::Rect(0,0,cx,cy));    // Top-Left
    cv::Mat q1(magI, cv::Rect(cx,0,cx,cy));   // Top-Right
    cv::Mat q2(magI, cv::Rect(0,cy,cx,cy));   // Bottom-Left
    cv::Mat q3(magI, cv::Rect(cx,cy,cx,cy));  // Bottom-Right
    
    // exchange Top-Left and Bottom-Right
    cv::Mat tmp;
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);
    
    // exchange Top-Right and Bottom-Left
    q1.copyTo(tmp);
    q2.copyTo(q1);
    tmp.copyTo(q2);
    
    normalize(magI, magI, 0, 1, CV_MINMAX);
    
    imshow("Input image", image);
    imshow("Spectrum magnitude before shift frequency", _magI);
    imshow("Spectrum magnitude after shift frequency", magI);
    
    
    
    
    
    cvNamedWindow("fftCV", 1);
    IplImage *ImageFFTCV = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    ImageFFTCV->imageData = (char *)planes[0].data;
    cvShowImage("fftCV", ImageFFTCV);
    
    char *sub = new char[SIZE];
    for (int i=0; i<SIZE; i++) {
//        sub[i] = ImageFFTCV->imageData[i]-Image->imageData[i];
//        sub[i] = Image->imageData[i]-ImageFFTCV->imageData[i];
//        sub[i] = ImageFFTCV->imageData[i];
//        sub[i] = planes[0].data[i]-Image->imageData[i];

    }
    cvNamedWindow("sub", 1);
    IplImage *ImageSUB = cvCreateImage(cvSize(WIDTH,HEIGHT), IPL_DEPTH_8U, 1);
    ImageSUB->imageData = sub;
    cvShowImage("sub", ImageSUB);
    
    
    cvWaitKey();
    cvDestroyWindow("pic");
    cvReleaseImage(&ImagePic);
    cvDestroyWindow("ori");
    cvReleaseImage(&ImageOri);
    cvDestroyWindow("fft");
    cvReleaseImage(&Image);
    cvDestroyWindow("ifft");
    cvReleaseImage(&ImageIFFT);

    
    
    //IDFT(dst, NULL, SIZE);
    getchar();
}
#elif 0
#include<opencv2/opencv.hpp>
// chap 6 DFT and IDFT
int main()
{
    IplImage* src=cvLoadImage("/Users/admin/Desktop/beautyTest/dilireba.png",CV_LOAD_IMAGE_GRAYSCALE);   // src 8UC1
    IplImage* temp=cvCreateImage(cvGetSize(src),8,1);    // 中间变量，用于显示图像
    CvMat* srcMat=cvCreateMat(src->height,src->width,CV_64FC1);
    cvConvert(src,srcMat);              // 类型转换：图像指针转矩阵
    cvDFT(srcMat,srcMat,CV_DXT_FORWARD);    // DFT
    cvConvert(srcMat,temp);                 // 将矩阵转换为图像，方便显示DFT结果
    cvNamedWindow("DFT");
    cvShowImage("DFT",temp);
    
    cvZero(temp);           // 清零，方便下次使用
    cvDFT(srcMat,srcMat,CV_DXT_INV_SCALE); // IDFT
    cvConvert(srcMat,temp);             // 将 IDFT 结果矩阵转换为图像指针
    cvNamedWindow("IDFT");
    cvShowImage("IDFT",temp);           // IDFT
    cvWaitKey(0);
}
#elif 0
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


using namespace std;
using namespace  cv;


int main()
{
    Mat img = imread("/Users/admin/Desktop/beautyTest/dilireba.png",CV_LOAD_IMAGE_GRAYSCALE); //读入图像 转化为八位灰度图
    imshow("source",img); //原图像显示
    int h = cv::getOptimalDFTSize(img.rows); //图像尺寸的修正目的在于计算过程中可以优化，加快计算
    int w = cv::getOptimalDFTSize(img.cols);
    Mat padded;
    cv::copyMakeBorder(img,padded,0,h-img.rows,0,w-img.cols,BORDER_CONSTANT,::Scalar::all(0));
    //扩充成最优尺寸,注意：此处是复制了img的数据到padded，并非简单通过指针和img共享内存
    //但这种情况则没有复制: 因为dst本身就和img共享了内存，只是简单的扩充了边界
    //Mat img = imread("xxx.jpg");
    //Mat dst = img;
    //cv::copyMakeBorder(img,dst,2,2,2,2,BORDER_CONSTANT,Scalar::all(0));
    
    
    Mat plane[] = {Mat_<float>(padded),Mat::zeros(padded.size(),CV_32F)}; //实部为img，虚部填充0
    Mat complexImg;
    cv::merge(plane,2,complexImg); //可以理解为组合成2通道（实部+虚部）图像
    dft(complexImg,complexImg); //DFT变换后的数据复制到原处，没有另外开辟内存
    split(complexImg,plane); //实部plane[0] 虚部plane[1]
    imshow("lalala1",plane[0]);
    
    
    //--------------------------------------------------------------------------
    //dft变换后利用idft变换恢复原图像，显示效果与原图像相差无几，
    Mat frequencyImg;
    complexImg.copyTo(frequencyImg);
    Mat iDft[] = {Mat::zeros(padded.size(),CV_32F),Mat::zeros(padded.size(),CV_32F)};
    ::idft(frequencyImg,frequencyImg);
    ::split(frequencyImg,iDft);
    ::magnitude(iDft[0],iDft[1],iDft[0]);
    cv::normalize(iDft[0],iDft[0],1,0,CV_MINMAX);
    imshow("iDFT",iDft[0]); //dft idft 后的图像
    //-------------------------------------------------------------------------
    
    
    //-------------------------------------------------------------------------------
    //根据DFT的能力集中地特点，即图像信息集中在低频部分
    //先取出图像四个角（低频部分）的小块区域（其余灰度置为0），在idft反变换回去，和原图像进行比较
    //显示结果表明只取出反变换后的16%信息，在利用idft反变换，发现与原图相差不太大
    //这是图像处理中利用正交变换（包括傅里叶变换）的能量集中特性（图像有损压缩中有利用）
    Mat partFrequencyImg;
    complexImg.copyTo(partFrequencyImg);
    int nx1 = int(0.2f*padded.cols);
    int nx2 = int(0.8f*padded.cols);
    int ny1 = int(0.2f*padded.rows);
    int ny2 = int(0.8f*padded.rows);
    partFrequencyImg.colRange(nx1,nx2).setTo(Scalar::all(0));
    partFrequencyImg.rowRange(ny1,ny2).setTo(Scalar::all(0));
    Mat iPartDft[] = {Mat::zeros(padded.size(),CV_32F),Mat::zeros(padded.size(),CV_32F)};
    ::idft(partFrequencyImg,partFrequencyImg);
    ::split(partFrequencyImg,iPartDft);
    ::magnitude(iPartDft[0],iPartDft[1],iPartDft[0]);
    cv::normalize(iPartDft[0],iPartDft[0],1,0,CV_MINMAX);
    imshow("iPartDFT",iPartDft[0]); //取dft变换后的部分低频部分，利用idft反变换后的图像
    //--------------------------------------------------------------------------
    
    
    ::split(complexImg,plane); //实部plane[0] 虚部plane[1]
    imshow("lalala",plane[0]);
    
    
    ::magnitude(plane[0],plane[1],plane[0]);//计算幅频 结果保存到plane[0]
    imshow("dft",plane[0]); //显示出来的幅频图，看不出任何有效信息，需要做增强处理
    
    
    plane[0] += Scalar::all(1); //整体加1的目的在于log计算后非负
    cv::log(plane[0],plane[0]); //log可起到拉伸，增大对比度的作用
    cv::normalize(plane[0],plane[0],1,0,CV_MINMAX);//整体映射到0~1,也称归一化处理
    imshow("enhance",plane[0]); //增强后的图像，显示效果明显改善
    
    
    //保证plane[0]的行列均为偶数 以下将低频部分往中心移动
    int cx = plane[0].cols / 2;
    int cy = plane[0].rows / 2;
    Mat m1(plane[0],cv::Rect(0,0,cx,cy)); //左上部分
    Mat m2(plane[0],cv::Rect(cx,0,cx,cy)); //右上部分
    Mat m3(plane[0],cv::Rect(0,cy,cx,cy)); //左下部分
    Mat m4(plane[0],cv::Rect(cx,cy,cx,cy)); //右下部分
    Mat temp;
    m1.copyTo(temp); //左上与右下交换
    m4.copyTo(m1);
    temp.copyTo(m4);
    m2.copyTo(temp);// 右上与左下交换
    m3.copyTo(m2);
    temp.copyTo(m3);
    imshow("shift",plane[0]); //低频中心化的图像
    waitKey(0);
}
#endif
