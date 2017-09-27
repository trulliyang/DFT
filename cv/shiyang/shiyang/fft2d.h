typedef struct{
    float real;
    float image;
}complex;

void shiyangfft2d();
complex Add(complex, complex);
complex Sub(complex, complex);
complex Mul(complex, complex);
int calculate_M(int);
void reverse(int,int);
void readData();
void fft(int,int);
void Ifft();
void printResult_fft();
void printResult_Ifft();
