//
//  color_adjust.cpp
//  test
//
//  Created by 郭秀艳 on 2017/6/21.
//  Copyright © 2017年 郭秀艳. All rights reserved.
//

#include "color_adjust.hpp"
#include <math.h>
#include <memory.h>
void get_single_curve(unsigned char* target_curve,
                      int x1, int x2,
                      int gain1, int gain2);

int color_adjust_RGB(unsigned char *pSrc, unsigned char *pOut,
                     int nWidth, int nHeight,
                     int Strie, int* average_brightness)
{
    int nRet = 0;
    int i, j;
    //int currrent_brightness = BRIGHTNESS;
    unsigned char target_curve[3][256];
    memset(target_curve, 0, 256);
    
    getCurve(&target_curve[0][0], average_brightness);
    for(i=0; i<3; i++) {
        for(j=0; j<256; j++) {
            printf("%d char is :%d \n", j, target_curve[i][j]);
        }
    }
    unsigned char *pLineHeader = pSrc;
    unsigned char *pR = pLineHeader + 2;
    unsigned char *pG = pLineHeader + 1;
    unsigned char *pB = pLineHeader;
    for(i=0; i<nHeight; i++)
    {
        for(j=0; j<nWidth; j++)
        {
            if (j==12) {
                printf("before *pR=%d *pG=%d *pB=%d", *pR, *pG, *pB);
            }
            
            *pR = 255;//target_curve[0][*pR];
            *pG = 0;//target_curve[1][*pG];
            *pB = 0;//target_curve[2][*pB];
            if (j==12) {
                printf(" after *pR=%d *pG=%d *pB=%d\n", *pR, *pG, *pB);
            }
            pR += 3;
            pG += 3;
            pB += 3;
        }
    }
    return nRet;
}

void getCurve(unsigned char *target_curve, int *average_brightness)
{
    //int currrent_brightness = BRIGHTNESS;

    int x1 = average_brightness[0] * 0.3;
    int x2 = average_brightness[0] + (255 - average_brightness[0]) * 0.5;
    
    //////first///////////////
    int gainr1 = x1/8;
    gainr1 = MAX2(gainr1, MIN_GAIN);
    
    ///////////second/////////////
    int gainr2 = gainr1 * (1.1 + (255 - average_brightness[0]) / 255.0);
    gainr2 = MIN2(gainr2, MAX_GAIN);
    
    get_single_curve(target_curve, x1, x2, gainr1, gainr2);
    
    x1 = average_brightness[1] * 0.3;
    x2 = average_brightness[1] + (255 - average_brightness[1]) * 0.5;
    
    int gaing2 = gainr1 * (1.1 + (255 - average_brightness[1]) / 255.0);
    gaing2 = MIN2(gaing2, gainr2 * 0.8);
    get_single_curve(target_curve + 256, x1, x2, gainr1, gaing2);
    
    x1 = average_brightness[2] * 0.3;
    x2 = average_brightness[2] + (255 - average_brightness[2]) * 0.5;
    int gainb2 = gainr1 * (1.1 + (255 - average_brightness[2]) / 255);
    gainb2 = MIN2(gainb2, gainr2);
    get_single_curve(target_curve + 512, x1, x2, gainr1, gainb2);
    return;

}

void get_single_curve(unsigned char* target_curve, int x1, int x2, int gain1, int gain2)
{
    float a, b, c;
    int i;
    //////first///////////////
    a = 4.0 * gain1 / (x1 * x1);
    b = 1 - a * x1;
    //c = 0;
    for(i=0; i<x1; i++)
    {
        target_curve[i] = a * i * i + b * i;// + c;
    }
    
    float dis = (x2 - x1) * 0.5;
    a = gain2 / (float)(dis * dis - gain2 * gain2);
    b = 1 - a * (x1 +x2);
    c = a * x1 * x2;
    target_curve[x1] = x1;
    int resualt;
    for(i=x1 + 1; i<x2; i++)
    {
        resualt = (sqrt(b * b - 4 * a * (c - i)) - b) / (2 * a) + 0.5;
        target_curve[i] = CLIP3(resualt, 0, 255);
    }
    
    for(i = x2; i<256; i++)
    {
        target_curve[i] = i;
    }
    
    
    return;


}
