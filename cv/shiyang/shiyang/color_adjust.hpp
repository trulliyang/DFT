//
//  color_adjust.hpp
//  test
//
//  Created by 郭秀艳 on 2017/6/21.
//  Copyright © 2017年 郭秀艳. All rights reserved.
//

#ifndef color_adjust_hpp
#define color_adjust_hpp
#define BRIGHTNESS          120
#define MIN_GAIN            4
#define MAX_GAIN            15
#define MIN2(a, b)          ((a > b) ? b : a)
#define MAX2(a, b)          ((a > b) ? a : b)
#define CLIP3(x, min, max)  ((x > min) ? (x < max ? x : max) : min)
#include <stdio.h>
void getCurve(unsigned char *target_curve, int* average_brightness);
int color_adjust_RGB(unsigned char *pSrc, unsigned char *pOut, int nWidth, int nHeight, int Strie, int* average_brightness);
#endif /* color_adjust_hpp */
