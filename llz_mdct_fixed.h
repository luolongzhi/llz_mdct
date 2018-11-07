/*
  llzlab - luolongzhi algorithm lab 
  Copyright (C) 2013 luolongzhi 罗龙智 (Chengdu, China)

  This program is part of llzlab, all copyrights are reserved by luolongzhi. 

  filename: llz_mdct_fixed.h 
  time    : 2012/07/19 - 2012/07/20  
  author  : luolongzhi ( luolongzhi@gmail.com )
*/

#ifndef _LLZ_MDCT_FIXED_H
#define _LLZ_MDCT_FIXED_H


#include <stdint.h>
//typedef unsigned uintptr_t;

/*
    origin: the naive mdct using the original formular, this mdct is leading to you learning mdct
    fft   : normally called fast mdct, use fft transform to compute mdct
    fft4  : the wildly used, N/4 point FFT to compute mdct
*/
enum{
    MDCT_FIXED_ORIGIN = 0,
    MDCT_FIXED_FFT,
    MDCT_FIXED_FFT4,
};

uintptr_t llz_mdct_fixed_init(int type, int len);
void      llz_mdct_fixed_uninit(uintptr_t handle);

void llz_mdct_fixed(uintptr_t handle, int *x, int *X);
void llz_imdct_fixed(uintptr_t handle, int *X, int *x);

#endif
