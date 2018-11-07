/*
  llzlab - luolongzhi algorithm lab 
  Copyright (C) 2013 luolongzhi 罗龙智 (Chengdu, China)

  This program is part of llzlab, all copyrights are reserved by luolongzhi. 

  filename: main.c 
  time    : 2012/07/14 22:14 
  author  : luolongzhi ( luolongzhi@gmail.com )
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "llz_mdct.h"
#include "llz_mdct_fixed.h"

#define N  256 //32
#define N2 (N >> 1)

int main(int argc, char *argv[])
{
	int i;
    int type = 0;

	double buf_in[N];
	double buf_out[N];
    double mdct_buf[N2];

	int   buf_in_int[N];
	int   buf_out_int[N];
    int   mdct_buf_int[N2];

	uintptr_t handle;
	uintptr_t handle_int;

	for (i = 0 ; i < N ; i++) {
		/*buf_in[i] = (1<<25)+i; //i*i;*/
		buf_in[i] = (1<<21)+i+1; //i*i;
		buf_in_int[i] = (1<<21)+i+1; //i*i;
    }

	printf("\n");

	handle     = llz_mdct_init(type, N);
	handle_int = llz_mdct_fixed_init(type, N);

    memset(mdct_buf, 0, sizeof(double)*N2);
    memset(mdct_buf_int, 0, sizeof(int)*N2);

	llz_mdct(handle, buf_in, mdct_buf);
	llz_mdct_fixed(handle_int, buf_in_int, mdct_buf_int);
    for(i = 0; i < N2; i++) 
        printf("mdct[%d] = %f\t mdct_int[%d] = %d\t err_mdct = %f%%\n", 
                i, mdct_buf[i], i, mdct_buf_int[i], 
                100*fabs(((double)(mdct_buf[i])-(double)(mdct_buf_int[i]))/mdct_buf[i]));

    memset(buf_out, 0, sizeof(double)*N);
    memset(buf_out_int, 0, sizeof(int)*N);
	llz_imdct(handle, mdct_buf, buf_out);
	llz_imdct_fixed(handle_int, mdct_buf_int, buf_out_int);

    for (i = 0 ; i < N ; i++)
        printf("double(%f\t%f)\t int(%d\t%d)\t err_imdct = %f%%\n",
                buf_in[i], buf_out[i], buf_in_int[i], buf_out_int[i], 
                100*fabs(((double)(buf_out[i])-(double)(buf_out_int[i]))/buf_out[i]));


	llz_mdct_uninit(handle);
	llz_mdct_fixed_uninit(handle_int);

    return 0;
}

