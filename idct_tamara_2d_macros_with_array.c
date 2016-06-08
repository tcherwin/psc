#define	OUTIMG 
// 4 byte reversible DCT from 2006
// gcc -O3 -o idct_tamara_2d idct_tamara_2d.c -ltiff -lpng -ljpeg -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <jpeglib.h>
#include <tiffio.h>
#include "/home/awetzel/bin/myimio.h"
struct image *im0, *im1;
#define	HS	4096 // big enough for largest histogram range



//can optimize math within here (ie, subtract 4 and add "height" multiplier at each step if necessary
#define NEXTBLOCK(height_add,ip, im0, pp) do {\
		ip = im0->pp;\
		ip +=i;\
		ip +=im0->ht*(j + height_add);\
		x0 = *ip++;\
		x1 = *ip++;\
		x2 = *ip++;\
		x3 = *ip++;\
	     } while(0)



#define BUTTERFLY(index) do {\
	     i0 = *block_array_x[index,0] + *block_array_x[index,3];\
	     i1 = *block_array_x[index,1] + *block_array_x[index,2];\
	     i2 = *block_array_x[index,1] - *block_array_x[index,2];\
	     i3 = *block_array_x[index,0] - *block_array_x[index,3];\
	     printf("%d",index);\
	} while(0)

#define TRANSFORM() do {\
		t0 = i0 + i1;\
		t2 = i0 - i1;\
		t1 = i2 + (i3<<1);\
		t3 = i3 - (i2<<1);\
	    } while(0)

#define REVERSETRANSFORM(index) do {\
		I0 = *block_array_t[index,0] + *block_array_t[index,2];\
		I1 = *block_array_t[index,0] - *block_array_t[index,2];\
		I0 += (I0 >> 2);\
		I1 += (I1 >> 2);\
		I2 = (*block_array_t[index,1]>>1) - *block_array_t[index,3];\
		I3 = *block_array_t[index,1] + (*block_array_t[index,3]>>1);\
		} while(0)
		
#define RECOVER(index) do {\
		r0 = (I0 + I3)/10;\
		r1 = (I1 + I2)/10;\
		r2 = (I1 - I2)/10;\
		r3 = (I0 - I3)/10;\
		} while(0)

//if I don't want to create this multitude of variables and instead store values in an array
/*
#define BLOCK_PLACE(currentrow, currentarray,var0,var1,var2,var3) do {\
		currentarray[currentrow,var0]\
		currentarray[currentrow
	    } while(0)

*/


int hist[HS];	// original image histogram
int hist0[HS];	// DCT 0...
int hist1[HS];
int hist2[HS];
int hist3[HS];
double imgent, dctent;

int block_array_x[4][4];
int block_array_i[4][4];
int block_array_t[4][4];
int block_array_t2[4][4];
int block_array_r[4][4];
int block_array_I[4][4];



double entropy(char *title, int *hp, int n) { // hp -> integerized historgram
	long i, count = 0;
	double p, e = 0.0;
	for(i = 0; i < n; i++)
		count += hp[i];
	for(i = 0; i < n; i++) {
		if(hp[i]) {
			p = (float)hp[i]/count;
			e += p*log(p);
		}
	}
	e /= log(2.0);
	if(title)
		fprintf(stderr, "%s\t%ld\t%f\n", title, count, e);
	return(e);
}


//if we choose to expand block size, should use matrices instead of all these named variables

main(int argc, char *argv[]) {
	int x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15; // 16 original pixels
	int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15; // 16 intermediate temp vars
	int t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15; // 16 DCT terms
	
	//int x2_0, x2_1, x2_2, x2_3, x2_4, x2_5, x2_6, x2_7, x2_8, x2_9, x2_10, x2_11, x2_12, x2_13, x2_14, x2_15; // (1D DCT terms for transform)
	int t2_0, t2_1, t2_2, t2_3, t2_4, t2_5, t2_6, t2_7, t2_8, t2_9, t2_10, t2_11, t2_12, t2_13, t2_14, t2_15; // 16 DCT terms (2D)

	int r0, r1, r2, r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15; // 4 recovered values
	int I0, I1, I2, I3, I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15; // recovery intermediate vars
	int i, n,j;
	
	

	unsigned char *ip;
	fprintf(stderr, "file %s\n", argv[1]);
	im0 = read_img(argv[1]);
	fprintf(stderr, "%d * %d  %d bpp\n", im0->wid, im0->ht, im0->bpp);
#ifdef OUTIMG
	printf("P5\n%d %d\n255\n", im0->wid, im0->ht);
#endif
	//ip = im0->pp; //unnecessary because of statement within for loop
	n = im0->wid * im0->ht;
	int x = im0->ht;
	int y = im0->wid;
	printf("determining divisibility by 4");
	int xmod = x%4; //will need to pad with 4 - xmod	
	int ymod = y%4; //will need to pad with 4 - ymod
	int numblocksx = ceil(im0->wid/4); //after padding, number of blocks will equal ceiling
	int numblocksy = ceil(im0-> ht/4);	
		
	printf("x mod 4 is\t %d", xmod);
	printf("y mod 4 is \t %d", ymod);
	printf("x blocks is  %d", numblocksx);
	printf("y blocks is  %d", numblocksy);
	//pad here

	//need to consider uneven images
	


	for(i = 0; i < x; i += 4) 
	{ // traverse width		
		//printf("new column");
		for (j = 0; j <y; j +=4)
		{
			
			// **** 1D *****			

			//get pixels
			NEXTBLOCK(0,ip,im0,pp);	
			block_array_x[0,x0];
			block_array_x[0,x1];
			block_array_x[0,x2];
			block_array_x[0,x3];
			
			//printf("\nrow 1 %d %d %d %d", x0, x1, x2, x3);
			NEXTBLOCK(1,ip,im0,pp);
			
			block_array_x[1,x0];
			block_array_x[1,x1];
			block_array_x[1,x2];
			block_array_x[1,x3];

			//printf("\nrow 2 %d %d %d %d", x4, x5, x6, x7);		
			NEXTBLOCK(2,ip,im0,pp);
			
			block_array_x[2,x0];
			block_array_x[2,x1];
			block_array_x[2,x2];
			block_array_x[2,x3];


			//printf("\nrow 3 %d %d %d %d", x8, x9, x10, x11);		
			NEXTBLOCK(3,ip,im0,pp);
			
			block_array_x[3,x0];
			block_array_x[3,x1];
			block_array_x[3,x2];
			block_array_x[3,x3];			

			//printf("\nrow 4 %d %d %d %d", x12, x13, x14, x15);
			
			//butterfly
			BUTTERFLY(0);
			
			//block_array_i[0,i0];
			//block_array_i[0,i1];
			//block_array_i[0,i2];
			//block_array_i[0,i3];

			BUTTERFLY(1);

			//block_array_i[1,i0];
			//block_array_i[1,i1];
			//block_array_i[1,i2];
			//block_array_i[1,i3];

			
			BUTTERFLY(2);
			
			//block_array_i[2,i0];
			//block_array_i[2,i1];
			//block_array_i[2,i2];
			//block_array_i[2,i3];

			BUTTERFLY(3);	

			//block_array_i[3,i0];
			//block_array_i[3,i1];
			//block_array_i[3,i2];
			//block_array_i[3,i3];



			
			//1D transform
			TRANSFORM(0);
			
			block_array_t[0,t0];
			block_array_t[0,t1];
			block_array_t[0,t2];
			block_array_t[0,t3];


			TRANSFORM(1);
			
			block_array_t[1,t0];
			block_array_t[1,t1];
			block_array_t[1,t2];
			block_array_t[1,t3];

			TRANSFORM(2);
			
			block_array_t[2,t0];
			block_array_t[2,t1];
			block_array_t[2,t2];
			block_array_t[2,t3];


			TRANSFORM(3);
			
			block_array_t[3,t0];
			block_array_t[3,t1];
			block_array_t[3,t2];
			block_array_t[3,t3];


			//reverse 1D

			// mpy terms by 2x to simplify later arith
		
			t0 += t0;
			t1 += t1;
			t2 += t2;
			t3 += t3;
			t4 += t4;
			t5 += t5;
			t6 += t6;
			t7 += t7;	
			t8 += t8;
			t9 += t9;
			t10 += t10;
			t11 += t11;
			t12 += t12;
			t13 += t13;
			t14 += t14;
			t15 += t15;

			REVERSETRANSFORM(1);
			
			block_array_I[0,I0];
			block_array_I[0,I1];
			block_array_I[0,I2];
			block_array_I[0,I3];
			
			REVERSETRANSFORM(2);
			
			

			REVERSETRANSFORM(3);
			


			REVERSETRANSFORM(4);

			//need to think about how to "unrotate" the matrix for recovery	
			//recover values from 2D to 1D

			RECOVER(1);
			


			RECOVER(2);



			RECOVER(3);



			RECOVER(4);

			//if(x0 != r0 || x1 != r1 || x2 != r2 || x3 != r3) // does it really match?
//printf("%d %d %d %d     %d %d %d %d    %d %d %d %d\n",x0, x1, x2, x3, t0, t1, t2, t3, r0, r1, r2, r3);

//if(x4 != r4 || x5 != r5 || x6 != r6 || x7 != r7) // does it really match?
//printf("%d %d %d %d     %d %d %d %d    %d %d %d %d\n",x4, x5, x6, x7, t4, t5, t6, t7, r4, r5, r6, r7);
			//if(t0 != r0 || t1 != r1 || t2 != r2 || t3 != r3 || t4 != r4 || t5 != r5 || t6 !=r6 || t7 !=r7 || t8 !=r8 || t9 !=r9 || t10 !=r10 || t11 != r11 || t12 !=r12 || t13 != r13 || t14 !=r14 || t15 !=r15)
	//printf("FAIL");			

			// ***** 2D **** 
	

			//convention: rotate such that 0  1  2  3
			//			       4  5  6  7 	
			//                             8  9  10 11 
			//			       12 13 14 15
			                             


			// becomes                     12  8  4  0
			//			       13  9  5  1
			//                             14  10 6  2
			//			       15  11 7  3 	

			//reuse intermediate temp variables

			/*
			BUTTERFLY(t12,t8,t4,t0,i12,i8,i4,i0);
			BUTTERFLY(t13,t9,t5,t1,i13,i9,i5,i1);
			BUTTERFLY(t14,t10,t6,t2,i14,i10,i6,i2);
			BUTTERFLY(t15,t11,t7,t3,i15,i11,i7,i3);
	
			TRANSFORM(i12,i8,i4,i0,t2_12,t2_8,t2_4,t2_12);
			TRANSFORM(i13,i9,i5,i1,t2_13,t2_9,t2_5,t2_1);
			TRANSFORM(i14,i10,i6,i2,t2_14,t2_10,t2_6,t2_2);
			TRANSFORM(i15,i11,i7,i3,t2_15,t2_11,t2_7,t2_3);
		
			//*** end 1D transform ****
	

			//*** attempting to reverse, *****

			//reversing 2D
		

			// mpy terms by 2x to simplify later arith
		
			t2_0 += t2_0;
			t2_1 += t2_1;
			t2_2 += t2_2;
			t2_3 += t2_3;
			t2_4 += t2_4;
			t2_5 += t2_5;
			t2_6 += t2_6;
			t2_7 += t2_7;	
			t2_8 += t2_8;
			t2_9 += t2_9;
			t2_10 += t2_10;
			t2_11 += t2_11;
			t2_12 += t2_12;
			t2_13 += t2_13;
			t2_14 += t2_14;
			t2_15 += t2_15;
		
			

			//reverse values from 2D to 1D			

			REVERSETRANSFORM(t2_12,t2_8,t2_4,t2_0,I12,I8,I4,I0);
			REVERSETRANSFORM(t2_13,t2_9,t2_5,t2_1,I13,I9,I5,I1);
			REVERSETRANSFORM(t2_14,t2_10,t2_6,t2_2,I14,I10,I6,I2);
			REVERSETRANSFORM(t2_15,t2_11,t2_7,t2_3,I15,I11,I7,I3);

			//need to think about how to "unrotate" the matrix for recovery	
			//recover values from 2D to 1D

			RECOVER(I12,I8,I4,I0,r12,r8,r4,r0);
			RECOVER(I13,I9,I5,I1,r13,r9,r5,r1);
			RECOVER(I14,I10,I6,I2,r14,r10,r6,r2);
			RECOVER(I15,I11,I7,I3,r15,r11,r7,r3);

			//printf("checking recovered 2D");

			//need to think about how to "unrotate" the matrix for recovery				
			

			if(t0!= r0)
			{
				//printf("FAIL");
				printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d     %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d      %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d  \n",t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t2_0,t2_1,t2_2,t2_3,t2_4,t2_5,t2_6,t2_7,t2_8,t2_9,t2_10,t2_11,t2_12,t2_13,t2_14,t2_15,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
			};			
			//if(t0 != r12 || t1 != r8 || t2 != r4 || t3 != r0 || t4 != r13 || t5 != r9 || t6 !=r5 || t7 !=r1 || t8 !=r14 || t9 !=r10 || t10 !=r6 || t11 != r2 || t12 !=r15 || t13 != r11 || t14 !=r7 || t15 !=r3)
			*/			

			//printf("FAIL");
			/*if(t0 != r0 || t1 != r1 || t2 != t2 || t3 != t3) // does it really match?
			printf("%d %d %d %d     %d %d %d %d    %d %d %d %d\n",
			t0, t1, t2, t3, t2_0, t2_1, t2_2, t2_3, r0, r1, r2, r3);
			*/


		}//end with for
		
	} //end height for (take this paranthesis out later)











		/*

		t0 += t0;
		t1 += t1;
		t2 += t2;
		t3 += t3;





		I0 = t0 + t2;
		I1 = t0 - t2;
		I0 += (I0 >> 2);
		I1 += (I1 >> 2);
		I2 = (t1>>1) - t3;
		I3 = t1 + (t3>>1);
		// recovered values r0-r3
		r0 = (I0 + I3)/10; // test mpy and shift alternatives to div
		r1 = (I1 + I2)/10;
		r2 = (I1 - I2)/10;
		r3 = (I0 - I3)/10;
		

	
		/*
hist[x0]++;
hist[x1]++;
hist[x2]++;
hist[x3]++;

		
		

hist0[t0]++;
hist1[t1+HS/2]++;
hist2[t2+HS/2]++;
hist3[t3+HS/2]++;
		
			// mpy terms by 2x to simplify later arith
		//t0 = t1 = t2 =t3=0;
		//t0 = 100;
		t0 += t0;
		t1 += t1;
		t2 += t2;
		t3 += t3;
		I0 = t0 + t2;
		I1 = t0 - t2;
		I0 += (I0 >> 2);
		I1 += (I1 >> 2);
		I2 = (t1>>1) - t3;
		I3 = t1 + (t3>>1);
		// recovered values r0-r3
		r0 = (I0 + I3)/10; // test mpy and shift alternatives to div
		r1 = (I1 + I2)/10;
		r2 = (I1 - I2)/10;
		r3 = (I0 - I3)/10;

if(x0 != r0 || x1 != r1 || x2 != r2 || x3 != r3) // does it really match?
printf("%d %d %d %d     %d %d %d %d    %d %d %d %d\n",
x0, x1, x2, x3, t0, t1, t2, t3, r0, r1, r2, r3);


//printf("%d %d %d %d\n", i0, i1, i2, i3);
//printf("%d %d %d %d\n", I0, I1, I2, I3);
#ifdef	OUTIMG
		// output the reconstructed image
		putchar(r0);
		putchar(r1);
		putchar(r2);
		putchar(r3);
#endif
	}
#ifndef OUTIMG
	for(i = 0; i < HS; i++)
		printf("%d: %d %d %d %d %d\n", i, hist[i],
			hist0[i], hist1[i], hist2[i], hist3[i]);
#endif
	imgent = entropy("image", hist, HS);
	dctent = entropy("DCT0", hist0, HS);
	dctent += entropy("DCT1", hist1, HS);
	dctent += entropy("DCT2", hist2, HS);
	dctent += entropy("DCT3", hist3, HS);
	fprintf(stderr, "sum of dct ents %g\n", dctent);
	fprintf(stderr, "full img ent %g  full dct ents %g\n",
		imgent*n, dctent*n/4);

*/
}
