#define	OUTIMG
// 4 byte reversible DCT from 2006
// gcc -g -O -o idct_macros idct_macros.c -ltiff -lpng -ljpeg -lm
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

//want to fiddle with token concatenation so that I can define the names of the variables below
//otherwise I could store these in a 2D matrix....
#define BUTTERFLY(x0, x1, x2, x3,var0,var1,var2,var3) do {\
	     var0 = x0 + x3;\
	     var1 = x1 + x2;\
	     var2 = x1 - x2;\
	     var3 = x0 - x3;\
           } while(0)

#define TRANSFORM(i0,i1,i2,i3,var0,var1,var2,var3) do {\
		var0 = i0 + i1;\
		var2 = i0 - i1;\
		var1 = i2 + (i3<<1);\
		var3 = i3 - (i2<<1);\
	    } while(0)

#define REVERSETRANSFORM(t0,t1,t2,t3,var0,var1,var2,var3) do {\
		var0 = t0 + t2;\
		var1 = t0 - t2;\
		var0 += (I0 >> 2);\
		var1 += (I1 >> 2);\
		var2 = (t1>>1) - t3;\
		var3 = t1 + (t3>>1);\
		} while(0)

#define RECOVER(I0,I1,I2,I3,var0,var1,var2,var3) do {\
		r0 = (I0 + I3)/10;\
		r1 = (I1 + I2)/10;\
		r2 = (I1 - I2)/10;\
		r3 = (I0 - I3)/10;\
		} while(0)


int hist[HS];	// original image histogram
int hist0[HS];	// DCT 0...
int hist1[HS];
int hist2[HS];
int hist3[HS];
double imgent, dctent;

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

main(int argc, char *argv[]) {
	int x0, x1, x2, x3; // 4 original pixels
	int i0, i1, i2, i3; // intermediate temp vars
	int t0, t1, t2, t3; // 4 DCT terms
	int r0, r1, r2, r3; // 4 recovered values
	int I0, I1, I2, I3; // recovery intermediate vars
	int i, n;
	unsigned char *ip;
	fprintf(stderr, "file %s\n", argv[1]);
	im0 = read_img(argv[1]);
	fprintf(stderr, "%d * %d  %d bpp\n", im0->wid, im0->ht, im0->bpp);
#ifdef OUTIMG
	printf("P5\n%d %d\n255\n", im0->wid, im0->ht);
#endif
	ip = im0->pp;
	n = im0->wid * im0->ht;
	for(i = 0; i < n; i += 4) { // width must be mpl of 4
		x0 = *ip++;
		x1 = *ip++;
		x2 = *ip++;
		x3 = *ip++;
hist[x0]++;
hist[x1]++;
hist[x2]++;
hist[x3]++;
		BUTTERFLY(x0,x1,x2,x3,i0,i1,i2,i3);		
		//i0 = x0 + x3;
		//i1 = x1 + x2;
		//i2 = x1 - x2;
		//i3 = x0 - x3;
		
		TRANSFORM(i0,i1,i2,i3,t0,t1,t2,t3);
		//t0 = i0 + i1;
		//t2 = i0 - i1;
		//t1 = i2 + (i3<<1);
		//t3 = i3 - (i2<<1);
hist0[t0]++;
hist1[t1+HS/2]++;
hist2[t2+HS/2]++;
hist3[t3+HS/2]++;
		// mpy terms by 2x to simplify later arith
		t0 += t0;
		t1 += t1;
		t2 += t2;
		t3 += t3;
		
		REVERSETRANSFORM(t0,t1,t2,t3,I0,I1,I2,I3);
		//I0 = t0 + t2;
		//I1 = t0 - t2;
		//I0 += (I0 >> 2);
		//I1 += (I1 >> 2);
		//I2 = (t1>>1) - t3;
		//I3 = t1 + (t3>>1);
		
		// recovered values r0-r3
		RECOVER(I0,I1,I2,I3,r0,r1,r2,r3);		
		//r0 = (I0 + I3)/10; // test mpy and shift alternatives to div
		//r1 = (I1 + I2)/10;
		//r2 = (I1 - I2)/10;
		//r3 = (I0 - I3)/10;
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
	imgent += entropy("image", hist, HS);
	dctent += entropy("DCT0", hist0, HS);
	dctent += entropy("DCT1", hist1, HS);
	dctent += entropy("DCT2", hist2, HS);
	dctent += entropy("DCT3", hist3, HS);
	fprintf(stderr, "sum of dct ents %g\n", dctent);
	fprintf(stderr, "full img ent %g  full dct ents %g\n",
		imgent*n, dctent*n/4);
}
