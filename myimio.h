#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <jpeglib.h>
#include <tiffio.h>
#include <png.h>

#define uchar unsigned char
#define ushort unsigned short

struct image {
	unsigned char *pp;
	int wid, ht, ydelta, bpp;
};

struct image *newimage(int wid, int ht, int bpp) {
	struct image *ip = (struct image *)malloc(sizeof(struct image));
	ip->wid = wid;
	ip->ht = ht;
	ip->bpp = bpp;
	ip->ydelta = wid;
	ip->pp = (unsigned char *)malloc((ht+1)*(long)wid*bpp); // extra yguard
//fprintf(stderr, "newimage 0x%lx  %d %d % d\n", ip->pp, wid, ht, bpp);
	if(!ip->pp) {
		fprintf(stderr, "image malloc failed %d %d %d\n", ht, wid, bpp);
		return(NULL);
	}
	return(ip);
}

int tiffout = 0;
int ltiff_reader(char *fname, int *widp, int *htp, int *bpp, uchar **bp) {
	TIFF *timg;
	uint16 bps, spp, phot;
	unsigned char *buff;
	int i, wid, ht, stsize, stmax, strips, bytesperpixel;
	unsigned long qread, got;
	unsigned long buffsz, count;
	char pgmhdr[100];

	TIFFSetWarningHandler(NULL);
	if((timg = TIFFOpen(fname, "r")) == NULL) {
		fprintf(stderr, "Could not open input: %s\n", fname);
		return(-1);
	}

	if(TIFFGetField(timg, TIFFTAG_BITSPERSAMPLE, &bps) == 0 ||
		(bps != 8 && bps != 16)) {
		fprintf(stderr, "bad bps %d\n", bps);
		goto end;
	}
	bytesperpixel = bps/8;

	if(TIFFGetField(timg, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0 || (spp != 1 && spp != 3)) {
		fprintf(stderr, "bad spp %d - assume 1\n", spp);
		spp = 1;
		//goto end;
	}
//else fprintf(stderr, "got spp %d\n", spp);

	if(TIFFGetField(timg, TIFFTAG_IMAGEWIDTH, &wid) == 0) {
		fprintf(stderr, "No wid\n");
		goto end;
	}

	if(TIFFGetField(timg, TIFFTAG_IMAGELENGTH, &ht) == 0) {
		fprintf(stderr, "No ht\n");
		goto end;
	}

	stsize = TIFFStripSize(timg);
	stmax = TIFFNumberOfStrips(timg);
	qread = 0;

	buffsz = stsize*bytesperpixel*(unsigned long)TIFFNumberOfStrips(timg);
	if((buff = (unsigned char *)malloc(buffsz)) == NULL) {
		fprintf(stderr, "malloc err (bytes = %ld)\n", buffsz);
		goto end;
	}

	for(strips = 0; strips < stmax; strips++) {
		if((got = TIFFReadEncodedStrip(timg, strips,
		    buff + qread, stsize)) == -1) {
			fprintf(stderr, "error strip %d\n", strips);
			goto end;
		}
		qread += got;
	}

	*widp = wid;
	*htp = ht;
	*bp = buff;
	//*bpp = spp;
	*bpp = bytesperpixel*spp;
//fprintf(stderr, "bpp %d = %d*%d\n", *bpp, bytesperpixel, spp);
	tiffout = 1;
end:	TIFFClose(timg);
/*
fprintf(stderr, "lt %d %d bpp %d 0x%lx\n", wid, ht, spp, buff);
{
int fd, nw;
char txt[100];
fprintf(stderr, "tgot\n");
if(spp == 3) {
fd = creat("tgot.ppm", 0666);
sprintf(txt, "P6\n%d %d\n255\n", wid, ht);
} else {
fd = creat("tgot.pgm", 0666);
sprintf(txt, "P5\n%d %d\n255\n", wid, ht);
}
nw = write(fd, txt, strlen(txt));
nw = write(fd, buff, wid*ht*spp);
close(fd);
}
*/
	return(0);
}

struct image *read_img(char *fname) {
	FILE *fp;
	char tmptxt[256];
	struct image *ip;
	int bpp = 1, wid, ht, range, c, nr;
	unsigned char *idata;
	if(strstr(fname, ".raw") || strstr(fname, ".RAW")
	    || strstr(fname, ".dat") || strstr(fname, ".DAT")) {
		int fd;
#define RWID	2560	// Davi temca2
#define RHT	2160
#define RBPP	1
#define	RHDR	0
fprintf(stderr, "Davi temca2 raw\n");
		ip = (struct image *)malloc(sizeof(struct image));
		ip->wid = RWID;
		ip->ht = RHT;
		ip->bpp = RBPP;
		ip->ydelta = RWID;
		ip->pp = (unsigned char *)malloc(ip->wid*ip->ht);
		fd = open(fname, 0);
fprintf(stderr, "raw fd %d\n", fd);
		if(RHDR)
			nr = read(fd, ip->pp, RHDR);
		nr = read(fd, ip->pp, ip->wid*ip->ht*ip->bpp);
		close(fd);
		return(ip);
	}
	if(strstr(fname, ".tif") || strstr(fname, ".TIF")) { // hits tiff too
//fprintf(stderr, "use ltiff\n");
		if(ltiff_reader(fname, &wid, &ht, &bpp, &idata) < 0)
			return(NULL);
//fprintf(stderr, "%d %d 0x%lx\n", wid, ht, idata);
		ip = (struct image *)malloc(sizeof(struct image));
		ip->wid = wid;
		ip->ht = ht;
		ip->bpp = bpp;
		ip->ydelta = wid;
		ip->pp = idata;
//fprintf(stderr, "TIFF wid %d ht %d bpp %d -> 0x%lx\n", wid, ht, bpp, idata);
		return(ip);
	}
	if(!(fp = fopen(fname, "rb"))) {
		//fprintf(stderr, "Can't open %s\n", fname);
		return(NULL);
	}
//fprintf(stderr, "fileno %d\n", fileno(fp));
	c = fgetc(fp);
	ungetc(c, fp);
	if(c == 0211) {		// assume its png
		int i, o, nr = 0, x, y, number_of_passes, rowbytes;
		png_bytep *row_ptr;
		uchar png_hdr[8];
		png_structp png_ptr;
		png_infop info_ptr;
		unsigned char *idata, color_type, bit_depth;
fprintf(stderr, "PNG\n");
		nr += fread(png_hdr, 1, 8, fp);
		if(png_sig_cmp(png_hdr, 0, 8))
			exit(1);
		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
fprintf(stderr, "png_ptr %p\n", png_ptr);
		png_set_sig_bytes(png_ptr, 8);
fprintf(stderr, "set done\n");
		info_ptr = png_create_info_struct(png_ptr);
		png_init_io(png_ptr, fp);
fprintf(stderr, "info_ptr %p\n", info_ptr);
		png_read_info(png_ptr, info_ptr);
fprintf(stderr, "png_ptr %p  info_ptr %p\n", png_ptr, info_ptr);
		ip = (struct image *)malloc(sizeof(struct image));
fprintf(stderr, "ip %p\n", ip);
		ip->wid = png_get_image_width(png_ptr, info_ptr);
		ip->ydelta = ip->wid;
		ip->ht = png_get_image_height(png_ptr, info_ptr);
		color_type = png_get_color_type(png_ptr, info_ptr);
		bit_depth = png_get_bit_depth(png_ptr, info_ptr);
		if(color_type == 0 && bit_depth == 8)
			ip->bpp = 1;
		else
			ip->bpp = 3;
fprintf(stderr, "%dx%d  %d %d\n", ip->wid, ip->ht, color_type, bit_depth);
		number_of_passes = png_set_interlace_handling(png_ptr);
fprintf(stderr, "%d passes\n", number_of_passes);
		png_read_update_info(png_ptr, info_ptr);
		row_ptr = (png_bytep *)malloc(sizeof(png_bytep)*ip->ht);
		if(bit_depth == 16)
			rowbytes = ip->wid*8;
		else
			rowbytes = ip->wid*4;
fprintf(stderr, "rowbytes %d\n", rowbytes);
		idata = malloc(ip->ht*ip->wid*ip->bpp);
		if(color_type == 0) {
			int rowlen = ip->wid*ip->bpp;
			for(y = 0; y < ip->ht; y++)
				row_ptr[y] = (png_byte *)&idata[y*rowlen];
			png_read_image(png_ptr, row_ptr);
		} else {
			int rowlen = ip->wid*ip->bpp;
			uchar *in, *out;
			for(y = 0; y < ip->ht; y++)
				row_ptr[y] = (png_byte *)malloc(rowbytes);
			png_read_image(png_ptr, row_ptr);
			for(y=0; y < ip->ht; y++) {
				in = row_ptr[y];
				out = &idata[y*rowlen];
				//bcopy(row_ptr[y], &idata[y*rowlen], ip->wid*ip->bpp);
				for(x = 0; x < ip->wid; x++) {
					*out++ = *in++; // R
					*out++ = *in++; // G
					*out++ = *in++; // B
					in++; // A
				}
			}
			for(y = 0; y < ip->ht; y++)
				free(row_ptr[y]);
		}
		free(row_ptr);
fprintf(stderr, "read done\n");
		ip->pp = idata;
	} else if(c == 0377) {		// assume its jpeg
		int i = 0, l = 0, row_stride;
		struct jpeg_decompress_struct cinfo;
		struct jpeg_error_mgr jerr;
		JSAMPARRAY row_buffer;
//fprintf(stderr, "jpg\n");
		//jread_n++;
		//jread_ti -= getticks();
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_decompress(&cinfo);
		jpeg_stdio_src(&cinfo, fp);
		jpeg_read_header(&cinfo, TRUE);
		jpeg_calc_output_dimensions(&cinfo);
		row_stride = cinfo.output_width * cinfo.output_components;
		wid = cinfo.output_width;
		ht = cinfo.output_height;
		//jread_pix += wid*ht;
		bpp = cinfo.output_components;
//fprintf(stderr, "whb %d %d %d\n", wid, ht, bpp);
		if((ip = newimage(wid, ht, bpp))) {
//fprintf(stderr, "ip 0x%lx pp 0x%lx whb %d %d %d\n", ip, ip->pp, ip->wid, ip->ht, ip->bpp);
			row_buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr)
				&cinfo, JPOOL_IMAGE, row_stride, 1);
			jpeg_start_decompress(&cinfo);
			for(i = 0; cinfo.output_scanline < cinfo.output_height; i++) {
				jpeg_read_scanlines(&cinfo, row_buffer, 1);
				memcpy(ip->pp + i*row_stride, row_buffer[0],
					row_stride);
			}
			jpeg_finish_decompress(&cinfo);
			jpeg_destroy_decompress(&cinfo);
			//jread_ti += getticks();
		} else {
			//jread_ti += getticks();
			return(NULL);
		}
	} else if(c == 'I') {	// assume its a tiff
		int i, o, nr = 0;
		ushort tiff_hdr[4], ntags, *tags;
		unsigned char *idata;
		nr += fread(tiff_hdr, 1, 8, fp);
		if(tiff_hdr[0] != 0x4949) {
			fprintf(stderr, "only II for now\n");
			exit(1);
		}
		o = *(int *)(&tiff_hdr[2]);
		if(o > 8) {
			idata = (unsigned char *)malloc(o-8);
			nr += fread(idata, 1, o-8, fp);
		}
		nr += fread(&ntags, 1, 2, fp);
		tags = (unsigned short *)malloc(12*ntags);
		nr += fread(tags, 1, 12*ntags, fp);
		for(i = 0; i < ntags; i++) {
//fprintf(stderr, "tag %d 0x%x\n", i, tags[6*i]);
			if(tags[6*i] == 0x100)
				wid = tags[6*i+4];
			else if(tags[6*i] == 0x101)
				ht = tags[6*i+4];
			else if(tags[6*i] == 0x111)
				o = (tags[6*i+3] << 16) + tags[6*i+4];
			/* else ignore all those other worthless tags!!! */
		}
//fprintf(stderr, "w %d h %d o %d nr %d img 0x%x\n", wid, ht, o, nr, idata);
		if(nr <= o) {
			idata = (unsigned char *)malloc(wid*ht);
			if(nr < o)
				nr += fread(idata, 1, o-nr, fp);
			nr += fread(idata, 1,  wid*ht, fp);
		} else if(o > 8) {
			// reshift data to start of malloc area
			bcopy(idata+o-8, idata, wid*ht);
		}
		free(tags);
		ip = (struct image *)malloc(sizeof(struct image));
		ip->wid = wid;
		ip->ht = ht;
		ip->bpp = bpp;
		ip->ydelta = wid;
		ip->pp = idata;
	} else { // ppm or pgm
		int ignore = fscanf(fp, "%[^\n] ", tmptxt);
		if(strcmp(tmptxt, "P5")) {
			if(strcmp(tmptxt, "P6")) {
				fprintf(stderr, "Not P5 or P6\n");
				return(NULL);
			}
			bpp = 3;
		}
//fprintf(stderr, "PPM/PGM %d\n", bpp);
		while((c = fgetc(fp)) == '#') {	// skip annoying XV comments
//fprintf(stderr, "SKIP CMT\n");
			int t;
			while((t = fgetc(fp)) != '\n')
				/*fprintf(stderr, "cmt %c\n", t)*/;
		}
//fprintf(stderr, "first ftell %d\n", ftell(fp));
		ungetc(c, fp);
//fprintf(stderr, "second ftell %d\n", ftell(fp));
		ignore = fscanf (fp, "%d%d", &wid, &ht);
		ignore = fscanf(fp, "%d", &range); // XXX itermittent "%d\n"
		if(range == 65535)
			bpp = 2;
//fprintf(stderr, "range %d bpp %d\n", range, bpp);
//fprintf(stderr, "third ftell %d\n", ftell(fp));
		c = fgetc(fp);
//fprintf(stderr, "w %d ht %d range %d c %d\n", wid, ht, range, c);
		if(range != 255 && range != 65535) {
			fprintf(stderr, "bad PPM range %d\n", range);
			free(ip);
			return(NULL);
		}
//fprintf(stderr, "fourth ftell %ld\n", ftell(fp));
		ip = newimage(wid, ht, bpp);
		if(!ip)
			return(NULL);
//fprintf(stderr, "fifth ftell %ld\n", ftell(fp));
		ignore = fread(ip->pp, 1, (long)ht*wid*bpp, fp);
/*
fprintf(stderr, "ignore %d  pp %p  %ld\n", ignore, ip->pp, (long)ht*wid*bpp);
{
static int beenhere = 0;
char hdr[50];
int fd, nw;
fprintf(stderr, "readcheck beenhere %d\n", beenhere);
if(beenhere++)
fd = creat("readcheck1.pgm", 0666);
else
fd = creat("readcheck0.pgm", 0666);
sprintf(hdr, "P%d\n%d %d\n%d\n", bpp==3?6:5, wid, ht, bpp==2?65535:255);
nw = write(fd, hdr, strlen(hdr));
nw = write(fd, ip->pp, (long)ht*wid*bpp);
close(fd);
}
*/
	}
	fclose(fp);
//#define	INGRID 64
#ifdef	INGRID
{
int x, y;
fprintf(stderr, "do grid.....\n");
//for(y = 0; y < wid*ht; y++) ip->pp[y] = 32;
for(y = 0; y < ht; y += INGRID) for(x = 0; x < wid; x++) ip->pp[y*wid+x] = 255;
for(y = 0; y < ht; y++) for(x = 0; x < wid; x += INGRID) ip->pp[y*wid+x] = 255;
}
#endif
	return(ip);
}

// www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
ushort tiff_hdr[] = { 0x4949, 42, 0x0, 0x0 };	// II 42 - offset to mfd
ushort tiff_mfd[] = { 7,			// ntags
0x100, 3, 1, 0, 0, 0, // width
0x101, 3, 1, 0, 0, 0, // height
0x102, 3, 1, 0, 8, 0, // bits per sample
0x106, 3, 1, 0, 1, 0, // photometric interpretation  0=white0, 1=black0,2=RGB
0x111, 4, 1, 0, 0, 0, // strip offsets
0x115, 4, 1, 0, 1, 0, // samples per pixel  1=gray 3=RGB
0x117, 4, 1, 0, 0, 0, // strip byte count; qiv+display warn but don't need
};

write_img(char *fname, struct image *ip) {
	char hdr[100];
	int fd, nw, rowbytes, color_type = 2;
	long qw, tqw = 0, want;
	unlink(fname);
	//fd = creat(fname, 0666);
fprintf(stderr, "write_img %d %s  %d %d %d\n",
fd, fname, ip->wid, ip->ht, ip->bpp);
	if(strstr(fname, "png") || strstr(fname, "PNG")) {
// see libpng-short-example.c
		png_structp png_ptr;
		png_infop info_ptr;
		png_bytep *row_ptr;
		FILE *fp = fopen(fname, "wb");
fprintf(stderr, "writePNG\n");
		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		info_ptr = png_create_info_struct(png_ptr);
		png_init_io(png_ptr, fp);
// PNG_COLOR_TYPE_GRAY_ALPHA
		if(ip->bpp == 1)
			color_type = PNG_COLOR_TYPE_GRAY;
		png_set_IHDR(png_ptr, info_ptr, ip->wid, ip->ht,
			8, color_type, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
		png_write_info(png_ptr, info_ptr);
		row_ptr = (png_bytep *)malloc(sizeof(png_bytep)*ip->ht);
		if(1 /* ip->bpp == 3*/) {
			int y, rowlen = ip->wid*ip->bpp;
			for(y = 0; y < ip->ht; y++)
				row_ptr[y] = (png_byte *)&ip->pp[y*rowlen];
			png_write_image(png_ptr, row_ptr);
		}
		free(row_ptr);
		png_write_end(png_ptr, NULL);
		fclose(fp);
	} else if(strstr(fname, "tif") || strstr(fname, "TIF")) {
		fd = creat(fname, 0666);
		*(int *)(&tiff_hdr[2]) = ip->wid*ip->ht*ip->bpp+8;
		nw = write(fd, tiff_hdr, 8);
		nw = write(fd, ip->pp, ip->wid*ip->ht*ip->bpp);
		tiff_mfd[5] = ip->wid;
		tiff_mfd[11] = ip->ht;
		if(ip->bpp == 3)
			*(int *)(&tiff_mfd[23]) = 2; // RGB
		*(int *)(&tiff_mfd[29]) = 8; // data offset
		*(int *)(&tiff_mfd[35]) = ip->bpp; // samples per pixel
		*(int *)(&tiff_mfd[41]) = ip->wid*ip->ht*ip->bpp; // byte count
		nw = write(fd, tiff_mfd, sizeof(tiff_mfd));
		close(fd);
	} else {
		fd = creat(fname, 0666);
		if(ip->bpp == 3)
			sprintf(hdr, "P6\n%d %d\n255\n", ip->wid, ip->ht);
		else if(ip->bpp == 2)
			sprintf(hdr, "P5\n%d %d\n65535\n", ip->wid, ip->ht);
		else
			sprintf(hdr, "P5\n%d %d\n255\n", ip->wid, ip->ht);
		nw = write(fd, hdr, strlen(hdr));
		//write(fd, ip->pp, (ip->ht*ip->wid*ip->bpp+7)/8);
		want = ip->ht*(unsigned long)ip->wid*ip->bpp;
		while(tqw < want) {
			qw = write(fd, ip->pp+tqw, want-tqw);
			if(qw <= 0) {
				fprintf(stderr, "write_img exit\n");
				exit(1);
			}
			tqw += qw;
		}
		close(fd);
	}
}