#include <stdio.h>
#include "tga.h"
#include "model.h"
void swap(int *a, int *b);
int abs(int a);

/*
* Using Digital Differential Analyzer algorihm
* to draw interval connecting (x0, y0) with (x1, y1)
* on image using color
*/
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);

int main(int argc, char **argv)
{	tgaImage * image = tgaNewImage(800, 800, RGB);
	Model*mod=loadFromObj("cat.obj");
    	int nf=mod->nface;
	int nv=mod->nvert;
	int i;int j;
	for (i=0;i<nf;++i){;
		for(j=0;j<3;++j){
	
	Vec3 *v_0=getVertex(mod,i,j);
	double x=((*v_0)[0]+1.)*400;
	double y=((*v_0)[1]+1.)*400;
	
	Vec3 *v_1=getVertex(mod,i,(j+1)%3);
	double x1=((*v_1)[0]+1.)*400;
	double y1=((*v_1)[1]+1.)*400;

	line(image, x,y,x1,y1, tgaRGB(0, 255, 0));
}
}
tgaFlipVertically(image);
int rv = 0;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }

    if (-1 == tgaSaveToFile(image, argv[1])) {
        perror("tgaSateToFile");
        rv = -1;
    }

    tgaFreeImage(image);    
    return rv;
}

void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color)
{
    int steep = 0;
    if (abs(y1 - y0) > abs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }

    int x;
    double y;
    double k = ((double)(y1 - y0))/(x1 - x0);
    for (x = x0, y = y0; x <= x1; ++x, y += k) {
        if (steep) {
            tgaSetPixel(image, y, x, color);
        } else {
            tgaSetPixel(image, x, y, color);
        }
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int abs(int a) {
    return (a >= 0) ? a : -a;
}
