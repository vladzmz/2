#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tga.h"
#include "model.h"

typedef int Vector[3];
typedef double Mat4x4[4][4];
typedef double Mat4x1[4];
 
void swap(int *a, int *b);

int abs(int a);
double d_abs(double a);
int Round(double a);
int c_length(int a, int b, int c, int *s);
void product_vec3(Vec3 A, Vec3 B, Vec3 *W);
double v_length(Vec3 A);
double product_dot(Vec3 A, Vec3 B);
double intension (Vec3* v0, Vec3* v1, Vec3* v2, Vec3 light);
void normal_vec3(Vec3* A, double l);
void product_mat(Mat4x4 A, Mat4x1 B, Mat4x1* C);


void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);
void triangle(tgaImage *image, Model* model, Vec3 a, Vec3 b, Vec3 c, double I);

int main(int argc, char **argv)
{
    int rv = 0;
    if (argc < 4) {
        fprintf(stderr, "Usage: %s directory to obj or outfile \n ", argv[0]);
        return -1;
    }
    tgaImage * image = tgaNewImage(1000, 1000, RGB);
    Model *model = loadFromObj(argv[1]);
	int DiffuseMap = loadDiffuseMap(model, argv[2]); printf("%d", DiffuseMap);
    int i,j;
    double I;
	Vec3 light = { 2.0, 0.0, -1.0 };
    for(j = 0; j < model->nface; ++j) {
	Vec3 *v0 = getVertex(model, j, 0);
	Vec3 *v1 = getVertex(model, j, 1);
	Vec3 *v2 = getVertex(model, j, 2);
	Vec3 v0s, v1s, v2s;
	v0s[0]=((*v0)[0]+1)*500; 
	v1s[0]=((*v1)[0]+1)*500;		
	v2s[0]=((*v2)[0]+1)*500;
	v0s[1]=((*v0)[1]+1)*500; 
	v1s[1]=((*v1)[1]+1)*500;		
	v2s[1]=((*v2)[1]+1)*500;
	v0s[2]=(*v0)[2]; 
	v1s[2]=(*v1)[2];		
	v2s[2]=(*v2)[2];
// Pokazat intensivnjosti tekuwego poligona
	I=intension(v0, v1, v2, light);
	if(I>0){
	//line (image, v0s[0], v0s[1], v1s[0], v1s[1],tgaRGB(0,0,0));
	//line (image, v1s[0], v1s[1], v2s[0], v2s[1],tgaRGB(0,0,0));
	//line (image, v0s[0], v0s[1], v2s[0], v2s[1],tgaRGB(0,0,0));
	triangle(image, model, v0s, v1s, v2s, d_abs(I));
   } 
}
    tgaFlipVertically(image);
    if (-1 == tgaSaveToFile(image, argv[3])) {
        perror("tgaSateToFile");
        rv = -1;
    }
    tgaFreeImage(image); 
	freeModel(model);
    return rv;
}

void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color)
{  
    int flag = 0;
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    if (dx < dy) {
	    swap(&x0,&y0);
	    swap(&x1,&y1);
	    swap(&dx,&dy);
        flag = 1;
    }
    if (x0 > x1) {
        swap(&x0,&x1);
	    swap(&y0,&y1);
    }
    int e = 0;
    int de = 2*abs(y1 - y0);
    int i;
    int x = x0;
    int y = y0;	
    for (i = 0;i <= dx; ++i) {
	if (e >= dx) {
        e -= 2*dx;
	    if (y1 < y0) { 	
            y -= 1;
	    } else {
	        y += 1;
	    }
        }
	if (flag) {
	    tgaSetPixel(image, y, x, color);
	} else {
	    tgaSetPixel(image, x, y, color);
	}
	x += 1;
	e += de;
    }
}

void triangle(tgaImage *image, Model* model, Vec3 a, Vec3 b, Vec3 c, double I) {
    int sh,t;
    Vec3 A,B;
    int i,j,height,half;
    double alpha,beta;
    if (a[1]>b[1]) 
	{t=a[0];a[0]=b[0];b[0]=t;t=a[1];a[1]=b[1];b[1]=t;}
    if (a[1]>c[1]) 
	{t=a[0];a[0]=c[0];c[0]=t;t=a[1];a[1]=c[1];c[1]=t;}
    if (b[1]>c[1]) 
	{t=b[0];b[0]=c[0];c[0]=t;t=b[1];b[1]=c[1];c[1]=t;}
    height=c[1]-a[1];
    for (i=0; i<height; i++) {
// flag polovin6I
        if ((i>b[1]-a[1])||(b[1]==a[1])) half=1; else half=0;
        sh = half ? c[1]-b[1] : b[1]-a[1];
        alpha = (float)i/height;
        beta  = (float)(i-(half ? b[1]-a[1] : 0))/sh; 
        A[0] = a[0] + (c[0]-a[0])*alpha;
	A[1] = a[1] + (c[1]-a[1])*alpha;
        B[0] = half ? b[0] + (c[0]-b[0])*beta : a[0] + (b[0]-a[0])*beta;
	B[1] = half ? b[1] + (c[1]-b[1])*beta : a[1] + (b[1]-a[1])*beta;
        if (A[0]>B[0]) 
		{t=A[0];A[0]=B[0];B[0]=t;t=A[1];A[1]=B[1];B[1]=t;};
        for ( j=A[0]; j<=B[0]; j++) {
		tgaSetPixel(image, j, a[1]+i, tgaRGB(I*255, I*255, I*255) ); 
        	}
              }
    }

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}


int c_length(int A, int B, int C, int *s) {
    int a = A;
    int b = B;
    int c = C;
    if (a >= b) {
        swap(&a, &b);
    }
    if (a >= c) {
        swap(&a, &c);
    }
    if (b >= c) {
        swap(&b, &c);
    }
    *s = a;
    return c ;
}

int abs(int a) {
    return (a > 0) ? a : -a;
}

double d_abs(double a) {
    return (a > 0) ? a : -a;
}

int Round(double a) {
    int b = a;
    if ((a - b) >= 0.5) {
        return (b + 1);
    } else {
	return b;
    }
}

void product_vec3(Vec3 A, Vec3 B, Vec3* W) {
    (*W)[0] = A[1]*B[2] - A[2]*B[1];
    (*W)[1] = B[0]*A[2] - A[0]*B[2];
    (*W)[2] = A[0]*B[1] - A[1]*B[0];
}

double product_dot(Vec3 A, Vec3 B) {
    double dot = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return dot;
} 

double v_length(Vec3 A) {
    double length = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    return length;
}
double intension (Vec3* v0, Vec3* v1, Vec3* v2, Vec3 light) {
    Vec3 AB, AC;
    AB[0] = (*v1)[0] - (*v0)[0];
    AB[1] = (*v1)[1] - (*v0)[1];
    AB[2] = (*v1)[2] - (*v0)[2];
    AC[0] = (*v2)[0] - (*v0)[0];
    AC[1] = (*v2)[1] - (*v0)[1];
    AC[2] = (*v2)[2] - (*v0)[2];
    Vec3 W;
    product_vec3(AB, AC, &W);
    W[0] = -W[0];
    W[1] = -W[1];
    W[2] = -W[2];
    normal_vec3(&W, v_length(W));
    return product_dot(light, W); 
}

void normal_vec3(Vec3* A, double l) {
    (*A)[0] = (*A)[0]/l;
    (*A)[1] = (*A)[1]/l;
    (*A)[2] = (*A)[2]/l;
}

void product_mat(Mat4x4 A, Mat4x1 B, Mat4x1* C) {
    int i,j;
    Mat4x1 V = {0.0, 0.0, 0.0, 0.0};
    for(i=0; i < 4; ++i) {
        for(j=0; j < 4; ++j) {
            V[i] = V[i] + A[i][j]*B[j];
        }
        (*C)[i] = V[i];
    }
}


 
