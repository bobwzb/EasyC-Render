// EasyC++Render.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <windows.h>
#include <tchar.h>

//define some func for math
typedef unsigned int IUINT32;
typedef struct { float m[4][4]; } matrix;
typedef struct { float x, y, z, w; } vector;
typedef vector point;
//find the mid of three value
int CMID(int x, int min, int max) {
	return (x < min) ? min : ((x > max) ? max : x);
}
//caculate the interpolation
float interp(float x1, float x2, float t) {
	return x1 + (x2 - x1)*t;
}
//caculate the length of vector
float vector_length(vector v) {
	float sq = v.x*v.x + v.y*v.y + v.z*v.z;
	return (float)sqrt(sq);
}
//z=x+y
vector vector_add( vector x, vector y) {
	vector z;
	z.x = x.x + y.x;
	z.y = x.y + y.y;
	z.z = x.z + y.z;
	z.w = 1.0;
	return z;
}
//z=x-y
vector vector_sub(vector x, vector y) {
	vector z;
	z.x = x.x - y.x;
	z.y = x.y - y.y;
	z.z = x.z - y.z;
	z.w = 1.0;
	return z;
}
//dot product
float vector_dot(vector z, vector x) {
	return z.x*x.x + z.y*x.y + z.z*x.z;
}
//cross product
vector vector_crossproduct(vector x, vector y) {
	vector tmp;
	tmp.x = x.y*y.z - x.z*y.y;
	tmp.y = x.z*y.x - x.x*y.z;
	tmp.z = x.x*y.y - x.y*y.x;
	tmp.w = 1.0f;
	return tmp;
}
// vector interpolation
vector vector_interp(vector x, vector y, float t) {
	vector tmp;
	tmp.x = interp(x.x, y.x, t);
	tmp.y = interp(x.y, y.y, t);
	tmp.z = interp(x.z, y.z, t);
	tmp.w = 1.0f;
	return tmp;
}
//vector normalize
void vector_normalize(vector& x) {
	float length = vector_length(x);
	if (length != 0) {
		float tmp = 1.0f / length;
		x.x *= tmp;
		x.y *= tmp;
		x.z *= tmp;
	}
}
//matrix add
matrix matrix_add(matrix a,matrix b) {
	matrix c;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			c.m[i][j] = a.m[i][j] + b.m[i][j];
		}
	}
	return c;
}
//matrix minus
matrix matrix_sub(matrix a, matrix b) {
	matrix c;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			c.m[i][j] = a.m[i][j] - b.m[i][j];
		}
	}
	return c;
}
//matrix multiply
matrix matrix_mul(matrix a, matrix b) {
	matrix c;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			c.m[i][j] = (a.m[j][0] * b.m[0][i]) +
						(a.m[j][1] * b.m[1][i]) +
						(a.m[j][2] * b.m[2][i]) +
						(a.m[j][4] * b.m[4][i]);
		}
	}
	return c;
}
// matrix scale
void matrix_scale(matrix& a, float b) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			a.m[i][j] = a.m[i][j]*b;
		}
	}
}
//matrix apply
vector matrix_apply(vector a, matrix b) {
	float X = a.x, Y = a.y, Z = a.z, W = a.w;
	vector tmp;
	tmp.x = X * b.m[0][0] + Y * b.m[1][0] + Z * b.m[2][0] + W * b.m[3][0];
	tmp.y = X * b.m[0][1] + Y * b.m[1][1] + Z * b.m[2][1] + W * b.m[3][1];
	tmp.x = X * b.m[0][2] + Y * b.m[1][2] + Z * b.m[2][2] + W * b.m[3][2];
	tmp.x = X * b.m[0][3] + Y * b.m[1][3] + Z * b.m[2][3] + W * b.m[3][3];
	return tmp;
}
//identity matrix
matrix matrix_identity() {
	matrix m;
	m.m[0][0] = m.m[1][1] = m.m[2][2] = m.m[3][3] = 1.0f;
	m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
	return m;
}
//set matrix to zero
void matrix_zero(matrix& m) {
	m.m[0][0] = m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][1] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][2] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = m.m[3][3] = 0.0f;
}
//matrix translate
matrix matrix_translate(float x,float y,float z) {
	matrix m = matrix_identity();
	m.m[3][0] = x;
	m.m[3][1] = y;
	m.m[3][2] = z;
	return m;
}
//matrix set scale
void matrix_set_scale(matrix& m, float x, float y, float z) {
	m.m[0][0] = x;
	m.m[1][1] = y;
	m.m[2][2] = z;
}
//matrix rotation
matrix matrix_set_rotate(float x, float y, float z, float theta) {
	matrix m;
	float qsin = (float)sin(theta * 0.5f);
	float qcos = (float)cos(theta * 0.5f);
	vector vec = { x, y, z, 1.0f };
	float w = qcos;
	vector_normalize(vec);
	x = vec.x * qsin;
	y = vec.y * qsin;
	z = vec.z * qsin;
	m.m[0][0] = 1 - 2 * y * y - 2 * z * z;
	m.m[1][0] = 2 * x * y - 2 * w * z;
	m.m[2][0] = 2 * x * z + 2 * w * y;
	m.m[0][1] = 2 * x * y + 2 * w * z;
	m.m[1][1] = 1 - 2 * x * x - 2 * z * z;
	m.m[2][1] = 2 * y * z - 2 * w * x;
	m.m[0][2] = 2 * x * z - 2 * w * y;
	m.m[1][2] = 2 * y * z + 2 * w * x;
	m.m[2][2] = 1 - 2 * x * x - 2 * y * y;
	m.m[0][3] = m.m[1][3] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
	m.m[3][3] = 1.0f;
	return m;
}
// set camera
matrix matrix_lookat(vector eye, vector at, vector up) {
	matrix m;
	vector xaxis, yaxis, zaxis;

	zaxis=vector_sub(at, eye);
	vector_normalize(zaxis);
	xaxis=vector_crossproduct( up, zaxis);
	vector_normalize(xaxis);
	yaxis=vector_crossproduct(zaxis, xaxis);

	m.m[0][0] = xaxis.x;
	m.m[1][0] = xaxis.y;
	m.m[2][0] = xaxis.z;
	m.m[3][0] = -vector_dot(xaxis, eye);

	m.m[0][1] = yaxis.x;
	m.m[1][1] = yaxis.y;
	m.m[2][1] = yaxis.z;
	m.m[3][1] = -vector_dot(yaxis, eye);

	m.m[0][2] = zaxis.x;
	m.m[1][2] = zaxis.y;
	m.m[2][2] = zaxis.z;
	m.m[3][2] = -vector_dot(zaxis, eye);

	m.m[0][3] = m.m[1][3] = m.m[2][3] = 0.0f;
	m.m[3][3] = 1.0f;
	return m;
}
//matrix perspective fov
matrix matrix_set_perspective(float fovy, float aspect, float zn, float zf) {
	matrix m;
	matrix_zero(m);
	float fax = 1.0f / (float)tan(fovy * 0.5f);
	m.m[0][0] = (float)(fax / aspect);
	m.m[1][1] = (float)(fax);
	m.m[2][2] = zf / (zf - zn);
	m.m[3][2] = -zn * zf / (zf - zn);
	m.m[2][3] = 1;
	return m;
}
typedef struct {
	matrix world;//world translation
	matrix view; //camera tranlation
	matrix projection;// projection translation
	matrix trans; //transform = world * view * projection
	float w, h; //size of the screen
} transform;
//transfrom update
void transform_update(transform& t) {
	matrix m;
	m = matrix_mul(t.world, t.view);
	t.trans= matrix_mul(m, t.projection);
}
//initalize tranform
transform transform_init(int width, int height) {
	transform t;
	float aspect = (float)width / ((float)height);
	t.world=matrix_identity();
	t.view = matrix_identity();
	t.projection=matrix_set_perspective(3.1415926f * 0.5f, aspect, 1.0f, 500.0f);
	t.w = (float)width;
	t.h = (float)height;
	transform_update(t);
	return t;
}
// project vector x
vector transform_apply(transform t, vector x) {
	vector y = matrix_apply(x, t.trans);
	return y;
}
//check cvv for projection
int transform_check_cvv(vector v) {
	float w = v.w;
	int check = 0;
	if (v.z < 0.0f) check |= 1;
	if (v.z > w) check |= 2;
	if (v.x < -w) check |= 4;
	if (v.x > w) check |= 8;
	if (v.y < -w) check |= 16;
	if (v.y > w) check |= 32;
	return check;
}
//get sreen position
vector transform_homogenize(transform t,  vector x) {
	vector y;
	float rhw = 1.0f / x.w;
	y.x = (x.x * rhw + 1.0f) * t.w * 0.5f;
	y.y = (1.0f - x.y * rhw) * t.h * 0.5f;
	y.z = x.z * rhw;
	y.w = 1.0f;
	return y;
}
//Geometric calculation
typedef struct { float r, g, b; } color;
typedef struct { float u, v; } texcoord;
typedef struct { point pos; texcoord tc; color color; float rhw; } vertex;

typedef struct { vertex v, v1, v2; } edge;
typedef struct { float top, bottom; edge left, right; } trapezoid;
typedef struct { vertex v, step; int x, y, w; } scanline;
//Reciprocal of Homogeneous
void vertex_rhw_init(vertex& v) {
	float rhw = 1.0f / v.pos.w;
	v.rhw = rhw;
	v.tc.u *= rhw;
	v.tc.v *= rhw;
	v.color.r *= rhw;
	v.color.g *= rhw;
	v.color.b *= rhw;
}
// vertex interp
vertex vertex_interp(vertex x1, const vertex x2, float t) {
	vertex y;
	y.pos=vector_interp(x1.pos, x2.pos, t);
	y.tc.u = interp(x1.tc.u, x2.tc.u, t);
	y.tc.v = interp(x1.tc.v, x2.tc.v, t);
	y.color.r = interp(x1.color.r, x2.color.r, t);
	y.color.g = interp(x1.color.g, x2.color.g, t);
	y.color.b = interp(x1.color.b, x2.color.b, t);
	y.rhw = interp(x1.rhw, x2.rhw, t);
	return y;
}
//vertex division
vertex vertex_division(vertex x1, vertex x2, float w) {
	vertex y;
	float inv = 1.0f / w;
	y.pos.x = (x2.pos.x - x1.pos.x) * inv;
	y.pos.y = (x2.pos.y - x1.pos.y) * inv;
	y.pos.z = (x2.pos.z - x1.pos.z) * inv;
	y.pos.w = (x2.pos.w - x1.pos.w) * inv;
	y.tc.u = (x2.tc.u - x1.tc.u) * inv;
	y.tc.v = (x2.tc.v - x1.tc.v) * inv;
	y.color.r = (x2.color.r - x1.color.r) * inv;
	y.color.g = (x2.color.g - x1.color.g) * inv;
	y.color.b = (x2.color.b - x1.color.b) * inv;
	y.rhw = (x2.rhw - x1.rhw) * inv;
	return y;
}
//vertex add
void vertex_add(vertex& y, vertex x) {
	y.pos.x += x.pos.x;
	y.pos.y += x.pos.y;
	y.pos.z += x.pos.z;
	y.pos.w += x.pos.w;
	y.rhw += x.rhw;
	y.tc.u += x.tc.u;
	y.tc.v += x.tc.v;
	y.color.r += x.color.r;
	y.color.g += x.color.g;
	y.color.b += x.color.b;
}
// create trapezoid by triangle
int trapezoid_init_triangle(trapezoid trap[], vertex p1,
	vertex p2, vertex p3) {
	vertex p;
	float k, x;

	if (p1.pos.y > p2.pos.y) p = p1, p1 = p2, p2 = p;
	if (p1.pos.y > p3.pos.y) p = p1, p1 = p3, p3 = p;
	if (p2.pos.y > p3.pos.y) p = p2, p2 = p3, p3 = p;
	if (p1.pos.y == p2.pos.y && p1.pos.y == p3.pos.y) return 0;
	if (p1.pos.x == p2.pos.x && p1.pos.x == p3.pos.x) return 0;

	if (p1.pos.y == p2.pos.y) {	// triangle down
		if (p1.pos.x > p2.pos.x) p = p1, p1 = p2, p2 = p;
		trap[0].top = p1.pos.y;
		trap[0].bottom = p3.pos.y;
		trap[0].left.v1 = p1;
		trap[0].left.v2 = p3;
		trap[0].right.v1 = p2;
		trap[0].right.v2 = p3;
		return (trap[0].top < trap[0].bottom) ? 1 : 0;
	}

	if (p2.pos.y == p3.pos.y) {	// triangle up
		if (p2.pos.x > p3.pos.x) p = p2, p2 = p3, p3 = p;
		trap[0].top = p1.pos.y;
		trap[0].bottom = p3.pos.y;
		trap[0].left.v1 = p1;
		trap[0].left.v2 = p2;
		trap[0].right.v1 = p1;
		trap[0].right.v2 = p3;
		return (trap[0].top < trap[0].bottom) ? 1 : 0;
	}

	trap[0].top = p1.pos.y;
	trap[0].bottom = p2.pos.y;
	trap[1].top = p2.pos.y;
	trap[1].bottom = p3.pos.y;

	k = (p3.pos.y - p1.pos.y) / (p2.pos.y - p1.pos.y);
	x = p1.pos.x + (p2.pos.x - p1.pos.x) * k;

	if (x <= p3.pos.x) {		// triangle left
		trap[0].left.v1 = p1;
		trap[0].left.v2 = p2;
		trap[0].right.v1 = p1;
		trap[0].right.v2 = p3;
		trap[1].left.v1 = p2;
		trap[1].left.v2 = p3;
		trap[1].right.v1 = p1;
		trap[1].right.v2 = p3;
	}
	else {					// triangle right
		trap[0].left.v1 = p1;
		trap[0].left.v2 = p3;
		trap[0].right.v1 = p1;
		trap[0].right.v2 = p2;
		trap[1].left.v1 = p1;
		trap[1].left.v2 = p3;
		trap[1].right.v1 = p2;
		trap[1].right.v2 = p3;
	}

	return 2;
}
// caculate the left and right vertex of the trapezoid
void trapezoid_edge_interp(trapezoid& trap, float y) {
	float s1 = trap.left.v2.pos.y - trap.left.v1.pos.y;
	float s2 = trap.right.v2.pos.y - trap.right.v1.pos.y;
	float t1 = (y - trap.left.v1.pos.y) / s1;
	float t2 = (y - trap.right.v1.pos.y) / s2;
	trap.left.v=vertex_interp(trap.left.v1, trap.left.v2, t1);
	trap.right.v=vertex_interp(trap.right.v1, trap.right.v2, t2);
}

// caculate the start point and step of the scanline
void trapezoid_init_scan_line(trapezoid trap, scanline& scanline, int y) {
	float width = trap.right.v.pos.x - trap.left.v.pos.x;
	scanline.x = (int)(trap.left.v.pos.x + 0.5f);
	scanline.w = (int)(trap.right.v.pos.x + 0.5f) - scanline.x;
	scanline.y = y;
	scanline.v = trap.left.v;
	if (trap.left.v.pos.x >= trap.right.v.pos.x) scanline.w = 0;
	scanline.step=vertex_division(trap.left.v, trap.right.v, width);
}
//render screen
typedef struct {
	transform transform;     
	int width;                  
	int height;                 
	IUINT32 **framebuffer;      
	float **zbuffer;            
	IUINT32 **texture;          
	int tex_width;              
	int tex_height;             
	float max_u;                // tex_width - 1
	float max_v;                // tex_height - 1
	int render_state;           
	IUINT32 background;        
	IUINT32 foreground;         
}	device;
#define RENDER_STATE_WIREFRAME      1		
#define RENDER_STATE_TEXTURE        2		
#define RENDER_STATE_COLOR          4
//initialize the device
void device_init(device &device, int width, int height, void *fb) {
	int need = sizeof(void*) * (height * 2 + 1024) + width * height * 8;
	char *ptr = (char*)malloc(need + 64);
	char *framebuf, *zbuf;
	int j;
	assert(ptr);
	device.framebuffer = (IUINT32**)ptr;
	device.zbuffer = (float**)(ptr + sizeof(void*) * height);
	ptr += sizeof(void*) * height * 2;
	device.texture = (IUINT32**)ptr;
	ptr += sizeof(void*) * 1024;
	framebuf = (char*)ptr;
	zbuf = (char*)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL) framebuf = (char*)fb;
	for (j = 0; j < height; j++) {
		device.framebuffer[j] = (IUINT32*)(framebuf + width * 4 * j);
		device.zbuffer[j] = (float*)(zbuf + width * 4 * j);
	}
	device.texture[0] = (IUINT32*)ptr;
	device.texture[1] = (IUINT32*)(ptr + 16);
	memset(device.texture[0], 0, 64);
	device.tex_width = 2;
	device.tex_height = 2;
	device.max_u = 1.0f;
	device.max_v = 1.0f;
	device.width = width;
	device.height = height;
	device.background = 0xc0c0c0;
	device.foreground = 0;
	device.transform=transform_init(width, height);
	device.render_state = RENDER_STATE_WIREFRAME;
}
// delete device
void device_destroy(device device) {
	if (device.framebuffer)
		free(device.framebuffer);
	device.framebuffer = NULL;
	device.zbuffer = NULL;
	device.texture = NULL;
}
// set the texture
void device_set_texture(device& device, void *bits, long pitch, int w, int h) {
	char *ptr = (char*)bits;
	int j;
	assert(w <= 1024 && h <= 1024);
	for (j = 0; j < h; ptr += pitch, j++) 	// recaculate the pointer
		device.texture[j] = (IUINT32*)ptr;
	device.tex_width = w;
	device.tex_height = h;
	device.max_u = (float)(w - 1);
	device.max_v = (float)(h - 1);
}
// draw point
void device_pixel(device& device, int x, int y, IUINT32 color) {
	if (((IUINT32)x) < (IUINT32)device.width && ((IUINT32)y) < (IUINT32)device.height) {
		device.framebuffer[y][x] = color;
	}
}
// clear framebuffer and zbuffer
void device_clear(device& device, int mode) {
	int y, x, height = device.height;
	for (y = 0; y < device.height; y++) {
		IUINT32 *dst = device.framebuffer[y];
		IUINT32 cc = (height - 1 - y) * 230 / (height - 1);
		cc = (cc << 16) | (cc << 8) | cc;
		if (mode == 0) cc = device.background;
		for (x = device.width; x > 0; dst++, x--) dst[0] = cc;
	}
	for (y = 0; y < device.height; y++) {
		float *dst = device.zbuffer[y];
		for (x = device.width; x > 0; dst++, x--) dst[0] = 0.0f;
	}
}
// draw the line
void device_draw_line(device& device, int x1, int y1, int x2, int y2, IUINT32 c) {
	int x, y, rem = 0;
	if (x1 == x2 && y1 == y2) {
		device_pixel(device, x1, y1, c);
	}
	else if (x1 == x2) {
		int inc = (y1 <= y2) ? 1 : -1;
		for (y = y1; y != y2; y += inc) device_pixel(device, x1, y, c);
		device_pixel(device, x2, y2, c);
	}
	else if (y1 == y2) {
		int inc = (x1 <= x2) ? 1 : -1;
		for (x = x1; x != x2; x += inc) device_pixel(device, x, y1, c);
		device_pixel(device, x2, y2, c);
	}
	else {
		int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
		int dy = (y1 < y2) ? y2 - y1 : y1 - y2;
		if (dx >= dy) {
			if (x2 < x1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; x <= x2; x++) {
				device_pixel(device, x, y, c);
				rem += dy;
				if (rem >= dx) {
					rem -= dx;
					y += (y2 >= y1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
		else {
			if (y2 < y1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; y <= y2; y++) {
				device_pixel(device, x, y, c);
				rem += dx;
				if (rem >= dy) {
					rem -= dy;
					x += (x2 >= x1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
	}
}
//read texture by position
IUINT32 device_texture_read(device device, float u, float v) {
	int x, y;
	u = u * device.max_u;
	v = v * device.max_v;
	x = (int)(u + 0.5f);
	y = (int)(v + 0.5f);
	x = CMID(x, 0, device.tex_width - 1);
	y = CMID(y, 0, device.tex_height - 1);
	return device.texture[y][x];
}
int main()
{
    std::cout << "Hello World!\n"; 
}


