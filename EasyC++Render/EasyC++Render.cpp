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
vector vector_add(vector z, vector x, vector y) {
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
matrix matrix_set_scale(matrix& m, float x, float y, float z) {
	matrix m;
	m.m[0][0] = x;
	m.m[1][1] = y;
	m.m[2][2] = z;
	return m;
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
int main()
{
    std::cout << "Hello World!\n"; 
}


