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
int main()
{
    std::cout << "Hello World!\n"; 
}


