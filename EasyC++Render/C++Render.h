#pragma once
#ifndef _CPP_RENDER_
#define _CPP_RENDER_
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <map>
#include <initializer_list>
#include <stdexcept>
#include <functional>
#include <ostream>
#include <sstream>
#include <iostream>

//vector: N is the dimension, T is datatype
template<size_t N, typename T> struct Vector {
	int* hp;
	T m[N];
	inline Vector() {
		for (size_t i = 0; i < N; i++)
			m[i] = T();
	}
	inline Vector(const T *ptr) { 
		for (size_t i = 0; i < N; i++) 
			m[i] = ptr[i]; 
	}
	inline Vector(const Vector<N, T> &u) { 
		for (size_t i = 0; i < N; i++) 
			m[i] = u.m[i];
	}
	inline Vector(const std::initializer_list<T> &u) {
		auto it = u.begin(); 
		for (size_t i = 0; i < N; i++) 
			m[i] = *it++;
	}
	inline T& operator[] (size_t i) { 
		assert(i < N, "Out of bound");
		return m[i]; 
	}
	inline void load(const T *ptr) { 
		for (size_t i = 0; i < N; i++) 
			m[i] = ptr[i]; 
	}
	inline void save(T *ptr) { 
		for (size_t i = 0; i < N; i++) 
			ptr[i] = m[i]; 
	}
};

//specialize 2D vector
template <typename T> struct Vector<2, T> {
	union {
		struct { T x, y; };    
		struct { T u, v; };    
		T m[2];                
	};
	inline Vector() : x(T()), y(T()) {}
	inline Vector(T X, T Y) : x(X), y(Y) {}
	inline Vector(const Vector<2, T> &u) : x(u.x), y(u.y) {}
	inline Vector(const T *ptr) : x(ptr[0]), y(ptr[1]) {}
	inline const T& operator[] (size_t i) const { 
		assert(i < 2, "Out of bound");
		return m[i]; 
	}
	inline T& operator[] (size_t i) { 
		assert(i < 2, "Out of bound");
		return m[i]; 
	}
	inline void load(const T *ptr) { 
		for (size_t i = 0; i < 2; i++) 
			m[i] = ptr[i]; 
	}
	inline void save(T *ptr) { 
		for (size_t i = 0; i < 2; i++) 
			ptr[i] = m[i]; 
	}
	inline Vector<2, T> xy() const { 
		return *this; 
	}
	inline Vector<3, T> xy1() const { 
		return Vector<3, T>(x, y, 1); 
	}
	inline Vector<4, T> xy11() const { 
		return Vector<4, T>(x, y, 1, 1); 
	}
};

//specialize 3D vector
template <typename T> struct Vector<3, T> {
	union {
		struct { T x, y, z; };    
		struct { T r, g, b; };    
		T m[3];                   
	};
	inline Vector() : x(T()), y(T()), z(T()) {}
	inline Vector(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
	inline Vector(const Vector<3, T> &u) : x(u.x), y(u.y), z(u.z) {}
	inline Vector(const T *ptr) : x(ptr[0]), y(ptr[1]), z(ptr[2]) {}
	inline const T& operator[] (size_t i) const { 
		assert(i < 3, "Out of bound");
		return m[i]; 
	}
	inline T& operator[] (size_t i) { 
		assert(i < 3, "Out of bound");
	    return m[i];
	}
	inline void load(const T *ptr) { 
		for (size_t i = 0; i < 3; i++) 
			m[i] = ptr[i]; 
	}
	inline void save(T *ptr) { 
		for (size_t i = 0; i < 3; i++) 
			ptr[i] = m[i]; 
	}
	inline Vector<2, T> xy() const { 
		return Vector<2, T>(x, y); 
	}
	inline Vector<3, T> xyz() const { 
		return *this; 
	}
	inline Vector<4, T> xyz1() const { 
		return Vector<4, T>(x, y, z, 1); 
	}
};

//specialize 4D vector
template <typename T> struct Vector<4, T> {
	union {
		struct { T x, y, z, w; };    
		struct { T r, g, b, a; };    
		T m[4];                      
	};
	inline Vector() : x(T()), y(T()), z(T()), w(T()) {}
	inline Vector(T X, T Y, T Z, T W) : x(X), y(Y), z(Z), w(W) {}
	inline Vector(const Vector<4, T> &u) : x(u.x), y(u.y), z(u.z), w(u.w) {}
	inline Vector(const T *ptr) : x(ptr[0]), y(ptr[1]), z(ptr[2]), w(ptr[3]) {}
	inline const T& operator[] (size_t i) const { 
		assert(i < 4, "Out of bound");
		return m[i]; 
	}
	inline T& operator[] (size_t i) { 
		assert(i < 4, "Out of bound"); 
		return m[i]; }
	inline void load(const T *ptr) { 
		for (size_t i = 0; i < 4; i++) 
			m[i] = ptr[i]; 
	}
	inline void save(T *ptr) { 
		for (size_t i = 0; i < 4; i++) 
			ptr[i] = m[i]; 
	}
	inline Vector<2, T> xy() const { 
		return Vector<2, T>(x, y); 
	}
	inline Vector<3, T> xyz() const { 
		return Vector<3, T>(x, y, z); 
	}
	inline Vector<4, T> xyzw() const { 
		return *this; 
	}
};
// = (+a)
template <size_t N, typename T>
inline Vector<N, T> operator + (const Vector<N, T>& a) {
	return a;
}

// = (-a)
template <size_t N, typename T>
inline Vector<N, T> operator - (const Vector<N, T>& a) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) b[i] = -a[i];
	return b;
}

// = (a == b) ? true : false
template <size_t N, typename T>
inline bool operator == (const Vector<N, T>& a, const Vector<N, T>& b) {
	for (size_t i = 0; i < N; i++) if (a[i] != b[i]) return false;
	return true;
}

// = (a != b)? true : false
template <size_t N, typename T>
inline bool operator != (const Vector<N, T>& a, const Vector<N, T>& b) {
	return !(a == b);
}

// = a + b
template <size_t N, typename T>
inline Vector<N, T> operator + (const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = a[i] + b[i];
	return c;
}

// = a - b
template <size_t N, typename T>
inline Vector<N, T> operator - (const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = a[i] - b[i];
	return c;
}

// = a * b，not dot product, use for color multiple
template <size_t N, typename T>
inline Vector<N, T> operator * (const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = a[i] * b[i];
	return c;
}

// = a / b, divide each element
template <size_t N, typename T>
inline Vector<N, T> operator / (const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = a[i] / b[i];
	return c;
}

// = a * x
template <size_t N, typename T>
inline Vector<N, T> operator * (const Vector<N, T>& a, T x) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) b[i] = a[i] * x;
	return b;
}

// = x * a
template <size_t N, typename T>
inline Vector<N, T> operator * (T x, const Vector<N, T>& a) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) b[i] = a[i] * x;
	return b;
}

// = a / x
template <size_t N, typename T>
inline Vector<N, T> operator / (const Vector<N, T>& a, T x) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) b[i] = a[i] / x;
	return b;
}

// = x / a
template <size_t N, typename T>
inline Vector<N, T> operator / (T x, const Vector<N, T>& a) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) b[i] = x / a[i];
	return b;
}

// a += b
template <size_t N, typename T>
inline Vector<N, T>& operator += (Vector<N, T>& a, const Vector<N, T>& b) {
	for (size_t i = 0; i < N; i++) a[i] += b[i];
	return a;
}

// a -= b
template <size_t N, typename T>
inline Vector<N, T>& operator -= (Vector<N, T>& a, const Vector<N, T>& b) {
	for (size_t i = 0; i < N; i++) a[i] -= b[i];
	return a;
}

// a *= b
template <size_t N, typename T>
inline Vector<N, T>& operator *= (Vector<N, T>& a, const Vector<N, T>& b) {
	for (size_t i = 0; i < N; i++) a[i] *= b[i];
	return a;
}

// a /= b
template <size_t N, typename T>
inline Vector<N, T>& operator /= (Vector<N, T>& a, const Vector<N, T>& b) {
	for (size_t i = 0; i < N; i++) a[i] /= b[i];
	return a;
}

// a *= x
template <size_t N, typename T>
inline Vector<N, T>& operator *= (Vector<N, T>& a, T x) {
	for (size_t i = 0; i < N; i++) a[i] *= x;
	return a;
}

// a /= x
template <size_t N, typename T>
inline Vector<N, T>& operator /= (Vector<N, T>& a, T x) {
	for (size_t i = 0; i < N; i++) a[i] /= x;
	return a;
}

template<size_t N1, size_t N2, typename T>
inline Vector<N1, T> vector_convert(const Vector<N2, T>& a, T fill = 1) {
	Vector<N1, T> b;
	for (size_t i = 0; i < N1; i++)
		b[i] = (i < N2) ? a[i] : fill;
	return b;
}

// = |a| ^ 2
template<size_t N, typename T>
inline T vector_length_square(const Vector<N, T>& a) {
	T sum = 0;
	for (size_t i = 0; i < N; i++) sum += a[i] * a[i];
	return sum;
}

// = |a|
template<size_t N, typename T>
inline T vector_length(const Vector<N, T>& a) {
	return sqrt(vector_length_square(a));
}

// = |a|,specialize for float
template<size_t N>
inline float vector_length(const Vector<N, float>& a) {
	return sqrtf(vector_length_square(a));
}

// = a / |a|
template<size_t N, typename T>
inline Vector<N, T> vector_normalize(const Vector<N, T>& a) {
	return a / vector_length(a);
}

//transfer to different dimension
template<size_t N1, size_t N2, typename T>
inline Vector<N1, T> vector_convert(const Vector<N2, T>& a, T fill = 1) {
	Vector<N1, T> b;
	for (size_t i = 0; i < N1; i++)
		b[i] = (i < N2) ? a[i] : fill;
	return b;
}

// = |a| ^ 2
template<size_t N, typename T>
inline T vector_length_square(const Vector<N, T>& a) {
	T sum = 0;
	for (size_t i = 0; i < N; i++) sum += a[i] * a[i];
	return sum;
}

// = |a|
template<size_t N, typename T>
inline T vector_length(const Vector<N, T>& a) {
	return sqrt(vector_length_square(a));
}

// = |a|, specialize for float
template<size_t N>
inline float vector_length(const Vector<N, float>& a) {
	return sqrtf(vector_length_square(a));
}

// = a / |a|
template<size_t N, typename T>
inline Vector<N, T> vector_normalize(const Vector<N, T>& a) {
	return a / vector_length(a);
}

// dot product
template<size_t N, typename T>
inline T vector_dot(const Vector<N, T>& a, const Vector<N, T>& b) {
	T sum = 0;
	for (size_t i = 0; i < N; i++) sum += a[i] * b[i];
	return sum;
}

// = a + (b - a) * t
template<size_t N, typename T>
inline Vector<N, T> vector_lerp(const Vector<N, T>& a, const Vector<N, T>& b, float t) {
	return a + (b - a) * t;
}

// 2D cross product
template<typename T>
inline T vector_cross(const Vector<2, T>& a, const Vector<2, T>& b) {
	return a.x * b.y - a.y * b.x;
}

// 3D cross product
template<typename T>
inline Vector<3, T> vector_cross(const Vector<3, T>& a, const Vector<3, T>& b) {
	return Vector<3, T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

// 4D cross product, leave the last dimension
template<typename T>
inline Vector<4, T> vector_cross(const Vector<4, T>& a, const Vector<4, T>& b) {
	return Vector<4, T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, a.w);
}

//set element to maxium
template<size_t N, typename T>
inline Vector<N, T> vector_max(const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = (a[i] > b[i]) ? a[i] : b[i];
	return c;
}

// set elements to min
template<size_t N, typename T>
inline Vector<N, T> vector_min(const Vector<N, T>& a, const Vector<N, T>& b) {
	Vector<N, T> c;
	for (size_t i = 0; i < N; i++) c[i] = (a[i] < b[i]) ? a[i] : b[i];
	return c;
}

// set the value inside max and min
template<size_t N, typename T>
inline Vector<N, T> vector_between(const Vector<N, T>& minx, const Vector<N, T>& maxx, const Vector<N, T>& x) {
	return vector_min(vector_max(minx, x), maxx);
}

// detect whether the two vector is close enough
template<size_t N, typename T>
inline bool vector_near(const Vector<N, T>& a, const Vector<N, T>& b, T dist) {
	return (vector_length_square(a - b) <= dist);
}

// detect whether the two vector is same under single
template<size_t N>
inline bool vector_near_equal(const Vector<N, float>& a, const Vector<N, float>& b, float e = 0.0001) {
	return vector_near(a, b, e);
}

// detect whether the two vector is same under double
template<size_t N>
inline bool vector_near_equal(const Vector<N, double>& a, const Vector<N, double>& b, double e = 0.0000001) {
	return vector_near(a, b, e);
}

// set the value inside scope, default inside[0,1]
template<size_t N, typename T>
inline Vector<N, T> vector_clamp(const Vector<N, T>& a, T minx = 0, T maxx = 1) {
	Vector<N, T> b;
	for (size_t i = 0; i < N; i++) {
		T x = (a[i] < minx) ? minx : a[i];
		b[i] = (x > maxx) ? maxx : x;
	}
	return b;
}

// print the value to os
template<size_t N, typename T>
inline std::ostream& operator << (std::ostream& os, const Vector<N, T>& a) {
	os << "[";
	for (size_t i = 0; i < N; i++) {
		os << a[i];
		if (i < N - 1) os << ", ";
	}
	os << "]";
	return os;
}

// print the value
template<size_t N, typename T>
inline std::string vector_repr(const Vector<N, T>& a) {
	std::stringstream ss;
	ss << a;
	return ss.str();
}

template<size_t ROW, size_t COL, typename T> struct Matrix {
	T m[ROW][COL];

	inline Matrix() {}

	inline Matrix(const Matrix<ROW, COL, T>& src) {
		for (size_t r = 0; r < ROW; r++) {
			for (size_t c = 0; c < COL; c++)
				m[r][c] = src.m[r][c];
		}
	}

	//use vector to initialize matrix
	inline Matrix(const std::initializer_list<Vector<COL, T>> &u) {
		auto it = u.begin();
		for (size_t i = 0; i < ROW; i++) SetRow(i, *it++);
	}

	//get element
	inline const T* operator [] (size_t row) const { assert(row < ROW); return m[row]; }
	inline T* operator [] (size_t row) { assert(row < ROW); return m[row]; }

	// get row
	inline Vector<COL, T> Row(size_t row) const {
		assert(row < ROW);
		Vector<COL, T> a;
		for (size_t i = 0; i < COL; i++) a[i] = m[row][i];
		return a;
	}

	// get col
	inline Vector<ROW, T> Col(size_t col) const {
		assert(col < COL);
		Vector<ROW, T> a;
		for (size_t i = 0; i < ROW; i++) a[i] = m[i][col];
		return a;
	}

	// set row
	inline void SetRow(size_t row, const Vector<COL, T>& a) {
		assert(row < ROW);
		for (size_t i = 0; i < COL; i++) m[row][i] = a[i];
	}

	// set col
	inline void SetCol(size_t col, const Vector<ROW, T>& a) {
		assert(col < COL);
		for (size_t i = 0; i < ROW; i++) m[i][col] = a[i];
	}

	// get minor matrix from delete a col && row
	inline Matrix<ROW - 1, COL - 1, T> GetMinor(size_t row, size_t col) const {
		Matrix<ROW - 1, COL - 1, T> ret;
		for (size_t r = 0; r < ROW - 1; r++) {
			for (size_t c = 0; c < COL - 1; c++) {
				ret.m[r][c] = m[r < row ? r : r + 1][c < col ? c : c + 1];
			}
		}
		return ret;
	}

	// get tranpose
	inline Matrix<COL, ROW, T> Transpose() const {
		Matrix<COL, ROW, T> ret;
		for (size_t r = 0; r < ROW; r++) {
			for (size_t c = 0; c < COL; c++)
				ret.m[c][r] = m[r][c];
		}
		return ret;
	}

	// get 0 matrix
	inline static Matrix<ROW, COL, T> GetZero() {
		Matrix<ROW, COL, T> ret;
		for (size_t r = 0; r < ROW; r++) {
			for (size_t c = 0; c < COL; c++)
				ret.m[r][c] = 0;
		}
		return ret;
	}

	// get identity
	inline static Matrix<ROW, COL, T> GetIdentity() {
		Matrix<ROW, COL, T> ret;
		for (size_t r = 0; r < ROW; r++) {
			for (size_t c = 0; c < COL; c++)
				ret.m[r][c] = (r == c) ? 1 : 0;
		}
		return ret;
	}

	//
};

//operator overload for matrix
template<size_t ROW, size_t COL, typename T>
inline bool operator == (const Matrix<ROW, COL, T>& a, const Matrix<ROW, COL, T>& b) {
	for (size_t r = 0; r < ROW; r++) {
		for (size_t c = 0; c < COL; c++) {
			if (a.m[r][c] != b.m[r][c]) return false;
		}
	}
	return true;
}

template<size_t ROW, size_t COL, typename T>
inline bool operator != (const Matrix<ROW, COL, T>& a, const Matrix<ROW, COL, T>& b) {
	return !(a == b);
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator + (const Matrix<ROW, COL, T>& src) {
	return src;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator - (const Matrix<ROW, COL, T>& src) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++)
			out.m[j][i] = -src.m[j][i];
	}
	return out;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator + (const Matrix<ROW, COL, T>& a, const Matrix<ROW, COL, T>& b) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++)
			out.m[j][i] = a.m[j][i] + b.m[j][i];
	}
	return out;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator - (const Matrix<ROW, COL, T>& a, const Matrix<ROW, COL, T>& b) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++)
			out.m[j][i] = a.m[j][i] - b.m[j][i];
	}
	return out;
}

template<size_t ROW, size_t COL, size_t NEWCOL, typename T>
inline Matrix<ROW, NEWCOL, T> operator * (const Matrix<ROW, COL, T>& a, const Matrix<COL, NEWCOL, T>& b) {
	Matrix<ROW, NEWCOL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < NEWCOL; i++) {
			out.m[j][i] = vector_dot(a.Row(j), b.Col(i));
		}
	}
	return out;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator * (const Matrix<ROW, COL, T>& a, T x) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++) {
			out.m[j][i] = a.m[j][i] * x;
		}
	}
	return out;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator / (const Matrix<ROW, COL, T>& a, T x) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++) {
			out.m[j][i] = a.m[j][i] / x;
		}
	}
	return out;
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator * (T x, const Matrix<ROW, COL, T>& a) {
	return (a * x);
}

template<size_t ROW, size_t COL, typename T>
inline Matrix<ROW, COL, T> operator / (T x, const Matrix<ROW, COL, T>& a) {
	Matrix<ROW, COL, T> out;
	for (size_t j = 0; j < ROW; j++) {
		for (size_t i = 0; i < COL; i++) {
			out.m[j][i] = x / a.m[j][i];
		}
	}
	return out;
}

//speicial case, return type is vector
template<size_t ROW, size_t COL, typename T>
inline Vector<COL, T> operator * (const Vector<ROW, T>& a, const Matrix<ROW, COL, T>& m) {
	Vector<COL, T> b;
	for (size_t i = 0; i < COL; i++)
		b[i] = vector_dot(a, m.Col(i));
	return b;
}

template<size_t ROW, size_t COL, typename T>
inline Vector<ROW, T> operator * (const Matrix<ROW, COL, T>& m, const Vector<COL, T>& a) {
	Vector<ROW, T> b;
	for (size_t i = 0; i < ROW; i++)
		b[i] = vector_dot(a, m.Row(i));
	return b;
}

// calculate determinant 1D
template<typename T>
inline T matrix_det(const Matrix<1, 1, T> &m) {
	return m[0][0];
}

// calculate determinant 2D
template<typename T>
inline T matrix_det(const Matrix<2, 2, T> &m) {
	return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

// calculate determinant
template<size_t N, typename T>
inline T matrix_det(const Matrix<N, N, T> &m) {
	T sum = 0;
	for (size_t i = 0; i < N; i++) sum += m[0][i] * matrix_cofactor(m, 0, i);
	return sum;
}

// cofactor 1D
template<typename T>
inline T matrix_cofactor(const Matrix<1, 1, T> &m, size_t row, size_t col) {
	return 0;
}

// cofactor
template<size_t N, typename T>
inline T matrix_cofactor(const Matrix<N, N, T> &m, size_t row, size_t col) {
	return matrix_det(m.GetMinor(row, col)) * (((row + col) % 2) ? -1 : 1);
}

// get adjugate matrix
template<size_t N, typename T>
inline Matrix<N, N, T> matrix_adjoint(const Matrix<N, N, T> &m) {
	Matrix<N, N, T> ret;
	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) 
			ret[j][i] = matrix_cofactor(m, i, j);
	}
	return ret;
}

// inverse matrix
template<size_t N, typename T>
inline Matrix<N, N, T> matrix_invert(const Matrix<N, N, T> &m) {
	Matrix<N, N, T> ret = matrix_adjoint(m);
	T det = vector_dot(m.Row(0), ret.Col(0));
	return ret / det;
}

// add to os 
template<size_t ROW, size_t COL, typename T>
inline std::ostream& operator << (std::ostream& os, const Matrix<ROW, COL, T>& m) {
	for (size_t r = 0; r < ROW; r++) {
		Vector<COL, T> row = m.Row(r);
		os << row << std::endl;
	}
	return os;
}

//tool function
template<typename T> inline T Abs(T x) { 
	return (x < 0) ? (-x) : x; 
}
template<typename T> inline T Max(T x, T y) { 
	return (x < y) ? y : x; 
}
template<typename T> inline T Min(T x, T y) { 
	return (x > y) ? y : x; 
}

template<typename T> inline bool NearEqual(T x, T y, T error) {
	return (Abs(x - y) < error);
}

template<typename T> inline T Between(T xmin, T xmax, T x) {
	return Min(Max(xmin, x), xmax);
}

// set varaible to [0,1]
template<typename T> inline T Saturate(T x) {
	return Between<T>(0, 1, x);
}

//define type name
typedef Matrix<4, 4, float> Mat4x4f;
typedef Matrix<3, 3, float> Mat3x3f;
typedef Matrix<4, 3, float> Mat4x3f;
typedef Matrix<3, 4, float> Mat3x4f;
typedef Vector<2, float>  Vec2f;
typedef Vector<2, double> Vec2d;
typedef Vector<2, int>    Vec2i;
typedef Vector<3, float>  Vec3f;
typedef Vector<3, double> Vec3d;
typedef Vector<3, int>    Vec3i;
typedef Vector<4, float>  Vec4f;
typedef Vector<4, double> Vec4d;
typedef Vector<4, int>    Vec4i;

// vector to color, use 4*8 Bytes to store it
inline static uint32_t vector_to_color(const Vec4f& color) {
	uint32_t r = (uint32_t)Between(0, 255, (int)(color.r * 255.0f));
	uint32_t g = (uint32_t)Between(0, 255, (int)(color.g * 255.0f));
	uint32_t b = (uint32_t)Between(0, 255, (int)(color.b * 255.0f));
	uint32_t a = (uint32_t)Between(0, 255, (int)(color.a * 255.0f));
	return (r << 16) | (g << 8) | b | (a << 24);
}

// vector to color without alpha
inline static uint32_t vector_to_color(const Vec3f& color) {
	return vector_to_color(color.xyz1());
}

// color to vector
inline static Vec4f vector_from_color(uint32_t rgba) {
	Vec4f out;
	out.r = ((rgba >> 16) & 0xff) / 255.0f;
	out.g = ((rgba >> 8) & 0xff) / 255.0f;
	out.b = ((rgba >> 0) & 0xff) / 255.0f;
	out.a = ((rgba >> 24) & 0xff) / 255.0f;
	return out;
}

// get zero matrix
inline static Mat4x4f matrix_set_zero() {
	Mat4x4f m;
	m.m[0][0] = m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][1] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][2] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = m.m[3][3] = 0.0f;
	return m;
}

// get identity matrix
inline static Mat4x4f matrix_set_identity() {
	Mat4x4f m;
	m.m[0][0] = m.m[1][1] = m.m[2][2] = m.m[3][3] = 1.0f;
	m.m[0][1] = m.m[0][2] = m.m[0][3] = 0.0f;
	m.m[1][0] = m.m[1][2] = m.m[1][3] = 0.0f;
	m.m[2][0] = m.m[2][1] = m.m[2][3] = 0.0f;
	m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
	return m;
}

// matrix for translate
inline static Mat4x4f matrix_set_translate(float x, float y, float z) {
	Mat4x4f m = matrix_set_identity();
	m.m[3][0] = x;
	m.m[3][1] = y;
	m.m[3][2] = z;
	return m;
}

// matrix for scale
inline static Mat4x4f matrix_set_scale(float x, float y, float z) {
	Mat4x4f m = matrix_set_identity();
	m.m[0][0] = x;
	m.m[1][1] = y;
	m.m[2][2] = z;
	return m;
}

// get rotation matrix
inline static Mat4x4f matrix_set_rotate(float x, float y, float z, float theta) {
	float qsin = (float)sin(theta * 0.5f);
	float qcos = (float)cos(theta * 0.5f);
	float w = qcos;
	Vec3f vec = vector_normalize(Vec3f(x, y, z));
	x = vec.x * qsin;
	y = vec.y * qsin;
	z = vec.z * qsin;
	Mat4x4f m;
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

// camera transform matrix：eye/camera position，at/position look at，up
inline static Mat4x4f matrix_set_lookat(const Vec3f& eye, const Vec3f& at, const Vec3f& up) {
	Vec3f zaxis = vector_normalize(at - eye);//vector point to position we look at
	Vec3f xaxis = vector_normalize(vector_cross(up, zaxis));//get a normal of the plane form by z && up
	Vec3f yaxis = vector_cross(zaxis, xaxis);//use  z and x to create a new y
	Mat4x4f m;
	m.SetCol(0, Vec4f(xaxis.x, xaxis.y, xaxis.z, -vector_dot(eye, xaxis)));
	m.SetCol(1, Vec4f(yaxis.x, yaxis.y, yaxis.z, -vector_dot(eye, yaxis)));
	m.SetCol(2, Vec4f(zaxis.x, zaxis.y, zaxis.z, -vector_dot(eye, zaxis)));
	m.SetCol(3, Vec4f(0.0f, 0.0f, 0.0f, 1.0f));
	return m;
}


// D3D Matrix Perspective FovLH
inline static Mat4x4f matrix_set_perspective(float fovy, float aspect, float zn, float zf) {
	float fax = 1.0f / (float)tan(fovy * 0.5f);
	Mat4x4f m = matrix_set_zero();
	m.m[0][0] = (float)(fax / aspect);
	m.m[1][1] = (float)(fax);
	m.m[2][2] = zf / (zf - zn);
	m.m[3][2] = -zn * zf / (zf - zn);
	m.m[2][3] = 1;
	return m;
}

//Bitmap class,rgba_8888
class Bitmap
{
public:
	inline virtual ~Bitmap() { 
		if (_bits) delete[]_bits; 
		_bits = NULL; 
	}

	inline Bitmap(int width, int height) : _w(width), _h(height) {
		_pitch = width * 4;
		_bits = new uint8_t[_pitch * _h];
		Fill(0);
	}

	inline Bitmap(const Bitmap& src) : _w(src._w), _h(src._h), _pitch(src._pitch) {
		_bits = new uint8_t[_pitch * _h];
		memcpy(_bits, src._bits, _pitch * _h);
	}

	inline Bitmap(const char *filename) {
		Bitmap *tmp = LoadFile(filename);
		if (tmp == NULL) {
			std::string msg = "load failed: ";
			msg.append(filename);
			throw std::runtime_error(msg);
		}
		_w = tmp->_w; _h = tmp->_h; _pitch = tmp->_pitch; _bits = tmp->_bits;
		tmp->_bits = NULL;
		delete tmp;
	}

public:
	inline int GetW() const { return _w; }
	inline int GetH() const { return _h; }
	inline int GetPitch() const { return _pitch; }
	inline uint8_t *GetBits() { return _bits; }
	inline const uint8_t *GetBits() const { return _bits; }
	inline uint8_t *GetLine(int y) { return _bits + _pitch * y; }
	inline const uint8_t *GetLine(int y) const { return _bits + _pitch * y; }

	inline void Fill(uint32_t color) {
		for (int j = 0; j < _h; j++) {
			uint32_t *row = (uint32_t*)(_bits + j * _pitch);
			for (int i = 0; i < _w; i++, row++)
				memcpy(row, &color, sizeof(uint32_t));
		}
	}

	inline void SetPixel(int x, int y, uint32_t color) {
		if (x >= 0 && x < _w && y >= 0 && y < _h) {
			memcpy(_bits + y * _pitch + x * 4, &color, sizeof(uint32_t));
		}
	}

	inline uint32_t GetPixel(int x, int y) const {
		uint32_t color = 0;
		if (x >= 0 && x < _w && y >= 0 && y < _h) {
			memcpy(&color, _bits + y * _pitch + x * 4, sizeof(uint32_t));
		}
		return color;
	}

	inline void DrawLine(int x1, int y1, int x2, int y2, uint32_t color) {
		int x, y;
		if (x1 == x2 && y1 == y2) {
			SetPixel(x1, y1, color);
			return;
		}
		else if (x1 == x2) {
			int inc = (y1 <= y2) ? 1 : -1;
			for (y = y1; y != y2; y += inc) SetPixel(x1, y, color);
			SetPixel(x2, y2, color);
		}
		else if (y1 == y2) {
			int inc = (x1 <= x2) ? 1 : -1;
			for (x = x1; x != x2; x += inc) SetPixel(x, y1, color);
			SetPixel(x2, y2, color);
		}
		else {
			int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
			int dy = (y1 < y2) ? y2 - y1 : y1 - y2;
			int rem = 0;
			if (dx >= dy) {
				if (x2 < x1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;//swap
				for (x = x1, y = y1; x <= x2; x++) {
					//draw the col
					SetPixel(x, y, color);
					rem += dy;
					//draw the line
					if (rem >= dx) { 
						rem -= dx; 
						y += (y2 >= y1) ? 1 : -1; 
						SetPixel(x, y, color); 
					}
				}
				//set destination color
				SetPixel(x2, y2, color);
			}
			else {
				if (y2 < y1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
				for (x = x1, y = y1; y <= y2; y++) {
					SetPixel(x, y, color);
					rem += dx;
					if (rem >= dy) { 
						rem -= dy; 
						x += (x2 >= x1) ? 1 : -1; 
						SetPixel(x, y, color); 
					}
				}
				SetPixel(x2, y2, color);
			}
		}
	}

	struct BITMAPINFOHEADER { // bmih,copy from doc
		uint32_t	biSize;
		uint32_t	biWidth;
		int32_t		biHeight;
		uint16_t	biPlanes;
		uint16_t	biBitCount;
		uint32_t	biCompression;
		uint32_t	biSizeImage;
		uint32_t	biXPelsPerMeter;
		uint32_t	biYPelsPerMeter;
		uint32_t	biClrUsed;
		uint32_t	biClrImportant;
	};

	// read image
	inline static Bitmap* LoadFile(const char *filename) {
		FILE *fp = fopen(filename, "rb");
		if (fp == NULL) return NULL;
		BITMAPINFOHEADER info;
		uint8_t header[14];
		int hr = (int)fread(header, 1, 14, fp);
		if (hr != 14) { fclose(fp); return NULL; }
		if (header[0] != 0x42 || header[1] != 0x4d) { fclose(fp); return NULL; }
		hr = (int)fread(&info, 1, sizeof(info), fp);
		if (hr != 40) { fclose(fp); return NULL; }
		if (info.biBitCount != 24 && info.biBitCount != 32) { fclose(fp); return NULL; }
		Bitmap *bmp = new Bitmap(info.biWidth, info.biHeight);
		uint32_t offset;
		memcpy(&offset, header + 10, sizeof(uint32_t));
		fseek(fp, offset, SEEK_SET);
		uint32_t pixelsize = (info.biBitCount + 7) / 8;
		uint32_t pitch = (pixelsize * info.biWidth + 3) & (~3);
		for (int y = 0; y < (int)info.biHeight; y++) {
			uint8_t *line = bmp->GetLine(info.biHeight - 1 - y);
			for (int x = 0; x < (int)info.biWidth; x++, line += 4) {
				line[3] = 255;
				fread(line, pixelsize, 1, fp);
			}
			fseek(fp, pitch - info.biWidth * pixelsize, SEEK_CUR);
		}
		fclose(fp);
		return bmp;
	}

	// save file
	inline bool SaveFile(const char *filename, bool withAlpha = false) const {
		FILE *fp = fopen(filename, "wb");
		if (fp == NULL) return false;
		BITMAPINFOHEADER info;
		uint32_t pixelsize = (withAlpha) ? 4 : 3;
		uint32_t pitch = (GetW() * pixelsize + 3) & (~3);
		info.biSizeImage = pitch * GetH();
		uint32_t bfSize = 54 + info.biSizeImage;
		uint32_t zero = 0, offset = 54;
		fputc(0x42, fp);
		fputc(0x4d, fp);
		fwrite(&bfSize, 4, 1, fp);
		fwrite(&zero, 4, 1, fp);
		fwrite(&offset, 4, 1, fp);
		info.biSize = 40;
		info.biWidth = GetW();
		info.biHeight = GetH();
		info.biPlanes = 1;
		info.biBitCount = (withAlpha) ? 32 : 24;
		info.biCompression = 0;
		info.biXPelsPerMeter = 0xb12;
		info.biYPelsPerMeter = 0xb12;
		info.biClrUsed = 0;
		info.biClrImportant = 0;
		fwrite(&info, sizeof(info), 1, fp);
		for (int y = 0; y < GetH(); y++) {
			const uint8_t *line = GetLine(info.biHeight - 1 - y);
			uint32_t padding = pitch - GetW() * pixelsize;
			for (int x = 0; x < GetW(); x++, line += 4) {
				fwrite(line, pixelsize, 1, fp);
			}
			for (int i = 0; i < (int)padding; i++) fputc(0, fp);
		}
		fclose(fp);
		return true;
	}

	// 双线性插值
	inline uint32_t SampleBilinear(float x, float y) const {
		int32_t fx = (int32_t)(x * 0x10000);
		int32_t fy = (int32_t)(y * 0x10000);
		int32_t x1 = Between(0, _w - 1, fx >> 16);
		int32_t y1 = Between(0, _h - 1, fy >> 16);
		int32_t x2 = Between(0, _w - 1, x1 + 1);
		int32_t y2 = Between(0, _h - 1, y1 + 1);
		int32_t dx = (fx >> 8) & 0xff;
		int32_t dy = (fy >> 8) & 0xff;
		if (_w <= 0 || _h <= 0) return 0;
		uint32_t c00 = GetPixel(x1, y1);
		uint32_t c01 = GetPixel(x2, y1);
		uint32_t c10 = GetPixel(x1, y2);
		uint32_t c11 = GetPixel(x2, y2);
		return BilinearInterp(c00, c01, c10, c11, dx, dy);
	}

	// 纹理采样
	inline Vec4f Sample2D(float u, float v) const {
		uint32_t rgba = SampleBilinear(u * _w + 0.5f, v * _h + 0.5f);
		return vector_from_color(rgba);
	}

	// 纹理采样：直接传入 Vec2f
	inline Vec4f Sample2D(const Vec2f& uv) const {
		return Sample2D(uv.x, uv.y);
	}

	// 按照 Vec4f 画点
	inline void SetPixel(int x, int y, const Vec4f& color) {
		SetPixel(x, y, vector_to_color(color));
	}

	// 上下反转
	inline void FlipVertical() {
		uint8_t *buffer = new uint8_t[_pitch];
		for (int i = 0, j = _h - 1; i < j; i++, j--) {
			memcpy(buffer, GetLine(i), _pitch);
			memcpy(GetLine(i), GetLine(j), _pitch);
			memcpy(GetLine(j), buffer, _pitch);
		}
		delete[]buffer;
	}

	// 水平反转
	inline void FlipHorizontal() {
		for (int y = 0; y < _h; y++) {
			for (int i = 0, j = _w - 1; i < j; i++, j--) {
				uint32_t c1 = GetPixel(i, y);
				uint32_t c2 = GetPixel(j, y);
				SetPixel(i, y, c2);
				SetPixel(j, y, c1);
			}
		}
	}

protected:

	// 双线性插值计算：给出四个点的颜色，以及坐标偏移，计算结果
	inline static uint32_t BilinearInterp(uint32_t tl, uint32_t tr,
		uint32_t bl, uint32_t br, int32_t distx, int32_t disty) {
		uint32_t f, r;
		int32_t distxy = distx * disty;
		int32_t distxiy = (distx << 8) - distxy;  /* distx * (256 - disty) */
		int32_t distixy = (disty << 8) - distxy;  /* disty * (256 - distx) */
		int32_t distixiy = 256 * 256 - (disty << 8) - (distx << 8) + distxy;
		r = (tl & 0x000000ff) * distixiy + (tr & 0x000000ff) * distxiy
			+ (bl & 0x000000ff) * distixy + (br & 0x000000ff) * distxy;
		f = (tl & 0x0000ff00) * distixiy + (tr & 0x0000ff00) * distxiy
			+ (bl & 0x0000ff00) * distixy + (br & 0x0000ff00) * distxy;
		r |= f & 0xff000000;
		tl >>= 16; tr >>= 16; bl >>= 16; br >>= 16; r >>= 16;
		f = (tl & 0x000000ff) * distixiy + (tr & 0x000000ff) * distxiy
			+ (bl & 0x000000ff) * distixy + (br & 0x000000ff) * distxy;
		r |= f & 0x00ff0000;
		f = (tl & 0x0000ff00) * distixiy + (tr & 0x0000ff00) * distxiy
			+ (bl & 0x0000ff00) * distixy + (br & 0x0000ff00) * distxy;
		r |= f & 0xff000000;
		return r;
	}

protected:
	int32_t _w;//width
	int32_t _h;//height
	int32_t _pitch;//how many Byts we need for each line
	uint8_t *_bits;
};

#endif // !_CPP_RENDER_
