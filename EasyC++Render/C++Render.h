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

// = a * b£¬not dot product, use for color multiple
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

#endif // !_CPP_RENDER_
