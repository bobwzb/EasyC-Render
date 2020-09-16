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
#endif // !_CPP_RENDER_
