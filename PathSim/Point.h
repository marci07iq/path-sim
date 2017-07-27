#pragma once

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <ctime>
#include <random>


using namespace std;

template<typename T> struct vec3;
typedef vec3<double>         dVec3;
typedef vec3<int>            iVec3;

template<typename T>
class vec3 {
public:
  enum Exception { WRONG_INDEX };

  T x;
  T y;
  T z;

  vec3() : x(0), y(0), z(0) {}
  vec3(T value) : x(value), y(value), z(value) {};
  vec3(T ax, T ay, T az) : x(ax), y(ay), z(az) {};
  /*vec3(proxy p) : x(*p.px), y(*p.py), z(*p.pz) {};*/
  template<typename U> vec3(vec3<U> c)
    : x((T)c.x), y((T)c.y), z((T)c.z) {};

  vec3<T> & operator+=(vec3<T> rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
  vec3<T> & operator-=(vec3<T> rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *this;
  }
  vec3<T> & operator*=(vec3<T> rhs) {
    x *= rhs.x; y *= rhs.y; z *= rhs.z;
    return *this;
  }
  vec3<T> & operator*=(T rhs) {
    x *= rhs; y *= rhs; z *= rhs;
    return *this;
  }
  vec3<T> & operator/=(vec3<T> rhs) {
    x /= rhs.x; y /= rhs.y; z /= rhs.z;
    return *this;
  }
  vec3<T> & operator%=(vec3<T> rhs) {
    x %= rhs.x; y %= rhs.y; z %= rhs.z;
    return *this;
  }
  vec3<T> & operator&=(vec3<T> rhs) {
    x &= rhs.x; y &= rhs.y; z &= rhs.z;
    return *this;
  }
  vec3<T> & operator|=(vec3<T> rhs) {
    x |= rhs.x; y |= rhs.y; z |= rhs.z;
    return *this;
  }
  vec3<T> & operator^=(vec3<T> rhs) {
    x ^= rhs.x; y ^= rhs.y; z ^= rhs.z;
    return *this;
  }

  vec3<T> operator+ (vec3<T> rhs) const { return vec3<T>(*this) += rhs; }
  vec3<T> operator- (vec3<T> rhs) const { return vec3<T>(*this) -= rhs; }
  vec3<T> operator* (vec3<T> rhs) const { return vec3<T>(*this) *= rhs; }
  vec3<T> operator* (T rhs) const { return vec3<T>(*this) *= rhs; }
  vec3<T> operator/ (vec3<T> rhs) const { return vec3<T>(*this) /= rhs; }
  vec3<T> operator% (vec3<T> rhs) const { return vec3<T>(*this) %= rhs; }
  vec3<T> operator& (vec3<T> rhs) const { return vec3<T>(*this) &= rhs; }
  vec3<T> operator| (vec3<T> rhs) const { return vec3<T>(*this) |= rhs; }
  vec3<T> operator^ (vec3<T> rhs) const { return vec3<T>(*this) ^= rhs; }
  /*vec3<T> operator<<(vec3<T> rhs) const { return vec3<T>(*this) <<= rhs; }
  vec3<T> operator >> (vec3<T> rhs) const { return vec3<T>(*this) >>= rhs; }*/

  vec3<T> operator+() const { return vec3<T>(+x, +y, +z); }
  vec3<T> operator-() const { return vec3<T>(-x, -y, -z); }
  vec3<T> operator~() const { return vec3<T>(~x, ~y, ~z); }
  vec3<T> operator!() const { return vec3<T>(!x, !y, !z); }

  vec3<T> & operator++() { ++x, ++y, ++z; return *this; }
  vec3<T> & operator--() { --x, --y, --z; return *this; }
  vec3<T> operator++(int) { vec3<T> cp(*this); ++x, ++y, ++z; return cp; }
  vec3<T> operator--(int) { vec3<T> cp(*this); --x, --y, --z; return cp; }

  T& operator[](int n)
  {
    switch (n)
    {
    case 0:
      return x;
      break;

    case 1:
      return y;
      break;

    case 2:
      return z;
      break;

    default:
      throw WRONG_INDEX;
    }
  }
  T sqrlen() const {
    return x*x + y*y + z*z;
  }
  double length() const {
    return sqrt(double(this->sqrlen()));
  }
  void norm() {
    *this /= length();
  }
};

template<typename T> std::ostream& operator<<(std::ostream& os, const vec3<T>& v)
{
  os << v.x << "\t" << v.y << "\t" << v.z;
  return os;
}
template<typename T> std::istream& operator>>(std::istream& is, vec3<T>& v)
{
  is >> v.x >> v.y >> v.z;
  return is;
}

template<typename T> inline T dot(vec3<T> lhs, vec3<T> rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
template<typename T> inline vec3<T> crs(vec3<T> lhs, vec3<T> rhs) {
  return vec3<T>(lhs.y*rhs.z - lhs.z*rhs.y,
    lhs.z*rhs.x - lhs.x*rhs.z,
    lhs.x*rhs.y - lhs.y*rhs.x);
}

template<typename T> inline double distance(vec3<T> lhs, vec3<T> rhs) {
  return (lhs-rhs).length();
}