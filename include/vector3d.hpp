#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <cmath>
#include <iostream>

using namespace std;

/* 3-D vector class */
class vector3d {
    friend istream & operator >> (istream & is, vector3d & v) {
        is >> v._x >> v._y >> v._z;
        return is;
    }
    friend ostream & operator << (ostream & os, const vector3d & v) {
        os << v._x << " " << v._y<< " " << v._z;
        return os;
    }
    friend vector3d operator + (const vector3d & v1, const vector3d & v2) {
        vector3d res(v1);
        res+=v2;
        return res;
    }
    friend vector3d operator - (const vector3d & v1, const vector3d & v2) {
        vector3d res(v1);
        res-=v2;
        return res;
    }
    friend vector3d operator * (const double & c, const vector3d & p) {
        vector3d res=p;
        res*=c;
        return res;
    }
    friend vector3d operator * (const vector3d & p, const double & c) {
        return c*p;
    }
    friend double Abs3d(const vector3d & v) {
        return abs(v._x) + abs(v._y) + abs(v._z);
    }
    friend vector3d Abs3d(const vector3d & v1, const vector3d & v2) {
        vector3d res;
        res._x = abs(v1._x - v2._x); res._y = abs(v1._y - v2._y); res._z = abs(v1._z - v2._z);
        return res;
    }
    friend double MagSq3d(const vector3d & v) {
        return v._x*v._x + v._y*v._y + v._z*v._z;
    }
    friend double Norm3d(const vector3d & v) {
        return sqrt(v._x*v._x + v._y*v._y + v._z*v._z);
    }
    friend double Distance3d(const vector3d & v1, const vector3d & v2) {
        return sqrt( pow(v1._x-v2._x, 2) + pow(v1._y-v2._y, 2) + pow(v1._z-v2._z, 2) );
    }
    friend double ScalProd3d(const vector3d & v1, const vector3d & v2) {
        return v1._x*v2._x + v1._y*v2._y + v1._z*v2._z;
    }
    friend vector3d CrossProd3d(const vector3d & v1, const vector3d & v2) {
        vector3d res;
        res._x = v1._y * v2._z - v1._z * v2._y;
        res._y = v1._z * v2._x - v1._x * v2._z;
        res._z = v1._x * v2._y - v1._y * v2._x;
        return res;
    }
    friend double InterAngle(const vector3d & v1, const vector3d & v2) { /* intersect angle */
        return acos( ScalProd3d(v1,v2) / ( Norm3d(v1) * Norm3d(v2) ) );
    }

    public:
    explicit vector3d(double x=0,double y=0,double z=0): _x(x), _y(y), _z(z){};

    double & x() {return _x;}
    double x() const {return _x;}
    double & y() {return _y;}
    double y() const {return _y;}
    double & z() {return _z;}
    double z() const {return _z;}

    const vector3d & operator += (const vector3d & p){
        _x+=p._x; _y+=p._y; _z+=p._z;
        return *this;
    }
    const vector3d & operator -= (const vector3d & p){
        _x-=p._x; _y-=p._y; _z-=p._z;
        return *this;
    }
    const vector3d & operator *= (const double & c){
        _x*=c; _y*=c; _z*=c;
        return *this;
    }
    void copy(const vector3d & v) {
        _x = v._x;
        _y = v._y;
        _z = v._z;
    }
    void Minus(const vector3d & v1, const vector3d & v2) {
        _x = v1._x - v2._x;
        _y = v1._y - v2._y;
        _z = v1._z - v2._z;
    }
    void Add(const vector3d & v1, const vector3d & v2) {
        _x = v1._x + v2._x;
        _y = v1._y + v2._y;
        _z = v1._z + v2._z;
    }
    void AddScale(const vector3d & v1, const vector3d & v2, const double & c) {
        _x = v1._x + v2._x * c;
        _y = v1._y + v2._y * c;
        _z = v1._z + v2._z * c;
    }
    void MinusScale(const vector3d & v1, const vector3d & v2, const double & c) {
        _x = v1._x - v2._x * c;
        _y = v1._y - v2._y * c;
        _z = v1._z - v2._z * c;
    }
    void SetVec(const double & c1, const double & c2, const double & c3){
        _x = c1; _y = c2; _z = c3;
    }
    void Setx(const double & c) {
        _x = c;
    }
    void Sety(const double & c) {
        _y = c;
    }
    void Setz(const double & c) {
        _z = c;
    }
    int octant() {
        if (_x >= 0 && _y >= 0 && _z >= 0)
            return 0;
        else if (_x < 0 && _y >= 0 && _z >= 0)
            return 1;
        else if (_x < 0 && _y < 0 && _z >= 0)
            return 2;
        else if (_x >= 0 && _y < 0 && _z >= 0)
            return 3;
        else if (_x >= 0 && _y < 0 && _z < 0)
            return 4;
        else if (_x < 0 && _y < 0 && _z < 0)
            return 5;
        else if (_x < 0 && _y >= 0 && _z < 0)
            return 6;
        else
            return 7;
    }
    private:
    double _x,_y,_z;
};

const vector3d null(0,0,0);
#endif
