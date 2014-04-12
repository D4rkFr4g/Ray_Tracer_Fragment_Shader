#include "GPoint.h"


GPoint::GPoint(void)
{
	_x = 0; _y = 0; _z =0;
}

GPoint::~GPoint(void)
{
}


GPoint::GPoint(GLdouble a, GLdouble b, GLdouble c)
   {_x = a; _y = b; _z =c;}
     
GPoint::GPoint(const GLdouble pt[])
   {_x = pt[0]; _y = pt[1]; _z = pt[2];}
                
GPoint::GPoint(const GPoint& p)
   {_x = p._x; _y = p._y; _z = p._z;}
         
      
GLdouble GPoint::x(){return _x;}
GLdouble GPoint::y(){return _y;}
GLdouble GPoint::z(){return _z;}
                  
void GPoint::set(GLdouble a, GLdouble b, GLdouble c)
   {_x = a; _y = b; _z =c;}
      
bool GPoint::isZero(){ return (_x == 0 && _y ==0 && _z==0);}
GLdouble GPoint::length(){return sqrt(_x*_x + _y*_y + _z*_z); }
void GPoint::normalize(){ GLdouble l = length(); _x /=l; _y/= l; _z /= l;}
      
//GPoint operator*(GLdouble scalar, const GPoint &other); //scalar products
//GPoint operator*(const GPoint &other, GLdouble scalar);         

GPoint& GPoint::operator*=(GLdouble scalar)
   {_x *= scalar; _y *= scalar; _z *= scalar; return *this;}         

GPoint& GPoint::operator/=(GLdouble scalar)
   {_x /= scalar; _y /= scalar; _z /= scalar; return *this;}         
      
GPoint GPoint::operator*(const GPoint& other) //cross product
   {return GPoint(_y*other._z - other._y*_z, _z*other._x - _x*other._z, _x*other._y - _y*other._x); }

GLdouble GPoint::operator&(const GPoint& other)//dot Product
   {return _x * other._x + _y * other._y + _z * other._z; }

GPoint GPoint::operator%(const GPoint& other)//Hadamard Product
   {return GPoint(_x * other._x, _y * other._y, _z * other._z); }
         
GPoint GPoint::operator+(const GPoint &other) const
   {return GPoint(_x + other._x, _y + other._y, _z + other._z); }
            
GPoint GPoint::operator-(const GPoint &other) const
   {return GPoint(_x - other._x, _y - other._y, _z - other._z); }
          
GPoint& GPoint::operator+=(const GPoint &other)
   {_x += other._x; _y += other._y; _z += other._z; return *this;}
            
GPoint& GPoint::operator-=(const GPoint &other)
   {_x -= other._x; _y -= other._y; _z -= other._z; return *this;}

GPoint& GPoint::operator=(const GPoint &other)
   {_x = other._x; _y = other._y; _z = other._z; return *this;}