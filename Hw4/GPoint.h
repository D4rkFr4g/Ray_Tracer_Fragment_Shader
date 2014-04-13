#pragma once

#include <math.h>
#include <GL/glew.h>

class GPoint
/*
PURPOSE: Used to encapsulate the properties and operations
   for 3 dimension points/vectors used to describe our scene
   
REMARK: Many of the operations for Points will be used very often
   during ray-tracing. Since they are short and we want them to be inlined
   we define them in the class defintion itself using a standard one-line
   format (not exactly like the coding guidelines)
*/
{
   private:
      GLdouble _x;
      GLdouble _y;
      GLdouble _z;

   public:
      GPoint();
      GPoint(GLdouble a, GLdouble b, GLdouble c);
      GPoint(const GLdouble pt[]);
      GPoint(const GPoint& p);
		~GPoint();
      
		GLdouble x();
      GLdouble y();
      GLdouble z();
      void set(GLdouble a, GLdouble b, GLdouble c);
      bool isZero();
      GLdouble length();
      void normalize();
      
      friend GPoint operator*(GLdouble scalar, const GPoint &other); //scalar products
      friend GPoint operator*(const GPoint &other, GLdouble scalar);         

      GPoint& operator*=(GLdouble scalar);
      GPoint& operator/=(GLdouble scalar);
      GPoint operator*(const GPoint& other); //cross product
      GLdouble operator&(const GPoint& other); //dot Product
      GPoint operator%(const GPoint& other);//Hadamard Product
      GPoint operator+(const GPoint &other) const;
      GPoint operator-(const GPoint &other) const;
      GPoint& operator+=(const GPoint &other);
      GPoint& operator-=(const GPoint &other);
      GPoint& operator=(const GPoint &other);
};
