/****************************************************** 
* Copyright (c):   2014, All Rights Reserved. 
* Project:         CS 116B Homework #4
* File:            MySdlAppliation.cpp 
* Purpose:         Be able to implement a ray-tracing algorithm.
* Start date:      4/09/14 
* Programmer:      Zane Melcho
* 
****************************************************** 
*/
#define SDL_MAIN_HANDLED
#include "MySdlApplication.h"
//#include "GPoint.h"

//   CONSTANTS
enum { LIGHT, TETRAHEDRON, CUBE, SPHERE, CYLINDER, CONE };
static const bool G_GL2_COMPATIBLE = false;
static const unsigned char* KB_STATE = NULL;
//static const float G_FRUST_MIN_FOV = 60.0;  //A minimal of 60 degree field of view
static const int G_NUM_SHADERS = 1;

static const char * const G_SHADER_FILES[G_NUM_SHADERS][2] = 
{
   {"./Shaders/basic-gl3.vshader", "./Shaders/solid-gl3.fshader"},
};
static const char * const G_SHADER_FILES_GL2[G_NUM_SHADERS][2] = 
{
   {"./Shaders/basic-gl2.vshader", "./Shaders/solid-gl2.fshader"},
};

const GLdouble WHITE[3] = {1.0, 1.0, 1.0}; //RGB for white
const GLdouble BLACK[3] = {0.0, 0.0, 0.0}; //RGB for black
const GLdouble RED[3] = {1.0, 0.0, 0.0}; //RGB for RED

const GLdouble ATTENUATION_FACTOR = 100000; 
// used in our lighting equations to model how light attenuates with distance

const GLdouble CAMERA_POSITION[3] = {0, 100, 200}; //initial position of camera
const GLdouble LOOK_AT_VECTOR[3] = {0, 0, -160}; // where camera is looking at
const GLdouble UP_VECTOR[3] = {0, 1, 0}; // what direction is up for the camera

const GLdouble BOARD_POSITION[3] = {0, 0, -160}; // where in the g_scene the board is positioned   
const GLdouble BOARD_EDGE_SIZE = 320.0; // how wide the board is
const GLdouble BOARD_HALF_SIZE = BOARD_EDGE_SIZE/2; // useful to half this size ready for some calculations
const unsigned int NUM_SQUARES = 8; //number of squares wide our chess board is
const GLdouble SQUARE_EDGE_SIZE = BOARD_EDGE_SIZE/NUM_SQUARES; // pixels/per square

const unsigned int MAX_DEPTH = 5; // maximum depth our ray-tracing tree should go to

const GLdouble SMALL_NUMBER = .0001; // used rather than check with zero to avoid round-off problems 

const GLdouble SUPER_SAMPLE_NUMBER = 16; // how many random rays per pixel

struct ShaderState 
{
   GlProgram program;

   // Handles to uniform variables
   GLint h_uProjMatrix;
   GLint h_uModelViewMatrix;
   GLint h_uNormalMatrix;
   /*
   GLint h_uLight, h_uLight2;
   GLint h_uColor;
   GLint h_uTexUnit0;
   GLint h_uTexUnit1;
   GLint h_uTexUnit2;
   GLint h_uSamples;
   GLint h_uSampledx;
   GLint h_uSampledy;
   */

   // Handles to vertex attributes
   GLint h_aPosition;
   GLint h_aNormal;
   /*
   GLint h_aTangent;
   GLint h_aTexCoord0;
   GLint h_aTexCoord1;
   GLint h_aTexCoord2;
   */

   /*-----------------------------------------------*/
   ShaderState(const char* vsfn, const char* fsfn) 
   {
      /*	PURPOSE:		Constructor for ShaderState Object 
      RECEIVES:	vsfn - Vertex Shader Filename
      fsfn - Fragement Shader Filename
      RETURNS:		ShaderState object
      REMARKS:		
      */

      readAndCompileShader(program, vsfn, fsfn);

      const GLuint h = program; // short hand

      // Retrieve handles to uniform variables
      h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
      h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
      //h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");

      /*
      h_uLight = safe_glGetUniformLocation(h, "uLight");
      h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
      h_uColor = safe_glGetUniformLocation(h, "uColor");
      h_uTexUnit0 = safe_glGetUniformLocation(h, "uTexUnit0");
      h_uTexUnit1 = safe_glGetUniformLocation(h, "uTexUnit1");
      h_uTexUnit2 = safe_glGetUniformLocation(h, "uTexUnit2");
      h_uSamples = safe_glGetUniformLocation(h, "uSamples");
      h_uSampledx = safe_glGetUniformLocation(h, "uSampledx");
      h_uSampledy = safe_glGetUniformLocation(h, "uSampledy");
      */

      // Retrieve handles to vertex attributes
      h_aPosition = safe_glGetAttribLocation(h, "aPosition");
      //h_aNormal = safe_glGetAttribLocation(h, "aNormal");
      /*
      h_aTangent = safe_glGetAttribLocation(h, "aTangent");
      h_aTexCoord0 = safe_glGetAttribLocation(h, "aTexCoord0");
      h_aTexCoord1 = safe_glGetAttribLocation(h, "aTexCoord1");
      h_aTexCoord2 = safe_glGetAttribLocation(h, "aTexCoord2");
      */

      if (!G_GL2_COMPATIBLE)
         glBindFragDataLocation(h, 0, "fragColor");
      checkGlErrors();
   }
   /*-----------------------------------------------*/
};
/*-----------------------------------------------*/
// PROTOTYPES
class RayObject;

//   CLASS DEFINITIONS
/*-----------------------------------------------*/
class Point
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
   Point()
   {_x = 0; _y = 0; _z =0;}

   Point(GLdouble a, GLdouble b, GLdouble c)
   {_x = a; _y = b; _z =c;}

   Point(const GLdouble pt[])
   {_x = pt[0]; _y = pt[1]; _z = pt[2];}

   Point(const Point& p)
   {_x = p._x; _y = p._y; _z = p._z;}


   GLdouble x(){return _x;}
   GLdouble y(){return _y;}
   GLdouble z(){return _z;}

   void set(GLdouble a, GLdouble b, GLdouble c)
   {_x = a; _y = b; _z =c;}

   bool isZero(){ return (_x == 0 && _y ==0 && _z==0);}
   GLdouble length(){return sqrt(_x*_x + _y*_y + _z*_z); }
   void normalize(){ GLdouble l = length(); _x /=l; _y/= l; _z /= l;}

   friend Point operator*(GLdouble scalar, const Point &other); //scalar products
   friend Point operator*(const Point &other, GLdouble scalar);         

   Point& operator*=(GLdouble scalar)
   {_x *= scalar; _y *= scalar; _z *= scalar; return *this;}         

   Point& operator/=(GLdouble scalar)
   {_x /= scalar; _y /= scalar; _z /= scalar; return *this;}         

   Point operator*(const Point& other) //cross product
   {return Point(_y*other._z - other._y*_z, _z*other._x - _x*other._z, _x*other._y - _y*other._x); }

   GLdouble operator&(const Point& other)//dot Product
   {return _x * other._x + _y * other._y + _z * other._z; }

   Point operator%(const Point& other)//Hadamard Product
   {return Point(_x * other._x, _y * other._y, _z * other._z); }


   Point operator+(const Point &other) const
   {return Point(_x + other._x, _y + other._y, _z + other._z); }

   Point operator-(const Point &other) const
   {return Point(_x - other._x, _y - other._y, _z - other._z); }

   Point& operator+=(const Point &other)
   {_x += other._x; _y += other._y; _z += other._z; return *this;}

   Point& operator-=(const Point &other)
   {_x -= other._x; _y -= other._y; _z -= other._z; return *this;}

   Point& operator=(const Point &other)
   {_x = other._x; _y = other._y; _z = other._z; return *this;}


};
/*-----------------------------------------------*/
class Light
   /*
   PURPOSE: Use to encapsulate information about a 
   light (basically its color and position) in our g_scene
   to be ray-traced

   REMARK:         
   */
{
private:
   Point _color;
   Point _position;

public:

   Light(const Point& c, const Point& p){ _color = c; _position = p;}
   Point color(){return _color;}
   Point position(){return _position;}  
};
/*-----------------------------------------------*/
class Line
   /*
   PURPOSE: Used to encapsulate information about lines in our g_scene

   REMARK:    
   A ray that will be ray traced will be such a line
   We will adopt the convention that the _startPoint of such a line is
   where it is coming from and the the _endPoint is mainly used to
   specify the direction of the line
   */
{
private:
   Point _startPt;
   Point _endPt;

public:
   Line(){_startPt.set(0.0, 0.0, 0.0); _endPt.set(0.0, 0.0, 0.0);}
   Line(const Point& p1, const Point& p2) {_startPt = p1; _endPt = p2;}

   void set(const Point& p1, const Point& p2){_startPt = p1; _endPt = p2;}

   Point startPoint() const {return _startPt;}
   Point endPoint() const {return _endPt;}

   Point direction() const
   {
      Point p =  _endPt - _startPt;
      p.normalize();
      return p;
   }

   GLdouble length() const
   {
      Point p = _endPt - _startPt;
      return p.length();
   }
};
/*-----------------------------------------------*/
class Material
   /*
   PURPOSE: Used to hold information about how a Shape
   will react to various kinds of light in our lighting model

   REMARK:
   We use the Phong lighting model as the basis for calculating
   lighting when we do ray-tracing
   */
{
private:
   Point _ambient;
   Point _diffuse;
   Point _specular;
   Point _transparency;
   GLdouble _refraction;


public:
   Material()
   { _ambient.set(0.0, 0.0, 0.0); _diffuse =_ambient; _specular = _ambient; _transparency = _ambient;
   _refraction = 1;}
   Material(const Point& a, const Point& d, const Point& s, const Point& t, GLdouble r)
   {_ambient = a; _diffuse = d; _specular = s; _transparency = t; _refraction = r;}
   Material(const Material& m)
   {_ambient = m._ambient; _diffuse = m._diffuse; _specular = m._specular; 
   _transparency = m._transparency; _refraction = m._refraction;}


   Point ambient(){ return _ambient;}
   Point diffuse(){ return _diffuse;}
   Point specular(){ return _specular;}
   Point transparency(){ return _transparency;}
   GLdouble refraction(){ return _refraction;}

};
/*-----------------------------------------------*/
class Intersection
   /*
   PURPOSE: used to store information about how a ray intersects with a RayObject.

   REMARK:
   A flag is used to say if there was any intersection at all. If there is an intersection
   objects of this class store what was the material of the object that was intersected as well
   as its normal, the transmitted ray and the reflected ray. 

   */
{
private:
   bool _intersects;
   Point _point;
   Point _normal;

   Material _material;
   Line _reflectedRay;
   Line _transmittedRay;

public:
   Intersection(){}
   Intersection(bool intersects, const Point& p, const Point& n, const Material& m, const Line& 
      r, const Line& t)
   {_intersects = intersects; _point = p; _normal = n; _material = m; _reflectedRay = r; _transmittedRay =
   t;}


   bool intersects(){return _intersects;}

   Point point(){return _point;}
   Point normal(){return _normal;}

   Material material(){return _material;}

   Line reflectedRay(){return _reflectedRay;}
   Line transmittedRay(){return _transmittedRay;}

   void setIntersect(bool i){_intersects = i;}
   void setMaterial(const Material& m){_material = m;}

   void setValues(bool intersects, const Point& p, const Point& n, const Material& m, const
      Line& r, const Line& t)
   {_intersects = intersects; _point = p; _normal = n; _material = m; _reflectedRay = r; _transmittedRay =
   t;}

   void setValues(const Intersection& in)
   {_intersects = in._intersects; _point = in._point; _normal = in._normal; 
   _material = in._material; _reflectedRay = in._reflectedRay; _transmittedRay = in._transmittedRay;}

};
/*-----------------------------------------------*/
class RayObject
   /*
   PURPOSE: abstract class which serves as a base for
   all objects to be drawn in our ray-traced g_scene

   REMARK:
   */
{
protected:
   Point _position;
   Material _material;
public:
   RayObject(const Point& p, const Material& m )
   {_position = p; _material = m;}
   virtual void intersection(const Line& l, const Point& positionOffset, Intersection& inter) = 0; 
   // by overriding intersection in different ways control how rays hit objects in our g_scene

};
/*-----------------------------------------------*/
class Triangle : public RayObject
   /*
   PURPOSE: Triangle object are used as one of the basic building blocks for Shape's in our
   scen

   REMARK: Triangle's and Shape's are used according to a Composite design pattern to 
   define objects in our g_scene         
   */
{
private:
   Point _vertex0;
   Point _vertex1;
   Point _vertex2;

   Point _u;
   Point _v;
   Point _n;

   GLdouble _uv;
   GLdouble _uu;
   GLdouble _vv;
   GLdouble _denominator;

   bool _degenerate;

public:
   Triangle(const Point& p, const Material& m, const Point& p1, const Point& p2, const
      Point& p3) : RayObject(p,m)
   {
      _vertex0 = p1;
      _vertex1 = p2; 
      _vertex2 = p3;
      //compute intersection with plane of triangle
      _u = _vertex1 - _vertex0;
      _v = _vertex2 - _vertex0;
      _n = _u*_v;

      //handle last degenerates case by saying we don't intersect
      if(_n.length() < SMALL_NUMBER) _degenerate = true;
      else _degenerate = false;

      _n.normalize(); 

      _uv = _u & _v;
      _uu = _u & _u;
      _vv = _v & _v;

      _denominator = _uv*_uv - _uu*_vv;

      if( abs(_denominator) < SMALL_NUMBER) _degenerate = true;

   }

   void intersection(const Line& l, const Point& positionOffset, Intersection& inter);

};
/*-----------------------------------------------*/
class Shape : public RayObject
   /*
   PURPOSE: Shape's are either composite objects which represent a part of the g_scene
   to be ray-traced or our primary objects in which case they are used to model sphere's
   in our g_scene

   REMARK:  Triangle's and Shape's are used according to a Composite design pattern to 
   define objects in our g_scene                 
   */
{
protected:
   GLdouble _radius;
   bool _amSphere;
   bool _canIntersectOnlyOneSubObject;

   vector<RayObject *> _subObjects;

public:
   Shape() : RayObject(Point(0,0,0), Material())
   {_radius = 0; _amSphere = false;}
   Shape(Point p, Material m, GLdouble radius, bool a, bool c = false) : RayObject(p,m)
   {_radius = radius; _amSphere = a; _canIntersectOnlyOneSubObject = c; _subObjects.clear();}

   ~Shape();

   void setRadius(GLdouble r){_radius = r;}

   void addRayObject(RayObject *objects)
   {_subObjects.push_back(objects);}

   void intersection(const Line& l, const Point& positionOffset, Intersection& inter);


};
/*-----------------------------------------------*/
class Quad : public Shape
   /*
   PURPOSE: encapsulate information about quadrilaterals
   used in our g_scene. Quadralaterals are used for
   both quick tests for ray intersection as well as
   subobjects of more complicated objects in our g_scene

   REMARK:         
   */
{
public:
   Quad(Point p, Material m, Point p1, Point p2, Point p3, Point p4);

};
/*-----------------------------------------------*/
class Tetrahedron :  public Shape
   /*
   PURPOSE: encapsoluates information about
   tetrahedrons to be drawn in our g_scene (in this case just one)

   REMARK:         
   */

{
public:
   Tetrahedron(Point p, GLdouble edgeSize);
};
/*-----------------------------------------------*/
class Sphere :  public Shape
   /*
   PURPOSE: encapsoluates information about
   spheres to be drawn in our g_scene (in this case just one)

   REMARK:         
   */

{
public:
   Sphere(Point p, GLdouble radius);
};
/*-----------------------------------------------*/
class Cube :  public Shape
   /*
   PURPOSE: encapsoluates information about
   cubes to be drawn in our g_scene (in this case just one)

   REMARK:         
   */
{
public:
   Cube(Point p, GLdouble edgeSize);
};
/*-----------------------------------------------*/
class CheckerBoard :  public Shape
   /*
   PURPOSE: encapsoluates information about
   checkerboards to be drawn in our g_scene (in this case just one)

   REMARK:         
   */
{
private:
   Quad _boundingSquare; /* this Quad is used for a quick test to see if a
                         ray intersects our chessboard */
public:
   CheckerBoard(Point p);
   void intersection(const Line& l, const Point& positionOffset, Intersection& inter);

};

//   GLOBAL VARIABLES
GLsizei g_windowWidth = 500, g_windowHeight = 500; // used for size of window
GLsizei g_initX = 100, g_initY = 100; // used for initial position of window

Point g_lightPosition(0.0, 0.0, 0.0); /* although the ray tracer actually supports
                                      giving it a vector of lights, this program only makes use of one
                                      light which is placed at g_lightPosition The value is later changed from this
                                      default value to a value on the chess board.*/
Point g_lightColor(WHITE); // color of the light

Point g_whiteColor(WHITE); // some abbreviations for various colors
Point g_blackColor(BLACK);
Point g_redColor(RED);

Material g_whiteSquare(.1*g_whiteColor, .5*g_whiteColor, g_whiteColor, g_blackColor, 1); 
// some materials used by objects in  the g_scene
Material g_blackSquare(g_blackColor, .1*g_whiteColor, g_blackColor, g_blackColor, 1);
Material g_sphereMaterial(g_blackColor, .1*g_whiteColor, g_whiteColor, g_blackColor, 1);
Material g_tetrahedronMaterial(g_blackColor, g_blackColor, .1*g_whiteColor, g_whiteColor, 2.0/3.0);
Material g_cubeMaterial(.1*g_redColor, .4*g_redColor, g_redColor, g_blackColor, 1);

Shape g_scene(BOARD_POSITION, Material(), sqrt((double)3)*BOARD_HALF_SIZE, false); // global shape for whole g_scene

static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states
unsigned char kbPrevState[SDL_NUM_SCANCODES] = {0};
static int g_activeShader = 0;
map <string, int> boardMap;

/*
static float g_frustFovY = G_FRUST_MIN_FOV; // FOV in y direction

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

///////////////// END OF G L O B A L S ///////////////////////
*/


//   IMPLEMENTATIONS

//Triangle Class Implementations
/*-----------------------------------------------*/
void Triangle::intersection(const Line& ray, const Point& positionOffset, Intersection& inter)
   /*
   PURPOSE: used to fill in an Intersection object with information about how the supplied ray
   intersects with the current Triangle based on the given positionOffset Point vector.
   RECEIVES:
   ray -- ray to intersect with this Triangle
   positionOffset -- where in the overall g_scene this Triangle lives
   inter -- Intersection object to fill in with information about the if and where the ray
   intersects
   RETURNS: 

   REMARKS:
   The idea for the following way to test for triangle intersections was
   gotten from
   http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm

   Later noticed was in books already had like 3D Computer Graphics by Sam Buss.
   The coding of it is my own.

   Degenerate cases handle by saying there is no intersection
   */
{ 
   if(_degenerate)
   {
      inter.setIntersect(false);
      return;
   } 

   //get coordinates of triangle given our position
   Point position = _position + positionOffset;
   Point v0 = position + _vertex0;
   Point v1 = position + _vertex1;
   Point v2 = position + _vertex2;

   Point p0 = ray.startPoint();
   Point p1 = ray.endPoint();
   Point diffP = p1 - p0;
   GLdouble ndiffP = _n&diffP;

   //handle another degenerate case by saying we don't intersect
   if( abs(ndiffP) < SMALL_NUMBER)
   {
      inter.setIntersect(false);
      return;
   }   

   GLdouble m = (_n & (v0 - p0))/ (_n & diffP);

   if( m < SMALL_NUMBER) //if m is negative thenwe don't intersect
   {
      inter.setIntersect(false);
      return;
   }

   Point p = p0 + m*diffP; //intersection point with plane

   Point w = p - v0;

   //Now check if in triangle      
   GLdouble wu = w & _u;
   GLdouble wv = w & _v;

   GLdouble s = (_uv*wv - _vv*wu)/_denominator;
   GLdouble t = (_uv*wu - _uu*wv)/_denominator;

   if( s >= 0 && t >=0 && s + t <= 1) // intersect
   {

      diffP.normalize(); // now u is as in the book
      Point u=diffP;

      Point r =  u - (2*(u & _n))*_n; 
      Line reflected(p, p + r);

      //Transmitted vector calculated using thin lens equations from book
      GLdouble refractionRatio = _material.refraction();

      Point t(0.0, 0.0, 0.0);

      GLdouble cosThetai = u & _n;
      GLdouble modulus = 1 - refractionRatio*refractionRatio*( 1- cosThetai*cosThetai);

      if( modulus > 0)
      {
         GLdouble cosThetar = sqrt(modulus);
         t = refractionRatio*u - ( cosThetar + refractionRatio*cosThetai)*_n;
      }      

      Line transmitted(p, p + t);
      inter.setValues(true, p, _n, _material, reflected, transmitted);

   }
   else // don't intersect
   {
      inter.setIntersect(false);
   }
}
//Shape Class Implementations
/*-----------------------------------------------*/
Shape::~Shape()
   /*
   PURPOSE: destructs this shape and gets rid of any sub-object on it
   RECEIVES: nothing
   RETURNS: nothing
   REMARKS:    
   */
{
   for(size_t i = 0; i < _subObjects.size() ; i++)
   {
      delete _subObjects[i];
   }
}
/*-----------------------------------------------*/
void Shape::intersection(const Line& ray, const Point& positionOffset, Intersection& inter)
   /*
   PURPOSE: used to fill in an Intersection object with information about how the supplied ray
   intersects with the current Shape based on the given positionOffset Point vector.
   RECEIVES:
   ray -- ray to intersect with this Shape
   positionOffset -- where in the overall g_scene this Shape lives
   inter -- Intersection object to fill in with information about the if and where the ray
   intersects
   RETURNS: 
   REMARKS: note this Shape could be a composite object so we recurse through its _subObjects vector
   */
{
   Point u = ray.direction();
   Point p0 = ray.startPoint();
   Point position = _position + positionOffset;
   Point deltaP = position - p0;


   /* check for sphere intersection if radius is > 0.
   if radius ==0 assume we are dealing with a object which has
   subObjects which do the testing
   */
   if(_radius > 0 || _amSphere)
   {
      GLdouble uDeltaP = u & deltaP;
      GLdouble discriminant = uDeltaP*uDeltaP - (deltaP & deltaP) + _radius*_radius;

      GLdouble s = uDeltaP - sqrt(discriminant); //other solution is on far side of sphere

      if( discriminant < 0 || abs(s) < SMALL_NUMBER )
      {
         inter.setIntersect(false);
         return;
      }

      //calculate point of intersection

      Point p = p0 + s*u;
      Point directionP0 = p - position ;

      if(_amSphere)
      {
         if(s < SMALL_NUMBER) // if not in front of ray then don't intersect
         {
            inter.setIntersect(false);
            return;
         }

         //reflected vector calculated using equations from book
         Point n(directionP0);
         n.normalize();

         Point r =  u - (2*(u & n))*n; 
         Line reflected(p, p + r);

         //Transmitted vector calculated using thin lens equations from book
         Point t(0.0, 0.0, 0.0);
         GLdouble refractionRatio = _material.refraction();
         GLdouble cosThetai = u & n;
         GLdouble modulus = 1 - refractionRatio*refractionRatio*( 1- cosThetai*cosThetai);

         if( modulus > 0)
         {
            GLdouble cosThetar = sqrt(modulus);
            t = refractionRatio*u - ( cosThetar + refractionRatio*cosThetai)*n;
         }
         Line transmitted(p, p + t);
         inter.setValues(true, p, n, _material, reflected, transmitted);
      }
   }

   if(!_amSphere)
   {
      inter.setIntersect(false);
      Intersection interTmp;

      GLdouble minDistance = -1.0;
      GLdouble distanceTmp;
      size_t size = _subObjects.size();
      for(size_t i = 0; i < size; i++)
      {

         _subObjects[i]->intersection(ray, position, interTmp);

         if(interTmp.intersects())
         {
            Point directionCur = interTmp.point() - p0;
            distanceTmp = directionCur.length();
            if(distanceTmp < minDistance|| minDistance < 0.0)
            {
               minDistance = distanceTmp;
               inter.setValues(interTmp);
               if(_canIntersectOnlyOneSubObject) return;
            }
         }
      }
   }

}
//Quad Class Implementations
/*-----------------------------------------------*/
Quad::Quad(Point p, Material m, Point p1, Point p2, Point p3, Point p4) : Shape(p, m, 0, false, true)
   /*
   PURPOSE: constructs a Quad object (for quadralateral) at the given offset position made of
   the given material and with the supplied points
   RECEIVES:
   p -- position offset into our g_scene
   m -Material Quad is made out of
   p1, p2, p3, p4 - four local coordinate points (relative to p) that make up the Quad.
   RETURNS: 
   REMARKS: can only intersect one of two sub-triangles unless ray is same plane as Quad 
   */
{
   Point zero(0.0, 0.0, 0.0);

   addRayObject(new Triangle(zero, m, p1, p2, p3));
   addRayObject(new Triangle(zero, m, p1, p3, p4));

}
//Sphere Class Implementations
/*-----------------------------------------------*/
Sphere::Sphere(Point p, GLdouble r) : Shape(p, g_sphereMaterial, r, true)
   /*
   PURPOSE: constructs a Sphere object at the given offset position and radius in our Scene
   RECEIVES:
   p -- position offset into our g_scene
   r -- radius of Sphere
   RETURNS: a Sphere object
   REMARKS:  g_sphereMaterial is a global in this file.
   The constructor mainly just calls the base constructor with the appropriate material
   and with the flag for Shape telling the Shape that it is a non-composite
   sphere (the last paramter true to the Shape constructor)
   */
{

}
//Tetrahedron Class Implementations
/*-----------------------------------------------*/
Tetrahedron::Tetrahedron(Point p, GLdouble edgeSize) : Shape(p, g_tetrahedronMaterial, sqrt((double)3)*edgeSize/2,
   false)
   /*
   PURPOSE: constructs a Tetrahedron at the given offset position and edgeSize in our Scene
   RECEIVES:
   p -- position offset into our g_scene
   edgeSize -- size of an edge of our cube
   RETURNS: a Tetrahedron object
   REMARKS:  note g_tetrahedronMaterial is a global in this file.
   Since the HW description didn't say the tetrahedron was regular, we took the tetrahedron to be the
   one obtained by slicing the cube from a top corner through the diagonal of the bottom face
   */
{
   Point zero(0.0, 0.0, 0.0);
   GLdouble halfEdge = edgeSize/2;

   //bottom
   addRayObject(new Triangle(zero, g_tetrahedronMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, -halfEdge, -halfEdge), 
      Point(-halfEdge, -halfEdge, halfEdge)));
   //back
   addRayObject(new Triangle(zero, g_tetrahedronMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(-halfEdge, -halfEdge, halfEdge), 
      Point(-halfEdge, halfEdge, -halfEdge)));

   //left
   addRayObject(new Triangle(zero, g_tetrahedronMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(-halfEdge, halfEdge, -halfEdge), 
      Point(-halfEdge, -halfEdge, halfEdge)));


   //front
   addRayObject(new Triangle(zero, g_tetrahedronMaterial,
      Point(-halfEdge, -halfEdge, halfEdge),
      Point(halfEdge, -halfEdge, -halfEdge),
      Point(-halfEdge, halfEdge, -halfEdge)));

}
//Cube Class Implementations
/*-----------------------------------------------*/
Cube::Cube(Point p, GLdouble edgeSize) : Shape(p, g_cubeMaterial, sqrt((double)3)*edgeSize/2, false)
   /*
   PURPOSE: constructs a Cube at the given offset position and edgeSize in our Scene
   RECEIVES:
   p -- position offset into our g_scene
   edgeSize -- size of an edge of our cube
   RETURNS: a cube object
   REMARKS:  note g_cubeMaterial is a global in this file
   */
{
   GLdouble halfEdge = edgeSize/2;

   Point zero(0.0, 0.0, 0.0);
   //top
   addRayObject(new Quad(zero, g_cubeMaterial, Point(-halfEdge, halfEdge, -halfEdge), 
      Point(halfEdge, halfEdge, -halfEdge), 
      Point(halfEdge, halfEdge, halfEdge), 
      Point(-halfEdge, halfEdge, halfEdge)));

   //bottom
   addRayObject(new Quad(zero, g_cubeMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, -halfEdge, halfEdge), 
      Point(-halfEdge, -halfEdge, halfEdge)));

   //left
   addRayObject(new Quad(zero, g_cubeMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(-halfEdge, halfEdge, -halfEdge), 
      Point(-halfEdge, halfEdge, halfEdge), 
      Point(-halfEdge, -halfEdge, halfEdge)));
   //right
   addRayObject(new Quad(zero, g_cubeMaterial, Point(halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, halfEdge, -halfEdge), 
      Point(halfEdge, halfEdge, halfEdge), 
      Point(halfEdge, -halfEdge, halfEdge)));
   //back
   addRayObject(new Quad(zero, g_cubeMaterial, Point(-halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, -halfEdge, -halfEdge), 
      Point(halfEdge, halfEdge, -halfEdge), 
      Point(-halfEdge, halfEdge, -halfEdge)));

   //front
   addRayObject(new Quad(zero, g_cubeMaterial, Point(-halfEdge, -halfEdge, halfEdge), 
      Point(halfEdge, -halfEdge, halfEdge), 
      Point(halfEdge, halfEdge, halfEdge), 
      Point(-halfEdge, halfEdge, halfEdge)));

}
//CheckerBoard Class Implementations
/*-----------------------------------------------*/
CheckerBoard::CheckerBoard(Point p) : Shape(),

   _boundingSquare(p, Material(),  Point(- BOARD_HALF_SIZE, 0, - BOARD_HALF_SIZE), 
   Point( BOARD_HALF_SIZE, 0, - BOARD_HALF_SIZE), 
   Point( BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE), 
   Point(- BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE)) /*initialize
                                                 bounding square to be used for quick tests if ray intersects
                                                 chess board */
                                                 /*
                                                 PURPOSE: constructs a Checkerboard object at the supplied offset position
                                                 RECEIVES:
                                                 p - offset position to put chess board at
                                                 RETURNS: 
                                                 REMARKS: This constructor makes use of some constant in this file concerning the
                                                 chessboard.
                                                 */
{

}
/*-----------------------------------------------*/
void CheckerBoard::intersection(const Line& ray, const Point& positionOffset, Intersection& inter) 
   /*
   PURPOSE: used to calculate how a ray intersects with a Chessboard object
   RECEIVES:
   ray to bounce off of chessboard
   positionoOffset - offset vector for location of chessboard
   inter -- an Intersection object to be filled in with details on if and where ray intersects 
   with Chessboard
   RETURNS: nothing 
   REMARKS:    
   */
{

   _boundingSquare.intersection(ray, positionOffset, inter);

   if(inter.intersects())
   {
      Point p = inter.point() - positionOffset + Point(BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE);
      int squareSum = int(p.x()/SQUARE_EDGE_SIZE) + int(p.z()/SQUARE_EDGE_SIZE);

      if((squareSum & 1) == 0) 
      {
         inter.setMaterial(g_whiteSquare);
      }
      else
      {
         inter.setMaterial(g_blackSquare);
      }
   }

}

//CLASSLESS FUNCTIONS

/*-----------------------------------------------*/
Point operator*(GLdouble scalar, const Point& p)
   /*
   PURPOSE: multiplies the supplied point vector by the scalar amount
   RECEIVES: 
   p - Point vector (x,y,z)
   scalar -- the scalar `a' to multiply by
   RETURNS: 
   the point vector
   (a*x, a*y, a*z)
   REMARKS:    
   */
{
   return Point(scalar*p._x, scalar*p._y, scalar*p._z);
}
/*-----------------------------------------------*/
Point operator*(const Point& p, GLdouble scalar)
   /*
   PURPOSE: multiplies the supplied point vector by the scalar amount
   RECEIVES: 
   p - Point vector (x,y,z)
   scalar -- the scalar `a' to multiply by
   RETURNS: 
   the point vector
   (a*x, a*y, a*z)
   REMARKS:    
   */
{
   return Point(scalar*p._x, scalar*p._y, scalar*p._z);
}
/*-----------------------------------------------*/
Point randomUnit()
   /*
   PURPOSE: generate a random vector of length 1
   RECEIVES: Nothing
   RETURNS: Nothing
   REMARKS:
   */
{

   //generate random point within unit sphere
   Point vec(0.0, 0.0, 0.0);

   while(vec.isZero())
   {
      vec = Point(double(rand())/(RAND_MAX+1.0) - .5, 
         double(rand())/(RAND_MAX+1.0) - .5,
         double(rand())/(RAND_MAX+1.0) - .5);
   }
   vec.normalize(); //push out to unit sphere

   return vec;
}
/*-----------------------------------------------*/
inline GLdouble attenuation(GLdouble distance)
   /*
   PURPOSE: calculates how much light intensities will decay with distance
   RECEIVES: distance -- to use
   RETURNS: decimal value between 0 and 1 by which an intensity at the given distance
   would be reduced
   REMARKS:    
   */
{

   return ATTENUATION_FACTOR/(ATTENUATION_FACTOR + distance*distance);
}
/*-----------------------------------------------*/
void rayTraceRay(Shape& g_scene, vector<Light> lights, const Line& ray, Point& color, unsigned int 
   depth)
   /*
   PURPOSE: Does ray tracing of a single ray in a g_scene according to the supplied lights to the perscribed
   depth
   RECEIVES:
   g_scene -- Shape to do ray-tracing one
   lights -- Light's which are lighting the g_scene
   ray -- to be used for ray-tracing consists of two points (starting point to do ray-tracing from plus
   another point which together give the direction of the initial ray.)
   color -- used to store the color returned by doing the ray tracing
   depth -- in terms of tree of sub-rays we calculate
   RETURNS:  Nothing
   REMARKS:    
   */
{
   Intersection intersection;
   g_scene.intersection(ray, Point(0.0, 0.0, 0.0), intersection);

   if(!intersection.intersects()) return;

   Point pt = intersection.point();
   Material material = intersection.material();
   Line reflectedRay = intersection.reflectedRay();
   Line transmittedRay = intersection.transmittedRay();

   Line shadowRay;
   Point lColor;

   size_t size = lights.size();
   for(size_t i = 0 ; i < size; i++)
   {
      shadowRay.set(pt, lights[i].position());
      Intersection shadowIntersection;

      g_scene.intersection(shadowRay, Point(0.0, 0.0, 0.0), shadowIntersection );

      if(!shadowIntersection.intersects() || !shadowIntersection.material().transparency().isZero())
      {
         lColor = attenuation(shadowRay.length())*lights[i].color();
         color += (material.ambient()% lColor) +
            abs(intersection.normal() & shadowRay.direction())*(material.diffuse()% lColor) +
            abs(ray.direction() & reflectedRay.direction())*(material.specular()% lColor);
      }
   }

   if(depth > 0)
   {
      Point transmittedColor(0.0, 0.0, 0.0);
      Point reflectedColor(0.0, 0.0, 0.0);

      Point transparency = material.transparency();
      Point opacity = Point(1.0, 1.0, 1.0) - transparency;

      if(!transparency.isZero() && transparency.length() > SMALL_NUMBER) //if not transparent then don't send ray
      {
         rayTraceRay(g_scene, lights, transmittedRay, transmittedColor, depth - 1);
         color += (transparency % transmittedColor);
      }
      if(!opacity.isZero()) // if completely transparent don't send reflect ray
      {
         rayTraceRay(g_scene, lights, reflectedRay, reflectedColor, depth - 1);
         color += (opacity % reflectedColor);
      }
   } 
}
/*-----------------------------------------------*/
void rayTraceScreen(Shape& g_scene, vector<Light>& lights, Point camera, Point lookAt, Point up, int
   bottomX, int bottomY, int width, int height)
   /*
   PURPOSE: Does the ray-tracing of a shape according to the supplied lights, camera dimension and screen dimensions

   RECEIVES:
   g_scene --  Shape to be ray-traced (Shapes use the Composite pattern so are made of sub-shapes
   light -- a vector of Light's used to light the g_scene
   camera -- location of the viewing position
   lookat -- where one is looking at from this position
   up -- what direction is up from this position
   bottomX -- how far to the left from the lookat point is the start of the screen 
   bottomy -- how far down from the lookat point is the start of the screen
   width -- width of the screen
   height -- height of the screen
   RETURNS:  Nothing
   REMARKS:    
   */
{
   Point lookDirection = lookAt - camera;
   Point right = lookDirection * up;

   right.normalize();
   Point rightOffset = width*right;

   up = right*lookDirection;
   up.normalize();

   Point screenPt = lookAt + bottomX*right + bottomY*up;

   Line ray;
   Point color(0.0, 0.0, 0.0);
   Point avgColor(0.0, 0.0, 0.0);

   Point weightedColor(0.0, 0.0, 0.0);
   Point oldWeightedColor(0.0, 0.0, 0.0);
   GLdouble k;

   glBegin(GL_POINTS);
   for(int j = 0; j < height; j++)
   {
      for(int i = 0; i < width; i++)
      {
         for(k = 0.0; k < SUPER_SAMPLE_NUMBER; k++)
         {
            ray.set(camera, screenPt + .5*randomUnit() );

            color.set(0.0, 0.0, 0.0);

            rayTraceRay(g_scene, lights, ray, color, MAX_DEPTH);

            oldWeightedColor = (k + 1.0)*avgColor;

            avgColor += color;

            weightedColor = k*avgColor;
            if((weightedColor - oldWeightedColor).length() < SMALL_NUMBER*k*(k+1)) break;

         } 

         avgColor /= k;            
         glColor3d(avgColor.x(), avgColor.y(), avgColor.z());

         glVertex2i(i, j);
         screenPt += right;

         //cout << "Color = [" << avgColor.x() << ", " << avgColor.y() << ", " << avgColor.z() << "]" << endl;
         //cout << "Pixel = [" << i << ", " << j << "]" << endl;
      }
      screenPt -= rightOffset;
      screenPt += up;
   }
   glEnd();
}
/*-----------------------------------------------*/
Point  convertStringCoordinate(string coordString)
   /*
   PURPOSE: Converts a string of two characters of the form: 
   letter row + number color (for example, b4) into coordinates for a piece of 
   our ray tracing work. 
   RECEIVES: coordString -- string to convert
   RETURNS: Nothing
   REMARKS:    
   */
{
   Point firstSquare(-BOARD_EDGE_SIZE/2, 0.0, BOARD_EDGE_SIZE/2);

   Point rowOffset(0.0, 0.0, -(double(coordString[0] - 'a') + .5)*SQUARE_EDGE_SIZE);
   //negative because farther back == higher row number
   Point colOffset((double(coordString[1] - '0' - 1) + .5)*SQUARE_EDGE_SIZE, 0.0, 0.0);
   Point heightOffset(0.0, 1.5*SQUARE_EDGE_SIZE, 0.0);

   Point square = firstSquare + rowOffset + colOffset + heightOffset;

   return square;
}
/*-----------------------------------------------*/
static void initShaders()
{
   /* PURPOSE:		Initializes Shaders to be used 
   RECEIVES:	 
   RETURNS:     
   REMARKS:     
   */

   g_shaderStates.resize(G_NUM_SHADERS);
   for (int i = 0; i < G_NUM_SHADERS; ++i) 
   {
      if (G_GL2_COMPATIBLE) 
      {
         g_shaderStates[i].reset(new ShaderState(G_SHADER_FILES_GL2[i][0],
            G_SHADER_FILES_GL2[i][1]));
      } 
      else 
      {
         g_shaderStates[i].reset(new ShaderState(G_SHADER_FILES[i][0],
            G_SHADER_FILES[i][1]));
      }
   }
}
/*-----------------------------------------------*/
static void initGLState()
{
   /*	PURPOSE:		Initializes OpenGL 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   glClearColor(0.0, 0.0, 0.0, 1.0);
   glViewport(0, 0, g_windowWidth, g_windowHeight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, g_windowWidth, 0, g_windowHeight); // Upside-Down
}
/*-----------------------------------------------*/
void MySdlApplication::initScene()
   /*	PURPOSE: gets locations of objects from users
   sets the background color to black,
   sets up display list for coordinate axes,
   initializes extrusions vector
   RECEIVES: Nothing
   RETURNS: Nothing
   REMARKS:    
   */
{

   Point boardPosition(0.0, 0.0, 0.0);
   CheckerBoard *checkerBoard = new CheckerBoard(boardPosition);
   g_scene.addRayObject(checkerBoard);

   string tmp;

   cout << "Please enter the position of the light:\n";
   //cin >> tmp;
   tmp ="b6"; //commented values are nice values to demonstrate the ray tracer
   g_lightPosition = Point(BOARD_POSITION) + Point(0.0, 3.5*SQUARE_EDGE_SIZE, 0.0) + convertStringCoordinate(tmp);
   //with what convertStringCoordinate gives makes 5 squares above board

   cout << "Please enter the position of the tetrahedron:\n";
   //cin >>tmp;
   tmp = "b4";
   Tetrahedron *tetrahedron = new Tetrahedron(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE);
   g_scene.addRayObject(tetrahedron);

   cout << "Please enter the position of the sphere:\n";
   //cin >> tmp;
   tmp = "d7";
   Sphere *sphere = new Sphere(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE/2);
   g_scene.addRayObject(sphere);

   cout << "Please enter the position of the cube:\n";
   //cin >> tmp;
   tmp = "a7";
   Cube *cube = new Cube(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE);   
   g_scene.addRayObject(cube);

}
/*-----------------------------------------------*/
void MySdlApplication::initScene2()
   /*	PURPOSE: gets locations of objects from users
   sets the background color to black,
   sets up display list for coordinate axes,
   initializes extrusions vector
   RECEIVES: Nothing
   RETURNS: Nothing
   REMARKS:    
   */
{
   // Initialize Checkerboard
   Point boardPosition(0.0, 0.0, 0.0);
   CheckerBoard *checkerBoard = new CheckerBoard(boardPosition);
   g_scene.addRayObject(checkerBoard);
   boardMap.clear();

   string tmp;
   bool isFinished = false;
   cout << endl;

   while (!isFinished)
   {
      bool hasAnswered = false;

      while (!hasAnswered)
      {
         cout << "Please select the type of object to add:" << endl;
         cout << "(a) light, (b) tetrahedron, (c) cube, (d) sphere, (e) cylinder, (f) cone" << endl;
         cin >> tmp;

         hasAnswered = true;
         if (tmp.size() > 1)
            hasAnswered = false;

         int type = tmp[0] - 'a';
         if (type >= 0 && type < 6)
         {
            cout << "Please enter the position: (A1-H8)" << endl;
            cin >> tmp;

            boardMap[tmp] = type;
         }
         else
            hasAnswered = false;
      }

      hasAnswered = false;
      while (!hasAnswered)
      {
         cout << "Would you like to add another object? (yes/no)" << endl;
         cin >> tmp;

         if (!tmp.compare("no"))
         {
            isFinished = true;
            hasAnswered = true;
         }
         else if (!tmp.compare("yes"))
            hasAnswered = true;
      }
   }

   loadScene();
}
/*-----------------------------------------------*/
void MySdlApplication::loadScene()
   /*	PURPOSE: Loads objects into scene
   RECEIVES: Nothing
   RETURNS: Nothing
   REMARKS:    
   */
{
   map<string, int>::iterator iter;

   for (iter = boardMap.begin(); iter != boardMap.end(); ++iter)
   {
      string position = iter->first;
      int type = iter->second;

      if (type == LIGHT)
      {
         g_lightPosition = Point(BOARD_POSITION) + Point(0.0, 3.5*SQUARE_EDGE_SIZE, 0.0) + convertStringCoordinate(position);
      }
      else if (type == TETRAHEDRON)
      {
         Tetrahedron *tetrahedron = new Tetrahedron(convertStringCoordinate(position), SQUARE_EDGE_SIZE);
         g_scene.addRayObject(tetrahedron);
      }
      else if (type == CUBE)
      {
         Cube *cube = new Cube(convertStringCoordinate(position), SQUARE_EDGE_SIZE);   
         g_scene.addRayObject(cube);
      }
      else if (type == SPHERE)
      {
         Sphere *sphere = new Sphere(convertStringCoordinate(position), SQUARE_EDGE_SIZE/2);
         g_scene.addRayObject(sphere);
      }
      else if (type == CYLINDER)
      {
      }
      else if (type == CONE)
      {
      }
   }
}
/*-----------------------------------------------*/
void draw()
   /*
   PURPOSE: Used to craw the complete ray-traced chessboard

   RECEIVES: Nothing
   RETURNS:  Nothing
   REMARKS:    
   */
{
   glClear(GL_COLOR_BUFFER_BIT);

   vector<Light> lights;

   lights.push_back(Light(g_lightColor, g_lightPosition));

   Point camera(CAMERA_POSITION);
   Point lookAt(LOOK_AT_VECTOR);   
   Point up(UP_VECTOR); 

   rayTraceScreen(g_scene, lights, camera, lookAt, up, -g_windowWidth/2, -g_windowHeight/2, g_windowWidth, g_windowHeight);

   glFlush();
}
/*-----------------------------------------------*/
void reshape(int newWidth, int newHeight)
   /*
   PURPOSE: To redraw g_scene when  window gets resized.
   RECEIVES: newWidth -- new window x size
   newHeight -- new window y size
   RETURNS: Nothing 
   REMARKS:    
   */
{
   glViewport(0, 0, newWidth, newHeight);

   g_windowWidth = newWidth;
   g_windowHeight = newHeight;

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, g_windowWidth, 0, g_windowHeight);

   glClear(GL_COLOR_BUFFER_BIT);
}
/*-----------------------------------------------*/
void MySdlApplication::keyboard()
{
   /* PURPOSE:		Handles keyboard presses by user 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		Player related controls are handled through player class
   */

   if (KB_STATE[SDL_SCANCODE_ESCAPE])
   {
      g_running = false;
   }
}
/*-----------------------------------------------*/
void MySdlApplication::onLoop()
{
   /*	PURPOSE:		Handles function calls that need to run once per SDL loop 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   // Logic goes here
   keyboard();
}
/*-----------------------------------------------*/
void MySdlApplication::onRender()
{
   /*	PURPOSE:		Handles all graphics related calls once per SDL loop 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   // All draw calls go here
   //glUseProgram(g_shaderStates[g_activeShader]->program);
   //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   // clear framebuffer color&depth

   draw();

   SDL_GL_SwapWindow(g_display);
   checkGlErrors();
}
/*-----------------------------------------------*/
int MySdlApplication::onExecute()
{
   /*	PURPOSE:		Main function loop of MySdlApplication 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   if(onInit() == false) 
      return -1;

   SDL_Event Event;
   while(g_running) 
   {
      memcpy (kbPrevState, KB_STATE, sizeof( kbPrevState ));

      while(SDL_PollEvent(&Event)) 
      {
         onEvent(&Event);
      }

      onLoop();
      onRender();
   }

   onCleanup();

   return 0;
}
/*-----------------------------------------------*/
bool MySdlApplication::onInit()
{
   /*	PURPOSE:		Initializes SDL 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   //if(SDL_Init(SDL_INIT_EVERYTHING) < 0) 
   if(SDL_Init(SDL_INIT_VIDEO) < 0)
   {
      return false;
   }
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
   //SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);

   /* Turn on double buffering with a 24bit Z buffer.
   * You may need to change this to 16 or 32 for your system */
   SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
   SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

   if((g_display = SDL_CreateWindow("Ray Trace",
      g_initX, g_initY, g_windowWidth, g_windowHeight,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE)) == NULL) 
   {
      return false;
   }

   /* Create our opengl context and attach it to our window */
   SDL_GLContext maincontext = SDL_GL_CreateContext(g_display);
   /* This makes our buffer swap syncronized with the 
   monitor's vertical refresh */
   SDL_GL_SetSwapInterval(1);

   GLenum glewError = glewInit();
   if( glewError != GLEW_OK ) 
   {
      SDL_Quit();
      return 1;
   }
   if( !GLEW_VERSION_1_5 ) 
   {
      SDL_Quit();
      return 1;
   }

   cout << (G_GL2_COMPATIBLE ?
      "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3")
      << endl;
   if ((!G_GL2_COMPATIBLE) && !GLEW_VERSION_3_0)
      throw runtime_error(
      "Error: does not support OpenGL Shading Language v1.3");
   else if (G_GL2_COMPATIBLE && !GLEW_VERSION_2_0)
      throw runtime_error(
      "Error: does not support OpenGL Shading Language v1.0");

   //initShaders();
   initScene2();
   initGLState();

   KB_STATE = SDL_GetKeyboardState(NULL);

   return true;
}
/*-----------------------------------------------*/
void MySdlApplication::onEvent(SDL_Event* event) 
{
   /*	PURPOSE:		Handles SDL events 
   RECEIVES:	event - SDL Event to be handled 
   RETURNS:		 
   REMARKS:		 
   */

   Uint32 type = event->type;

   if (type == SDL_QUIT)
      g_running = false;
   /*
   else if (type == SDL_MOUSEBUTTONDOWN)
   mouse(event->button);
   else if (type == SDL_MOUSEBUTTONUP)
   mouse(event->button);
   else if (type == SDL_MOUSEMOTION)
   motion(event->motion.x, event->motion.y);
   */
   else if (type == SDL_WINDOWEVENT)
      if (event->window.event == SDL_WINDOWEVENT_RESIZED)
         reshape(event->window.data1,event->window.data2);
}
/*-----------------------------------------------*/
void MySdlApplication::onCleanup()
{
   /*	PURPOSE:		Everything to be done before program closes 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   SDL_Quit();
}
/*-----------------------------------------------*/
MySdlApplication::MySdlApplication()
{
   /*	PURPOSE:		Constructor for MySdlApplication 
   RECEIVES:	 
   RETURNS:		 
   REMARKS:		 
   */

   g_running = true;
}
/*-----------------------------------------------*/
int main (int argc, char **argv) 
{
   /*	PURPOSE:		Main function of MySdlApplication 
   RECEIVES:	argc - number of arguments passed
   argv - array of arguments
   RETURNS:		int - Whether or not program ran sucessfully
   REMARKS:		
   */

   MySdlApplication application;
   return application.onExecute();

} 


//Coding Guidelines template (REMOVE before Submission)


/*-----------------------------------------------*/


/*	PURPOSE:		What does this function do? (must be present) 
RECEIVES:	List every argument name and explain each argument. 
(omit if the function has no arguments) 
RETURNS:		Explain the value returned by the function. 
(omit if the function returns no value) 
REMARKS:		Explain any special preconditions or postconditions. 
See example below. (omit if function is unremarkable) 
*/