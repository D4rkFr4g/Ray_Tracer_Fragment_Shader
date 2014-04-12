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
#include "GPoint.h"

//   CONSTANTS
/*
static const bool G_GL2_COMPATIBLE = false;
static const float G_FRUST_MIN_FOV = 60.0;  //A minimal of 60 degree field of view
static const unsigned char* KB_STATE = NULL;

static const int G_NUM_SHADERS = 2;
static const char * const G_SHADER_FILES[G_NUM_SHADERS][2] = 
{
	{"./Shaders/basic-gl3.vshader", "./Shaders/diffuse-gl3.fshader"},
	{"./Shaders/basic-gl3.vshader", "./Shaders/solid-gl3.fshader"}
};
static const char * const G_SHADER_FILES_GL2[G_NUM_SHADERS][2] = 
{
	{"./Shaders/basic-gl2.vshader", "./Shaders/diffuse-gl2.fshader"},
	{"./Shaders/basic-gl2.vshader", "./Shaders/solid-gl2.fshader"}
};
*/

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

/*
// Global variables
int g_windowWidth = 640;
int g_windowHeight = 480;
unsigned char kbPrevState[SDL_NUM_SCANCODES] = {0};

static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

static float g_frustFovY = G_FRUST_MIN_FOV; // FOV in y direction

static shared_ptr<GlTexture> g_tex0, g_tex1, g_tex2, g_tex3;

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// --------- Scene

static Cvec3 g_light1(2.0, 3.0, 14.0);
static Cvec3 g_light2(-2000.0, -3000.0, -5000.0);
// define light positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.0, 0.0)); // Initialized here but set in initCamera()
static RigTForm g_eyeRbt = g_skyRbt;
static Cvec3f g_objectColors[1] = {Cvec3f(1, 0, 0)};

///////////////// END OF G L O B A L S ///////////////////////

*/

// PROTOTYPES
class RayObject;

/*
   CLASS DEFINITIONS
*/

/*-----------------------------------------------*/


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
      GPoint _color;
      GPoint _position;

   public:
     
      Light(const GPoint& c, const GPoint& p){ _color = c; _position = p;}
      GPoint color(){return _color;}
      GPoint position(){return _position;}  
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
      GPoint _startPt;
      GPoint _endPt;

   public:
      Line(){_startPt.set(0.0, 0.0, 0.0); _endPt.set(0.0, 0.0, 0.0);}
      Line(const GPoint& p1, const GPoint& p2) {_startPt = p1; _endPt = p2;}
      
      void set(const GPoint& p1, const GPoint& p2){_startPt = p1; _endPt = p2;}

      GPoint startPoint() const {return _startPt;}
      GPoint endPoint() const {return _endPt;}

      GPoint direction() const
      {
         GPoint p =  _endPt - _startPt;
         p.normalize();
         return p;
      }
      
      GLdouble length() const
      {
         GPoint p = _endPt - _startPt;
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
      GPoint _ambient;
      GPoint _diffuse;
      GPoint _specular;
      GPoint _transparency;
      GLdouble _refraction;
            
      
   public:
      Material()
         { _ambient.set(0.0, 0.0, 0.0); _diffuse =_ambient; _specular = _ambient; _transparency = _ambient;
_refraction = 1;}
      Material(const GPoint& a, const GPoint& d, const GPoint& s, const GPoint& t, GLdouble r)
         {_ambient = a; _diffuse = d; _specular = s; _transparency = t; _refraction = r;}
      Material(const Material& m)
         {_ambient = m._ambient; _diffuse = m._diffuse; _specular = m._specular; 
          _transparency = m._transparency; _refraction = m._refraction;}


      GPoint ambient(){ return _ambient;}
      GPoint diffuse(){ return _diffuse;}
      GPoint specular(){ return _specular;}
      GPoint transparency(){ return _transparency;}
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
      GPoint _point;
      GPoint _normal;

      Material _material;
      Line _reflectedRay;
      Line _transmittedRay;
            
   public:
      Intersection(){}
      Intersection(bool intersects, const GPoint& p, const GPoint& n, const Material& m, const Line& 
r, const Line& t)
         {_intersects = intersects; _point = p; _normal = n; _material = m; _reflectedRay = r; _transmittedRay =
t;}
      
         
      bool intersects(){return _intersects;}

      GPoint point(){return _point;}
      GPoint normal(){return _normal;}
      
      Material material(){return _material;}
      
      Line reflectedRay(){return _reflectedRay;}
      Line transmittedRay(){return _transmittedRay;}
      
      void setIntersect(bool i){_intersects = i;}
      void setMaterial(const Material& m){_material = m;}
      
      void setValues(bool intersects, const GPoint& p, const GPoint& n, const Material& m, const
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
      GPoint _position;
      Material _material;
   public:
      RayObject(const GPoint& p, const Material& m )
         {_position = p; _material = m;}
      virtual void intersection(const Line& l, const GPoint& positionOffset, Intersection& inter) = 0; 
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
      GPoint _vertex0;
      GPoint _vertex1;
      GPoint _vertex2;
      
      GPoint _u;
      GPoint _v;
      GPoint _n;
            
      GLdouble _uv;
      GLdouble _uu;
      GLdouble _vv;
      GLdouble _denominator;

      bool _degenerate;
      
   public:
      Triangle(const GPoint& p, const Material& m, const GPoint& p1, const GPoint& p2, const
GPoint& p3) : RayObject(p,m)
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
      
      void intersection(const Line& l, const GPoint& positionOffset, Intersection& inter);
        
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
      Shape() : RayObject(GPoint(0,0,0), Material())
         {_radius = 0; _amSphere = false;}
      Shape(GPoint p, Material m, GLdouble radius, bool a, bool c = false) : RayObject(p,m)
         {_radius = radius; _amSphere = a; _canIntersectOnlyOneSubObject = c; _subObjects.clear();}
 
      ~Shape();
      
      void setRadius(GLdouble r){_radius = r;}
      
      void addRayObject(RayObject *objects)
         {_subObjects.push_back(objects);}
      
      void intersection(const Line& l, const GPoint& positionOffset, Intersection& inter);

        
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
      Quad(GPoint p, Material m, GPoint p1, GPoint p2, GPoint p3, GPoint p4);
      
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
      Tetrahedron(GPoint p, GLdouble edgeSize);
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
      Sphere(GPoint p, GLdouble radius);
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
      Cube(GPoint p, GLdouble edgeSize);
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
      CheckerBoard(GPoint p);
      void intersection(const Line& l, const GPoint& positionOffset, Intersection& inter);
      
};



//   GLOBALS

GLsizei g_windowWidth = 500, g_windowHeight = 500; // used for size of window
GLsizei g_initX = 50, g_initY = 50; // used for initial position of window

GPoint g_lightPosition(0.0, 0.0, 0.0); /* although the ray tracer actually supports
   giving it a vector of lights, this program only makes use of one
   light which is placed at g_lightPosition The value is later changed from this
   default value to a value on the chess board.*/
GPoint g_lightColor(WHITE); // color of the light

GPoint g_whiteColor(WHITE); // some abbreviations for various colors
GPoint g_blackColor(BLACK);
GPoint g_redColor(RED);


Material g_whiteSquare(.1*g_whiteColor, .5*g_whiteColor, g_whiteColor, g_blackColor, 1); 
   // some materials used by objects in  the g_scene
Material g_blackSquare(g_blackColor, .1*g_whiteColor, g_blackColor, g_blackColor, 1);
Material g_sphereMaterial(g_blackColor, .1*g_whiteColor, g_whiteColor, g_blackColor, 1);
Material g_tetrahedronMaterial(g_blackColor, g_blackColor, .1*g_whiteColor, g_whiteColor, 2.0/3.0);
Material g_cubeMaterial(.1*g_redColor, .4*g_redColor, g_redColor, g_blackColor, 1);


Shape g_scene(BOARD_POSITION, Material(), sqrt((double)3)*BOARD_HALF_SIZE, false); // global shape for whole g_scene

/*
   IMPLEMENTATIONS
*/

//Triangle Class Implementations

/*-----------------------------------------------*/
void Triangle::intersection(const Line& ray, const GPoint& positionOffset, Intersection& inter)
/*
PURPOSE: used to fill in an Intersection object with information about how the supplied ray
   intersects with the current Triangle based on the given positionOffset GPoint vector.
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
   GPoint position = _position + positionOffset;
   GPoint v0 = position + _vertex0;
   GPoint v1 = position + _vertex1;
   GPoint v2 = position + _vertex2;
   
   GPoint p0 = ray.startPoint();
   GPoint p1 = ray.endPoint();
   GPoint diffP = p1 - p0;
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
      
   GPoint p = p0 + m*diffP; //intersection point with plane
   
   GPoint w = p - v0;
   
   //Now check if in triangle      
   GLdouble wu = w & _u;
   GLdouble wv = w & _v;
   
   GLdouble s = (_uv*wv - _vv*wu)/_denominator;
   GLdouble t = (_uv*wu - _uu*wv)/_denominator;
   
   if( s >= 0 && t >=0 && s + t <= 1) // intersect
   {

      diffP.normalize(); // now u is as in the book
      GPoint u=diffP;
      
      GPoint r =  u - (2*(u & _n))*_n; 
      Line reflected(p, p + r);
      
      //Transmitted vector calculated using thin lens equations from book
      GLdouble refractionRatio = _material.refraction();

      GPoint t(0.0, 0.0, 0.0);
      
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
void Shape::intersection(const Line& ray, const GPoint& positionOffset, Intersection& inter)
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
   GPoint u = ray.direction();
   GPoint p0 = ray.startPoint();
   GPoint position = _position + positionOffset;
   GPoint deltaP = position - p0;


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

      GPoint p = p0 + s*u;
      GPoint directionP0 = p - position ;
   
      if(_amSphere)
      {
         if(s < SMALL_NUMBER) // if not in front of ray then don't intersect
         {
            inter.setIntersect(false);
            return;
         }
         
         //reflected vector calculated using equations from book
         GPoint n(directionP0);
         n.normalize();
         
         GPoint r =  u - (2*(u & n))*n; 
         Line reflected(p, p + r);
         
         //Transmitted vector calculated using thin lens equations from book
         GPoint t(0.0, 0.0, 0.0);
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
            GPoint directionCur = interTmp.point() - p0;
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
Quad::Quad(GPoint p, Material m, GPoint p1, GPoint p2, GPoint p3, GPoint p4) : Shape(p, m, 0, false, true)
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
   GPoint zero(0.0, 0.0, 0.0);
   
   addRayObject(new Triangle(zero, m, p1, p2, p3));
   addRayObject(new Triangle(zero, m, p1, p3, p4));
      
}

//Sphere Class Implementations

/*-----------------------------------------------*/
Sphere::Sphere(GPoint p, GLdouble r) : Shape(p, g_sphereMaterial, r, true)
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
Tetrahedron::Tetrahedron(GPoint p, GLdouble edgeSize) : Shape(p, g_tetrahedronMaterial, sqrt((double)3)*edgeSize/2,
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
   GPoint zero(0.0, 0.0, 0.0);
   GLdouble halfEdge = edgeSize/2;

   //bottom
   addRayObject(new Triangle(zero, g_tetrahedronMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, -halfEdge, -halfEdge), 
                     GPoint(-halfEdge, -halfEdge, halfEdge)));
    //back
    addRayObject(new Triangle(zero, g_tetrahedronMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(-halfEdge, -halfEdge, halfEdge), 
                     GPoint(-halfEdge, halfEdge, -halfEdge)));

   //left
   addRayObject(new Triangle(zero, g_tetrahedronMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(-halfEdge, halfEdge, -halfEdge), 
                     GPoint(-halfEdge, -halfEdge, halfEdge)));


   //front
   addRayObject(new Triangle(zero, g_tetrahedronMaterial,
                     GPoint(-halfEdge, -halfEdge, halfEdge),
                     GPoint(halfEdge, -halfEdge, -halfEdge),
                     GPoint(-halfEdge, halfEdge, -halfEdge)));
 
}

//Cube Class Implementations

/*-----------------------------------------------*/
Cube::Cube(GPoint p, GLdouble edgeSize) : Shape(p, g_cubeMaterial, sqrt((double)3)*edgeSize/2, false)
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

   GPoint zero(0.0, 0.0, 0.0);
   //top
   addRayObject(new Quad(zero, g_cubeMaterial, GPoint(-halfEdge, halfEdge, -halfEdge), 
                     GPoint(halfEdge, halfEdge, -halfEdge), 
                     GPoint(halfEdge, halfEdge, halfEdge), 
                     GPoint(-halfEdge, halfEdge, halfEdge)));

   //bottom
   addRayObject(new Quad(zero, g_cubeMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, -halfEdge, halfEdge), 
                     GPoint(-halfEdge, -halfEdge, halfEdge)));

   //left
   addRayObject(new Quad(zero, g_cubeMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(-halfEdge, halfEdge, -halfEdge), 
                     GPoint(-halfEdge, halfEdge, halfEdge), 
                     GPoint(-halfEdge, -halfEdge, halfEdge)));
   //right
   addRayObject(new Quad(zero, g_cubeMaterial, GPoint(halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, halfEdge, -halfEdge), 
                     GPoint(halfEdge, halfEdge, halfEdge), 
                     GPoint(halfEdge, -halfEdge, halfEdge)));
    //back
    addRayObject(new Quad(zero, g_cubeMaterial, GPoint(-halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, -halfEdge, -halfEdge), 
                     GPoint(halfEdge, halfEdge, -halfEdge), 
                     GPoint(-halfEdge, halfEdge, -halfEdge)));
                     
   //front
    addRayObject(new Quad(zero, g_cubeMaterial, GPoint(-halfEdge, -halfEdge, halfEdge), 
                     GPoint(halfEdge, -halfEdge, halfEdge), 
                     GPoint(halfEdge, halfEdge, halfEdge), 
                     GPoint(-halfEdge, halfEdge, halfEdge)));
   
}

//CheckerBoard Class Implementations

/*-----------------------------------------------*/
CheckerBoard::CheckerBoard(GPoint p) : Shape(),

   _boundingSquare(p, Material(),  GPoint(- BOARD_HALF_SIZE, 0, - BOARD_HALF_SIZE), 
                     GPoint( BOARD_HALF_SIZE, 0, - BOARD_HALF_SIZE), 
                     GPoint( BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE), 
                     GPoint(- BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE)) /*initialize
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
void CheckerBoard::intersection(const Line& ray, const GPoint& positionOffset, Intersection& inter) 
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
      GPoint p = inter.point() - positionOffset + GPoint(BOARD_HALF_SIZE, 0, BOARD_HALF_SIZE);
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
GPoint operator*(GLdouble scalar, const GPoint& p)
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
   return GPoint(scalar*p._x, scalar*p._y, scalar*p._z);
}

/*-----------------------------------------------*/
GPoint operator*(const GPoint& p, GLdouble scalar)
/*
PURPOSE: multiplies the supplied point vector by the scalar amount
RECEIVES: 
   p - GPoint vector (x,y,z)
   scalar -- the scalar `a' to multiply by
RETURNS: 
   the point vector
   (a*x, a*y, a*z)
REMARKS:    
*/
{
   return GPoint(scalar*p._x, scalar*p._y, scalar*p._z);
}

/*-----------------------------------------------*/
GPoint randomUnit()
/*
PURPOSE: generate a random vector of length 1
RECEIVES: Nothing
RETURNS: Nothing
REMARKS:
*/
{

   //generate random point within unit sphere
   GPoint vec(0.0, 0.0, 0.0);

   while(vec.isZero())
   {
      vec = GPoint(double(rand())/(RAND_MAX+1.0) - .5, 
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
void rayTraceRay(Shape& g_scene, vector<Light> lights, const Line& ray, GPoint& color, unsigned int 
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
   g_scene.intersection(ray, GPoint(0.0, 0.0, 0.0), intersection);
   
   if(!intersection.intersects()) return;
   
   GPoint pt = intersection.point();
   Material material = intersection.material();
   Line reflectedRay = intersection.reflectedRay();
   Line transmittedRay = intersection.transmittedRay();

   
   Line shadowRay;
   GPoint lColor;
   
   size_t size = lights.size();
   for(size_t i = 0 ; i < size; i++)
   {
      shadowRay.set(pt, lights[i].position());
      Intersection shadowIntersection;

      g_scene.intersection(shadowRay, GPoint(0.0, 0.0, 0.0), shadowIntersection );
      
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
      GPoint transmittedColor(0.0, 0.0, 0.0);
      GPoint reflectedColor(0.0, 0.0, 0.0);

      GPoint transparency = material.transparency();
      GPoint opacity = GPoint(1.0, 1.0, 1.0) - transparency;
            
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
void rayTraceScreen(Shape& g_scene, vector<Light>& lights, GPoint camera, GPoint lookAt, GPoint up, int
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
   GPoint lookDirection = lookAt - camera;
   GPoint right = lookDirection * up;
  
   right.normalize();
   GPoint rightOffset = width*right;

   up = right*lookDirection;
   up.normalize();
   
   GPoint screenPt = lookAt + bottomX*right + bottomY*up;
   
   Line ray;
   GPoint color(0.0, 0.0, 0.0);
   GPoint avgColor(0.0, 0.0, 0.0);
   
   GPoint weightedColor(0.0, 0.0, 0.0);
   GPoint oldWeightedColor(0.0, 0.0, 0.0);
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
            //if ( k <  SUPER_SAMPLE_NUMBER && k >1) cout <<"hello"<<k<<endl;      
    
            avgColor /= k;            
            glColor3d(avgColor.x(), avgColor.y(), avgColor.z());

            glVertex2i(i, j);
            screenPt += right;

         }
         screenPt -= rightOffset;
         screenPt += up;
      }
   glEnd();

}

/*-----------------------------------------------*/
GPoint  convertStringCoordinate(string coordString)
/*
PURPOSE: Converts a string of two characters of the form: 
   letter row + number color (for example, b4) into coordinates for a piece of 
   our ray tracing work. 
RECEIVES: coordString -- string to convert
RETURNS: Nothing
REMARKS:    
*/
{
      GPoint firstSquare(-BOARD_EDGE_SIZE/2, 0.0, BOARD_EDGE_SIZE/2);
            
      GPoint rowOffset(0.0, 0.0, -(double(coordString[0] - 'a') + .5)*SQUARE_EDGE_SIZE);
         //negative because farther back == higher row number
      GPoint colOffset((double(coordString[1] - '0' - 1) + .5)*SQUARE_EDGE_SIZE, 0.0, 0.0);
      GPoint heightOffset(0.0, 1.5*SQUARE_EDGE_SIZE, 0.0);
      
      GPoint square = firstSquare + rowOffset + colOffset + heightOffset;

      return square;
}

/*-----------------------------------------------*/
void MySdlApplication::initScene()
/*
PURPOSE: gets locations of objects from users
         sets the background color to black,
         sets up display list for coordinate axes,
                 initializes extrusions vector
                 
RECEIVES: Nothing
RETURNS: Nothing
REMARKS:    
*/
{

   GPoint boardPosition(0.0, 0.0, 0.0);
   CheckerBoard *checkerBoard = new CheckerBoard(boardPosition);
   g_scene.addRayObject(checkerBoard);

   string tmp;
   
   cout << "Please enter the position of the light:\n";
   cin >> tmp;
   //tmp ="b6"; //commented values are nice values to demonstrate the ray tracer
   g_lightPosition = GPoint(BOARD_POSITION) + GPoint(0.0, 3.5*SQUARE_EDGE_SIZE, 0.0) + convertStringCoordinate(tmp);
      //with what convertStringCoordinate gives makes 5 squares above board

   cout << "Please enter the position of the tetrahedron:\n";
   cin >>tmp;
   //tmp = "b4";
   Tetrahedron *tetrahedron = new Tetrahedron(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE);
   g_scene.addRayObject(tetrahedron);

   cout << "Please enter the position of the sphere:\n";
   cin >> tmp;
   //tmp = "d7";
   Sphere *sphere = new Sphere(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE/2);
   g_scene.addRayObject(sphere);

   cout << "Please enter the position of the cube:\n";
   cin >> tmp;
   //tmp = "a7";
   Cube *cube = new Cube(convertStringCoordinate(tmp), SQUARE_EDGE_SIZE);   
   g_scene.addRayObject(cube);
   
   

   glClearColor(0.0, 0.0, 0.0, 0.0);
        
   glViewport(0, 0, g_windowWidth, g_windowHeight);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, g_windowWidth, 0, g_windowHeight);


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
   
   GPoint camera(CAMERA_POSITION);
   GPoint lookAt(LOOK_AT_VECTOR);   
   GPoint up(UP_VECTOR); 
   
   rayTraceScreen(g_scene, lights, camera, lookAt, up, -g_windowWidth/2, -g_windowHeight/2, g_windowWidth, g_windowHeight);
 
   glFlush();
}

/*-----------------------------------------------*/
void reshapeFn(int newWidth, int newHeight)
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
	glUseProgram(g_shaderStates[g_activeShader]->program);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// clear framebuffer color&depth
	draw();

	SDL_GL_SwapWindow(display);
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
	while(running) 
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
	
	if(SDL_Init(SDL_INIT_EVERYTHING) < 0) 
	{
		return false;
   }
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	//SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);

	/* Turn on double buffering with a 24bit Z buffer.
	* You may need to change this to 16 or 32 for your system */
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	if((display = SDL_CreateWindow("Ray Trace",
	g_initX, g_initY, g_windowWidth, g_windowHeight,
		SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE)) == NULL) 
	{
		return false;
	}

	/* Create our opengl context and attach it to our window */
	SDL_GLContext maincontext = SDL_GL_CreateContext(display);
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

   initScene();
	initGLState();
	initShaders();
	//initGeometry();
	//initTextures();
	initCamera();

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
		running = false;
	else if (type == SDL_MOUSEBUTTONDOWN)
		mouse(event->button);
	else if (type == SDL_MOUSEBUTTONUP)
		mouse(event->button);
	else if (type == SDL_MOUSEMOTION)
		motion(event->motion.x, event->motion.y);
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

	running = true;
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

/****************************************************** 
* Copyright (c):   1994, All Rights Reserved. 
* Project:         CS 46A Homework #4 
* File:            sortcomp.cpp 
* Purpose:         compare timings for sort routines 
* Start date:      4/2/97 
* Programmer:      John Chen 
* 
****************************************************** 
*/


/*-----------------------------------------------*/



	/*	PURPOSE:		What does this function do? (must be present) 
		RECEIVES:	List every argument name and explain each argument. 
						(omit if the function has no arguments) 
		RETURNS:		Explain the value returned by the function. 
						(omit if the function returns no value) 
		REMARKS:		Explain any special preconditions or postconditions. 
						See example below. (omit if function is unremarkable) 
	*/