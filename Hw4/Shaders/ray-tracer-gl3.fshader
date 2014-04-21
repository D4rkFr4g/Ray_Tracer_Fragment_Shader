#version 400

int NUM_LIGHTS = 1; // should be set in C++ land
int NUM_SHAPES = 1; //should be set in C++ land
uniform float uEyePosition[3];
int LIGHT_STRIDE = 3;
int GEOMETRY_STRIDE = 6;

uniform float uLights[3]; // Should be NUM_LIGHTS * LIGHT_STRIDE
uniform float uGeometry[6]; // should be GEOMETRY_STRIDE*NUM_SHAPES
uniform mat4 uInverseModelViewMatrix;
uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;

const double SMALL_NUMBER = .0001;
const int TETRAHEDRON = 0;
const int CUBE = 1;
const int SPHERE = 2;
const int CYLINDER = 3;
const int CONE = 4;
const int CHECKERBOARD = 5;
const double ATTENUATION_FACTOR = 100000; 

const vec3 WHITE = vec3(1.0);
const vec3 BLACK = vec3(0.0);
const vec3 RED = vec3(1.0, 0.0, 0.0);
const vec4 dRED = vec4(1, 0, 0,1);
const vec4 dGREEN = vec4(0, 1, 0,1);
const vec4 dBLUE = vec4(0, 0, 1,1);
const vec3 UP_VECTOR = vec3(0.0, 1.0, 0.0);
const vec3 LOOK_AT_VECTOR = vec3(0.0, 0.0, -160.0);

const unsigned int NUM_SQUARES = 8; //number of squares wide our chess board is

in vec3 gPosition;
in vec3 vPosition;
out vec4 fragColor;

//vec3 CAMERA = vec3(uEyePosition[0], uEyePosition[1], uEyePosition[2]);
vec4 debugColor = vec4(0,0,0,1);

/*-----------------------------------------------*/
// structure to help with psuedo recursion
struct RayBouncer
{
   int bounces;
   vec3 startPt;
   vec3 endPt;
   vec4 color;
   vec3 transmittedColor;
   vec3 reflectedColor;
   vec3 transparency;
   bool isBouncing;
};
/*-----------------------------------------------*/
struct Light
{
   vec3 position;
   vec3 color;
};

struct Line
{
   vec3 startPt;
   vec3 endPt;
};
/*-----------------------------------------------*/
struct Material
{
   vec3 ambient;
   vec3 diffuse;
   vec3 specular;
   vec3 transparency;
   double refraction;
};
/*-----------------------------------------------*/
struct Intersection
{
   bool intersects;
   vec3 point;
   vec3 normal;

   Material material;
   Line reflectedRay;
   Line transmittedRay;
};
/*-----------------------------------------------*/
struct Shape
{
   float type;
   float radius;
   float edge;
   vec3 pos;
   float height;
   Material material;
   int squares;
   double board_half_size;
   double square_edge_size;
};
/*-----------------------------------------------*/
struct Triangle
{
   vec3 v0;
   vec3 v1;
   vec3 v2;
   vec3 n;
};
/*-----------------------------------------------*/
struct Side
{
   vec3 v0;
   vec3 v1;
   vec3 v2;
   vec3 v3;
};
/*-----------------------------------------------*/

Material whiteSquare = Material(.1*WHITE, .5*WHITE, WHITE, BLACK, 1);
Material blackSquare = Material(BLACK, .1*WHITE, BLACK, BLACK, 1);
Material sphereMaterial = Material(BLACK, .1*WHITE, WHITE, BLACK, 1);
Material tetrahedronMaterial = Material(BLACK, BLACK, .1*WHITE, WHITE, 2.0/3.0);
Material cubeMaterial = Material(.1*RED, .4*RED, RED, BLACK, 1);

/*
   Format in array for shapes:
   If type is:
   TETRAHEDRON: type, edge_length, c_x, c_y, c_z, dummy
   CUBE: type, edge_length, c_x, c_y, c_z, dummy
   SPHERE: type, radius, c_x, c_y, c_z, dummy
   CYLINDER: type, radius, c_x, c_y, c_z, height
   CONE: type, radius, c_x, c_y, c_z, height
   CHECKERBOARD: type, edge_length, c_x, c_y, c_z, squares
 */
float ATTENUATION = 0.3;
/*-----------------------------------------------*/

double length(vec3 v)
{
   return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

/*-----------------------------------------------*/
double length(Line l)
{
   vec3 p = l.endPt - l.startPt;
   return length(p);
}
/*-----------------------------------------------*/
bool isZero(vec3 v)
{
   return v.x == 0 && v.y == 0 && v.z == 0;
}
/*-----------------------------------------------*/
vec3 direction(Line line)
{
   vec3 p = line.endPt - line.startPt;
   p = normalize(p);
   return p;
}

double attenuation(double distance)
{
   return ATTENUATION_FACTOR/(ATTENUATION_FACTOR + distance*distance);
}
/*-----------------------------------------------*/
void distanceCompare(inout double minDistance, inout Intersection inter, Intersection interTmp, Line ray) 
{
   vec3 directionCur = inter.point - ray.startPt;
   double distanceTmp = length(directionCur);

   if (distanceTmp < minDistance || minDistance < 0.0)
   {
      minDistance = distanceTmp;
      inter = interTmp;
   }
}
/*-----------------------------------------------*/
void createIntersection(inout Intersection inter)
{
   inter.intersects = false;
   inter.point = vec3(0,0,0);
   inter.normal = vec3(0,0,0);
   inter.material = blackSquare;
   inter.reflectedRay = Line(vec3(0),vec3(0));
   inter.transmittedRay = Line(vec3(0),vec3(0));
}
/*-----------------------------------------------*/
void intersectionCircle(Shape circle, inout Line ray, inout Intersection inter)
{
   vec3 u = ray.endPt - ray.startPt;
   vec3 p0 = ray.startPt;
   vec3 position = circle.pos;
   vec3 deltaP = position - p0;

   if (circle.radius > 0)
   {
      double uDeltaP = dot(u,deltaP);
      double discriminant = uDeltaP*uDeltaP - dot(deltaP, deltaP) + circle.radius*circle.radius; //dot(deltaP, deltaP)
      double s = uDeltaP - sqrt(discriminant);

      if (discriminant < 0 || abs(s) < SMALL_NUMBER)
      {
         inter.intersects = false;
         return;
      }

      //calculate point of intersection
      vec3 su = vec3(s*u.x, s*u.y, s*u.z);
      vec3 p = p0 + su;
      vec3 directionP0 = p - position;

      if(s < SMALL_NUMBER) // if not in front of ray then don't intersect
      {
         inter.intersects = false;
         return;
      }

      //reflected vector calculated using equations from book
      vec3 n = directionP0;
      n = normalize(n);

      vec3 r = reflect(u, n);
      Line reflected = Line(p, p + r);
      
      
      //Transmitted vector calculated using thin lens equations from book
      vec3 t = vec3(0.0);
      double refractionRatio = circle.material.refraction;
      double cosThetai = dot(u, n);
      double modulus = 1 - refractionRatio*refractionRatio*( 1- cosThetai*cosThetai);

      if( modulus > 0)
      {
         double cosThetar = sqrt(modulus);
         double cT = (cosThetar + refractionRatio * cosThetai);
         vec3 nC = vec3(cT*n.x, cT*n.y, cT*n.z);
         vec3 rU = vec3(refractionRatio*u.x, refractionRatio*u.y, refractionRatio*u.z);
         t = rU - nC;
      }
      Line transmitted = Line(p, p + t);

      inter = Intersection(true, p, n, circle.material, reflected, transmitted);
   }
}
/*-----------------------------------------------*/
void intersectionTriangle(Triangle triangle, inout Line ray, inout Intersection inter)
{
   //get coordinates of triangle given our position
   //vec3 position = _position + positionOffset;
   vec3 v0 = triangle.v0;
   vec3 v1 = triangle.v1;
   vec3 v2 = triangle.v2;
   
   vec3 _u = v1 - v0;
   vec3 v = v2 - v0;
   vec3 n = cross(_u,v);

   //handle last degenerates case by saying we don't intersect
   if(length(n) < SMALL_NUMBER) 
      return;

   n = normalize(n);

   double uv = dot(_u,v);
   double uu = dot(_u,_u);
   double vv = dot(v,v);
   
   double denominator = (uv*uv) - (uu*vv);

   if( abs(denominator) < SMALL_NUMBER)
         return;

   vec3 p0 = ray.startPt;
   vec3 p1 = ray.endPt;
   vec3 diffP = p1 - p0;
   double ndiffP = dot(n, diffP);

   //handle another degenerate case by saying we don't intersect
   if( abs(ndiffP) < SMALL_NUMBER)
   {
      inter.intersects = false;
      return;
   }
   
   double m = dot(n, (v0 - p0))/ dot(n, diffP);
   if( m < SMALL_NUMBER) //if m is negative then we don't intersect
   {
      inter.intersects = false;
      return;
   }

   vec3 p = p0 + vec3(diffP.x*m, diffP.y*m, diffP.z*m); //intersection point with plane

   vec3 w = p - v0;

   //Now check if in triangle   
    
   double wu = dot(w, _u);
   double wv = dot(w, v);

   double s = (uv*wv - vv*wu)/denominator;
   double t = (uv*wu - uu*wv)/denominator;

   if( s >= 0 && t >=0 && s + t <= 1) // intersect
   {
      diffP = normalize(diffP);// now u is as in the book
      vec3 u = diffP;

      vec3 r = reflect(u,n);
      Line reflected = Line(p, p + r);

      /*
      //Transmitted vector calculated using thin lens equations from book
      double refractionRatio = material.refraction();

      Point t(0.0, 0.0, 0.0);

      GLdouble cosThetai = u & _n;
      GLdouble modulus = 1 - refractionRatio*refractionRatio*( 1- cosThetai*cosThetai);

      if( modulus > 0)
      {
         GLdouble cosThetar = sqrt(modulus);
         t = refractionRatio*u - ( cosThetar + refractionRatio*cosThetai)*_n;
      }
      */

      Line transmitted = Line(vec3(0), vec3(0));
      Material mat = whiteSquare;
      
      inter = Intersection(true, p, n, mat, reflected, transmitted);
   }
   else // don't intersect
   {
      inter.intersects = false;
   }
}
/*-----------------------------------------------*/
void intersectionSquare(Shape square, inout Line ray, inout Intersection inter)
{
   double halfSize = square.board_half_size;
   vec3 p = square.pos;
   vec3 v0 = p + vec3(-halfSize, 0, -halfSize);
   vec3 v1 = p + vec3( halfSize, 0, -halfSize);
   vec3 v2 = p + vec3( halfSize, 0,  halfSize);
   vec3 v3 = p + vec3(-halfSize, 0,  halfSize);
   
   // Convert wCoords to camCoords
   /*
   v0 = vec3(uModelViewMatrix * vec4(v0,1));
   v1 = vec3(uModelViewMatrix * vec4(v1,1));
   v2 = vec3(uModelViewMatrix * vec4(v2,1));
   v3 = vec3(uModelViewMatrix * vec4(v3,1));
   */

   Triangle t0 = Triangle(v0, v1, v2, UP_VECTOR);
   Triangle t1 = Triangle(v0, v2, v3, UP_VECTOR);

   Intersection inter0;
   Intersection inter1;
   createIntersection(inter0);
   createIntersection(inter1);
   
   intersectionTriangle(t0, ray, inter0);
   intersectionTriangle(t1, ray, inter1);

   if (inter0.intersects)
   {
      inter = inter0;
   }
   else if (inter1.intersects)
   {
      inter = inter1;
   }
   else
   {
      inter.intersects = false;
   }
}
/*-----------------------------------------------*/
void intersectionCheckerBoard(Shape board, inout Line ray, inout Intersection inter)
{   
   
   intersectionSquare(board, ray, inter);

   if(inter.intersects)
   {
      double edge = 5;
      //vec3 boardOffset = vec3(uModelViewMatrix * vec4(board.board_half_size,0,board.board_half_size,1));
      //vec3 p = inter.point - board.pos + boardOffset;
      
      vec3 p = inter.point - board.pos + vec3(board.board_half_size, 0, board.board_half_size);
      
      //int squareSum = int(p.x/(1/board.square_edge_size)) + int(p.z/(1/board.square_edge_size));
      int squareSum = int(p.x/(board.square_edge_size)) + int(p.z/(board.square_edge_size));
      //int squareSum = int(p.x/edge) + int(p.z/edge);

      if((squareSum & 1) == 0)
         inter.material = whiteSquare;
      else
         inter.material = blackSquare;
   }

}
/*-----------------------------------------------*/
void intersectionTetrahedron(Shape tetrahedron, inout Line ray, inout Intersection inter)
{
   double halfEdge = tetrahedron.edge / 2;
   vec3 pos = tetrahedron.pos;
   vec3 n = vec3(0,0,0);

   Triangle t0 = Triangle(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, -halfEdge, halfEdge), n);
   Triangle t1 = Triangle(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, -halfEdge, halfEdge), 
      pos + vec3(-halfEdge, halfEdge, -halfEdge), n);
   Triangle t2 = Triangle(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, -halfEdge, halfEdge), n);
   Triangle t3 = Triangle(pos + vec3(-halfEdge, -halfEdge, halfEdge), 
      pos + vec3(halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, halfEdge, -halfEdge), n);

   Intersection interT0;
   Intersection interT1;
   Intersection interT2;
   Intersection interT3;
   createIntersection(interT0);
   createIntersection(interT1);
   createIntersection(interT2);
   createIntersection(interT3);
   
   intersectionTriangle(t0, ray, interT0);
   intersectionTriangle(t1, ray, interT1);
   intersectionTriangle(t2, ray, interT2);
   intersectionTriangle(t3, ray, interT3);
   
   // When multiple intersections check distance for closest
   double minDistance = -1;
   double distanceTmp = 0;

   if (interT0.intersects)
      distanceCompare(minDistance, inter, interT0, ray);
   if (interT1.intersects)
      distanceCompare(minDistance, inter, interT1, ray);
   if (interT2.intersects)
      distanceCompare(minDistance, inter, interT2, ray);
   if (interT3.intersects)
      distanceCompare(minDistance, inter, interT3, ray);

   inter.material = tetrahedron.material;
}
/*-----------------------------------------------*/
void intersectionPlane(vec3 v0, vec3 n, Line ray, inout Intersection inter)
{
   vec3 p0 = ray.startPt;
   vec3 p1 = ray.endPt;
   vec3 diffP = p1 - p0;
   double ndiffP = dot(n, diffP);
   
   //handle another degenerate case by saying we don't intersect
   if( abs(ndiffP) < SMALL_NUMBER)
   {
      inter.intersects = false;
      return;
   } 

   double m = dot(n, (v0 - p0))/ dot(n, diffP);

   if( m < SMALL_NUMBER) //if m is negative thenwe don't intersect
   {
      inter.intersects = false;
      return;
   }
   vec3 p = p0 + vec3(diffP.x*m, diffP.y*m, diffP.z*m); //intersection point with plane
}
/*-----------------------------------------------*/
void inSide(inout Intersection inter, Side side)
{
   vec3 p = inter.point;

   //Find min maxs
   double minX, maxX;
   double minY, maxY;
   double minZ, maxZ;

   minX = min(side.v0.x, min(side.v1.x, min(side.v2.x, side.v3.x)));
   maxX = max(side.v0.x, max(side.v1.x, max(side.v2.x, side.v3.x)));
   minY = min(side.v0.y, min(side.v1.y, min(side.v2.y, side.v3.y)));
   maxY = max(side.v0.y, max(side.v1.y, max(side.v2.y, side.v3.y)));
   minZ = min(side.v0.z, min(side.v1.z, min(side.v2.z, side.v3.z)));
   maxZ = max(side.v0.z, max(side.v1.z, max(side.v2.z, side.v3.z)));

   if (p.x < minX || p.x > maxX)
      inter.intersects = false;
   if (p.y < minY || p.y > maxY)
      inter.intersects = false;
   if (p.z < minZ || p.z > maxZ)
      inter.intersects = false;
}
/*-----------------------------------------------*/
void intersectionCube(Shape cube, inout Line ray, inout Intersection inter)
{
   double halfEdge = cube.edge / 2;
   vec3 pos = cube.pos;
   vec3 n = vec3(0,0,0);
   vec3 p = vec3(0,0,0);

   Side xP = Side(pos + vec3(halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, halfEdge, -halfEdge), 
      pos + vec3(halfEdge, halfEdge, halfEdge), 
      pos + vec3(halfEdge, -halfEdge, halfEdge));
   Side xN = Side(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, halfEdge, halfEdge), 
      pos + vec3(-halfEdge, -halfEdge, halfEdge));
   Side yP = Side(pos + vec3(-halfEdge, halfEdge, -halfEdge), 
      pos + vec3(halfEdge, halfEdge, -halfEdge), 
      pos + vec3(halfEdge, halfEdge, halfEdge), 
      pos + vec3(-halfEdge, halfEdge, halfEdge));
   Side yN = Side(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, -halfEdge, halfEdge), 
      pos + vec3(-halfEdge, -halfEdge, halfEdge));
   Side zP = Side(pos + vec3(-halfEdge, -halfEdge, halfEdge), 
      pos + vec3(halfEdge, -halfEdge, halfEdge), 
      pos + vec3(halfEdge, halfEdge, halfEdge), 
      pos + vec3(-halfEdge, halfEdge, halfEdge));
   Side zN = Side(pos + vec3(-halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, -halfEdge, -halfEdge), 
      pos + vec3(halfEdge, halfEdge, -halfEdge), 
      pos + vec3(-halfEdge, halfEdge, -halfEdge));

   // +z Plane
   Intersection interZp;
   createIntersection(interZp);
   n = vec3(0,0,1);
   p = pos + vec3(-halfEdge, -halfEdge, halfEdge);
   intersectionPlane(p, n, ray, interZp);
   if (interZp.intersects)
      inSide(interZp, zP);

   // -z Plane
   Intersection interZn;
   createIntersection(interZn);
   n = vec3(0,0,-1);
   p = pos + vec3(-halfEdge, -halfEdge, -halfEdge);
   intersectionPlane(p, n, ray, interZn);
   if (interZn.intersects)
      inSide(interZn, zN);

   // +x Plane
   Intersection interXp;
   createIntersection(interXp);
   n = vec3(1,0,0);
   p = pos + vec3(halfEdge, -halfEdge, -halfEdge);
   intersectionPlane(p, n, ray, interXp);
   if (interXp.intersects)
      inSide(interXp, xP);

   // -x Plane
   Intersection interXn;
   createIntersection(interXn);
   n = vec3(-1,0,0);
   p = pos + vec3(-halfEdge, -halfEdge, -halfEdge);
   intersectionPlane(p, n, ray, interXn);
   if (interXn.intersects)
      inSide(interXn, xN);

   // +y Plane
   Intersection interYp;
   createIntersection(interYp);
   n = vec3(0,1,0);
   p = pos + vec3(-halfEdge, halfEdge, -halfEdge);
   intersectionPlane(p, n, ray, interYp);
   if (interYp.intersects)
      inSide(interYp, yP);

   // -y Plane
   Intersection interYn;
   createIntersection(interYn);
   n = vec3(0,-1,0);
   p = pos + vec3(-halfEdge, -halfEdge, -halfEdge);
   intersectionPlane(p, n, ray, interYn);
   if (interYn.intersects)
      inSide(interYn, yN);

   double minDistance = -1;
   double distanceTmp = 0;
   
   // When multiple intersections check distance for closest
   if (interZp.intersects)
      distanceCompare(minDistance, inter, interZp, ray);
   if (interZn.intersects)
      distanceCompare(minDistance, inter, interZn, ray);
   if (interXp.intersects)
      distanceCompare(minDistance, inter, interXp, ray);
   if (interXn.intersects)
      distanceCompare(minDistance, inter, interXn, ray);
   if (interYp.intersects)
      distanceCompare(minDistance, inter, interYp, ray);
   if (interYn.intersects)
      distanceCompare(minDistance, inter, interYn, ray);

   inter.material = cube.material;
}
/*-----------------------------------------------*/
void intersectionCylinder(Shape cylinder, inout Line ray, inout Intersection inter)
{
   Intersection interTop;
   Intersection interBottom;
   Intersection interSide;
   createIntersection(interTop);
   createIntersection(interBottom);
   createIntersection(interSide);

   vec3 p = cylinder.pos;

   // Top
   double delta = (ray.startPt.z - p.z) / ray.endPt.z;
   vec3 iP = vec3(ray.startPt.x + (delta * ray.endPt.x),ray.startPt.y + (delta * ray.endPt.y), p.z);

   if ((iP.x - p.x)*(iP.x - p.x) + (iP.y - p.y)*(iP.y - p.y) <= cylinder.radius * cylinder.radius)
   {
      interTop.intersects = true;
      interTop.point = iP;
      interTop.normal = vec3(0,1,0);
      interTop.material = tetrahedronMaterial;
   }

   // Bottom

   // Side
}
/*-----------------------------------------------*/
void intersectionCone(Shape cone, inout Line ray, inout Intersection inter)
{

}
/*-----------------------------------------------*/
void intersection(Shape shape, inout Line ray, inout Intersection inter)
{
   if (shape.type == SPHERE)
      intersectionCircle(shape, ray, inter);
   else if (shape.type == CHECKERBOARD)
      intersectionCheckerBoard(shape, ray, inter);
   else if (shape.type == TETRAHEDRON)
      intersectionTetrahedron(shape, ray, inter);
   else if (shape.type == CUBE)
      intersectionCube(shape, ray, inter);
   else if (shape.type == CYLINDER)
      intersectionCylinder(shape, ray, inter);
   else if (shape.type == CONE)
      intersectionCone(shape, ray, inter);
}
/*-----------------------------------------------*/
void createShape(inout Shape shape, int index)
{
   float type;
   float radius = 0;
   float edge = 0;
   vec3 pos;
   float height = 0;
   Material material;
   int squares = 0;
   double board_half_size = 0;
   double square_edge_size = 0;

   if (uGeometry[index] == SPHERE)
   {
      type = SPHERE;
      radius = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      height = 0;
      material = sphereMaterial;
   }
   else if (uGeometry[index] == CHECKERBOARD)
   {
      type = CHECKERBOARD;
      edge = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      squares = int(uGeometry[++index]);
      board_half_size = edge / 2;
      square_edge_size = edge / squares;
   }
   else if (uGeometry[index] == TETRAHEDRON)
   {
      type = TETRAHEDRON;
      edge = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      material = tetrahedronMaterial;
   }
   else if (uGeometry[index] == CUBE)
   {
      type = CUBE;
      edge = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      material = cubeMaterial;
   }
   else if (uGeometry[index] == CYLINDER)
   {
      type = CYLINDER;
      radius = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      height = uGeometry[++index];
      material = tetrahedronMaterial;
   }
   else if (uGeometry[index] == CONE)
   {
      type = CONE;
      radius = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      height = uGeometry[++index];
      material = tetrahedronMaterial;
   }

   shape = Shape(type, radius, edge, pos, height, material, squares, board_half_size, square_edge_size);
}
/*-----------------------------------------------*/
void createRayBouncer(inout RayBouncer ray)
{
   ray.bounces = 0;
   ray.color = vec4(0,0,0,0);
   ray.transmittedColor = vec3(0.0, 0.0, 0.0);
   ray.reflectedColor = vec3(0.0, 0.0, 0.0);
   ray.transparency = vec3(1.0, 1.0, 1.0);
   ray.isBouncing = true;
}
/*-----------------------------------------------*/
void ORD(inout RayBouncer ray)
{
   vec4 tmp; 

   if (ray.isBouncing)
   {
      Line rayLine = Line(ray.startPt, ray.endPt);
      ray.bounces++;

      double minDistance = -1.0;
      double distanceTmp;
      Intersection inter;
      createIntersection(inter);

      // Find closest intersection point
      //int i = 0;
      for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
      {
         Intersection interTmp;
         createIntersection(interTmp);

         Shape shape;
         createShape(shape, i);

         intersection(shape, rayLine, interTmp);

         if (interTmp.intersects)
         {
            vec3 directionCur = interTmp.point - ray.startPt;
            distanceTmp = length(directionCur);

            if (distanceTmp < minDistance || minDistance < 0.0)
            {
               minDistance = distanceTmp;
               inter = interTmp;
            }
         }
      }

      if (!inter.intersects)
      {
         //debugColor = vec4(0,0,1,1);
         return;
      }

      vec3 pt = inter.point;
      Material material = inter.material;
      Line reflectedRay = inter.reflectedRay;
      Line transmittedRay = inter.transmittedRay;

      Line shadowRay;
      vec3 lColor;

      //TODO Shadow intersection
      // Doesn't really matter if it hits the light it matters if it hits something else
      for (int i = 0; i < NUM_LIGHTS * LIGHT_STRIDE; i += LIGHT_STRIDE)
      {
         double minDistance = -1.0;
         double distanceTmp;
         Intersection shadowIntersection;
         createIntersection(shadowIntersection);

         Light light = Light(vec3(uLights[i], uLights[i+1], uLights[i+2]), WHITE);
         
         shadowRay = Line(pt, light.position);
      
         for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
         {
            Intersection shadowInterTmp;
            createIntersection(shadowInterTmp);

            Shape shape;
            createShape(shape, i);

            intersection(shape, shadowRay, shadowInterTmp);
            if (shadowInterTmp.intersects)
            {
               vec3 directionCur = shadowInterTmp.point - ray.startPt;
               distanceTmp = length(directionCur);

               if (distanceTmp < minDistance || minDistance < 0.0)
               {
                  minDistance = distanceTmp;
                  shadowIntersection = shadowInterTmp;
               }
            }
         }

         if (!shadowIntersection.intersects || !isZero(shadowIntersection.material.transparency))
         {
            double att = attenuation(length(shadowRay));
            lColor = vec3(light.color.x*att, light.color.y*att, light.color.z*att);
         
            vec3 one = material.ambient * lColor;
            float delta = abs(dot(inter.normal, direction(shadowRay)));
            vec3 diffuse = material.diffuse * lColor;
            vec3 two = vec3(diffuse.x*delta, diffuse.y*delta, diffuse.z*delta);
            delta = abs(dot(direction(rayLine), direction(reflectedRay)));
            vec3 specular = (material.specular * lColor);
            vec3 three = vec3(specular.x*delta, specular.y*delta, specular.z*delta);
            vec3 color = one + two + three;
            color *= ray.transparency;
            ray.color += vec4(color.x, color.y, color.z, 1);
         }
      }

      // Setup for next bounce
      vec3 transparency = material.transparency;
      vec3 opacity = vec3(1.0, 1.0, 1.0) - transparency;

      // If not transparent then don't send ray
      if (!isZero(transparency) && length(transparency) > SMALL_NUMBER)
      {
         ray.startPt = transmittedRay.startPt;
         ray.endPt = transmittedRay.endPt;
         ray.transparency = material.transparency;
      }
      else
         ray.isBouncing = false;
      /*
      if(!isZero(opacity))
      {
         // Would need a new function unfortunately to mimick recursion
      }
      */
      /*
      //////////////////////////////////////////////////////////////////
	   //dummy code till above implemented
      if (!ray.isBouncing)
      {
         float delta = 0.05 * ray.bounces;
	      if(ray.endPt.x * ray.endPt.x + ray.endPt.y * ray.endPt.y < (0.2 - delta)) 
         {
            float colorDelta = 0.2 * (ray.bounces);
            vec4 color = vec4(colorDelta, colorDelta, colorDelta, colorDelta); 
            tmp = vec4(0.5, 0.5, 0.5, 0.5) + color;
         }
	      else
		      tmp = vec4(0.0, 0.0, 0.0, 0.0); 
	
         if (ray.bounces > 1)
            tmp *= ATTENUATION;

         // Update ray fields
         ray.color += tmp;
      }
      else
         ray.color = vec4(0,1,0,1);
      */
   }
}
/*-----------------------------------------------*/
vec4 ORD(inout RayBouncer ray, int bounces)
{
   for (int i = 0; i < bounces; i++)
      ORD(ray);

   return ray.color;
}
/*-----------------------------------------------*/
void main() 
{
   vec3 eyePos = vec3(uEyePosition[0], uEyePosition[1], uEyePosition[2]);

   RayBouncer ray;
   createRayBouncer(ray);
   ray.startPt = eyePos;
   //ray.endPt = gPosition;
   ray.endPt = vPosition;

   /*
   fragColor = outgoingRadianceDepth1(eyePos, vPosition) 
	    + ATTENUATION * outgoingRadianceDepth2(eyePos, vPosition);
   */
   
   fragColor = ORD(ray,1); // outgoingRadianceDepth

   if (debugColor != vec4(0,0,0,1))
      fragColor = debugColor;
}
/*-----------------------------------------------*/