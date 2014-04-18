#version 400

int NUM_LIGHTS = 1; // should be set in C++ land
int NUM_SHAPES = 1; //should be set in C++ land
uniform float uEyePosition[3];
int LIGHT_STRIDE = 3;
int GEOMETRY_STRIDE = 6;

uniform float uLights[3]; // Should be NUM_LIGHTS * LIGHT_STRIDE
uniform float uGeometry[6]; // should be GEOMETRY_STRIDE*NUM_SHAPES

const double SMALL_NUMBER = .0001;
const int TETRAHEDRON = 0;
const int CUBE = 1;
const int SPHERE = 2;
const int CYLINDER = 3;
const int CONE = 4;
const double ATTENUATION_FACTOR = 100000; 

const vec3 WHITE = vec3(1.0);
const vec3 BLACK = vec3(0.0);
const vec3 RED = vec3(1.0, 0.0, 0.0);
const vec3 UP_VECTOR = vec3(0.0, 1.0, 0.0);
const vec3 LOOK_AT_VECTOR = vec3(0.0, 0.0, -160.0);

vec3 CAMERA = vec3(uEyePosition[0], uEyePosition[1], uEyePosition[2]);

in vec3 vPosition;
out vec4 fragColor;
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
   vec3 color;
   vec3 position;
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
   normalize(p);
   return p;
}

double attenuation(double distance)
{
   return ATTENUATION_FACTOR/(ATTENUATION_FACTOR + distance*distance);
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
      double discriminant = uDeltaP*uDeltaP - (deltaP * deltaP) + circle.radius*circle.radius; //dot(deltaP, deltaP)
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
void intersection(Shape shape, inout Line ray, inout Intersection inter)
{
   if (shape.type == SPHERE)
      intersectionCircle(shape, ray, inter);
   else  // Remove once rest are finished
      intersectionCircle(shape, ray, inter);
}
/*-----------------------------------------------*/
void createShape(inout Shape shape, int index)
{
   float type;
   float radius;
   float edge;
   vec3 pos;
   float height;
   Material material;

   if (uGeometry[index] == SPHERE)
   {
      type = SPHERE;
      radius = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      height = 0;
      material = sphereMaterial;
   }
   else // Remove later after all are defined
   {
      type = SPHERE;
      radius = uGeometry[++index];
      pos = vec3(uGeometry[++index], uGeometry[++index], uGeometry[++index]);
      height = 0;
      material = sphereMaterial;
   }

   shape = Shape(type, radius, edge, pos, height, material);
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
   if (ray.isBouncing)
   {
      Line rayLine = Line(ray.startPt, ray.endPt);
      ray.bounces++;
      vec4 tmp;   

      double minDistance = -1.0;
      double distanceTmp;
      Intersection inter;

      // Find closest intersection point
      int i = 0;
      for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
      {
         Intersection interTmp;
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
         return;

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

         Light light = Light(vec3(uLights[i], uLights[i+1], uLights[i+2]), WHITE);
         shadowRay = Line(pt, light.position);
      
         for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
         {
            Intersection shadowInterTmp;
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
      /*
      if(!isZero(opacity))
      {
         // Would need a new function unfortunately to mimick recursion
      }
      */
	
      //////////////////////////////////////////////////////////////////
	   //dummy code till above implemented
      float delta = 0.05 * (ray.bounces - 1);
	   if(ray.endPt.x * ray.endPt.x + ray.endPt.y * ray.endPt.y < (0.2 - delta)) 
      {
         float colorDelta = 0.2 * (ray.bounces - 1);
         vec4 color = vec4(colorDelta, colorDelta, colorDelta, colorDelta); 
         tmp = vec4(0.5, 0.5, 0.5, 0.5) + color;
      }
	   else
		   tmp = vec4(0.0, 0.0, 0.0, 0.0); 
	
      if (ray.bounces > 1)
         tmp *= ATTENUATION;

      // Update ray fields
      ray.color += tmp;
      /*
      ray.startRay = ray.endRay;
      ray.endRay = Normal Bounce Ray
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
	//this demonstrates a call to a function in a shader
   
   RayBouncer ray;
   createRayBouncer(ray);
   ray.startPt = eyePos;
   ray.endPt = vPosition;

   /*
   fragColor = outgoingRadianceDepth1(eyePos, vPosition) 
	    + ATTENUATION * outgoingRadianceDepth2(eyePos, vPosition);
   */
   
   fragColor = ORD(ray,3); // outgoingRadianceDepth
}
/*-----------------------------------------------*/