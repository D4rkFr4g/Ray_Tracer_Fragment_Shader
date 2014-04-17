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

const vec3 WHITE = vec3(1.0);
const vec3 BLACK = vec3(0.0);
const vec3 RED = vec3(1.0, 0.0, 0.0);

in vec3 vPosition;
out vec4 fragColor;

// structure to help with psuedo recursion
struct RayBouncer
{
   int bounces;
   vec3 startRay;
   vec3 endRay;
   vec4 output;
};

struct Line
{
   vec3 startPt;
   vec3 endPt;
};

struct Material
{
   vec3 ambient;
   vec3 diffuse;
   vec3 specular;
   vec3 transparency;
   double refraction;
};

struct Intersection
{
   bool intersects;
   vec3 point;
   vec3 normal;

   Material material;
   Line reflectedRay;
   Line transmittedRay;
};

struct Shape
{
   float type;
   float radius;
   float edge;
   vec3 pos;
   float height;
   Material material;
};


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

void intersectionCircle(Shape circle, inout Line ray, inout vec3 positionOffset, inout Intersection inter)
{
   vec3 u = ray.endPt - ray.startPt;
   vec3 p0 = ray.startPt;
   vec3 position = circle.pos + positionOffset;
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

RayBouncer createRayBouncer()
{
   RayBouncer ray;
   ray.bounces = 0;
   ray.output = vec4(0,0,0,0);
   return ray;
}

RayBouncer ORD(RayBouncer ray)
{
   ray.bounces++;
   vec4 tmp;

	// this is dummy code, you could start implementing your ray tracer here
	/*
	    int i;
	    for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
		{
		   //check for collision of object with ray etc
		}
	 */
	//dummy code till above implemented
   float delta = 0.05 * (ray.bounces - 1);
	if(ray.endRay.x * ray.endRay.x + ray.endRay.y * ray.endRay.y < (0.2 - delta)) 
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
   ray.output += tmp;
   /*
   ray.startRay = ray.endRay;
   ray.endRay = Normal Bounce Ray
   */
   
   return ray;
}

vec4 ORD(RayBouncer ray, int bounces)
{
   for (int i = 0; i < bounces; i++)
      ray = ORD(ray);

   return ray.output;
}

void main() 
{
	vec3 eyePos = vec3(uEyePosition[0], uEyePosition[1], uEyePosition[2]);
	//this demonstrates a call to a function in a shader
   
   RayBouncer ray = createRayBouncer();
   ray.startRay = eyePos;
   ray.endRay = vPosition;

   /*
   fragColor = outgoingRadianceDepth1(eyePos, vPosition) 
	    + ATTENUATION * outgoingRadianceDepth2(eyePos, vPosition);
   */

   fragColor = ORD(ray, 3); // outgoingRadianceDepth
}