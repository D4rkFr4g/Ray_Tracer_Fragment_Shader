#version 130

uniform float uEyePosition[3];
int NUM_LIGHTS = 1; // should be set in C++ land
int LIGHT_STRIDE = 3;
uniform float uLights[3]; // should be LIGHT_STRIDE*NUM_LIGHTS
int NUM_SHAPES = 1; //should be set in C++ land
int GEOMETRY_STRIDE = 6;
uniform float uGeometry[6]; // should be GEOMETRY_STRIDE*NUM_SHAPES

in vec3 vPosition;
out vec4 fragColor;

/*
   Format in array for shapes:
   If type is:
   TETRAHEDRON: type, edge_length, c_x, c_y, c_z, dummy
   CUBE: type, edge_length, c_x, c_y, c_z, dummy
   SPHERE: type, radius, c_x, c_y, c_z, dummy
   CYLINDER: type, radius, height, c_x, c_y, c_z
   CONE: type, radius, height, c_x, c_y, c_z 
 */
int TETRAHEDRON = 0;
int CUBE = 1;
int SPHERE = 2;
int CYLINDER = 3;
int CONE = 4;

float ATTENUATION = 0.3;

// Shaders don't support recursion, but we one need to do bounces to some fixed depth so fake

/*
   This function might handle one bounce reflections:
 */
vec4 outgoingRadianceDepth1(vec3 startRay, vec3 endRay)
{
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
	if(endRay.x*endRay.x + endRay.y*endRay.y < 0.2) {
		tmp = vec4(0.5, 0.5, 0.5, 0.5);

	} else {
		tmp = vec4(0.0, 0.0, 0.0, 0.0); 
	}
	return tmp;
}

/*
   This function might handle two bounce reflections:
 */
vec4 outgoingRadianceDepth2(vec3 startRay, vec3 endRay)
{
    vec4 tmp;
	// this is dummy code, you could start implementing your ray tracer here
	// this might involve a call to outgoingRadianceDepth1
	/*
	    int i;
	    for(int i = 0; i < NUM_SHAPES * GEOMETRY_STRIDE; i += GEOMETRY_STRIDE)
		{
		   //check for collision of object with ray etc
		}
	 */
	//dummy code till above implemented
	if(endRay.x*endRay.x + endRay.y*endRay.y < 0.1 ) {
		tmp = vec4(0.7, 0.7, 0.7, 0.7);
	} else {
		tmp = vec4(0.0, 0.0, 0.0, 0.0);
	}
	return tmp;
}

void main() 
{
	vec3 eyePos = vec3(uEyePosition[0], uEyePosition[1], uEyePosition[2]);
	//this demonstrates a call to a function in a shader
	fragColor = outgoingRadianceDepth1(eyePos, vPosition) 
	    + ATTENUATION * outgoingRadianceDepth2(eyePos, vPosition);
}