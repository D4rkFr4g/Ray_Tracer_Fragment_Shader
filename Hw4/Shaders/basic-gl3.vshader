#version 130

uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;

in vec3 aPosition;

out vec3 vPosition;
out vec3 gPosition;

void main() 
{
    // send position (eye coordinates) to fragment shader
    vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);
    vPosition = vec3(tPosition);
    gl_Position = uProjMatrix * tPosition;
	gPosition = aPosition;
}