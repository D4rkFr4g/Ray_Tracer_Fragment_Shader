#version 130

uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

void main()
{
	vec4 tPosition = uModelViewMatrix * vec4(gl_Position, 1.0);
	gl_Position = uProjMatrix * tPosition;
}