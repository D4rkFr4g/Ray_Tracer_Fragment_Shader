#version 130

uniform vec3 uLight;

varying vec4 fragColor;

void main()
{
	fragColor = gl_Color;
}