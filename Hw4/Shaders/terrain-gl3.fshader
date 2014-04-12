#version 130

uniform vec3 uLight;

out vec4 fragColor;

void main()
{
fragColor = gl_Color;
}