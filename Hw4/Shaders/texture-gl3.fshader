#version 130

uniform vec3 uLight;
uniform sampler2D uTexUnit0, uTexUnit1;

in vec3 vNormal;
in vec3 vPosition;
in vec2 vTexCoord0;
in vec2 vTexCoord1;

out vec4 fragColor;

void main() {
  vec3 tolight = normalize(uLight - vPosition);
  vec3 normal = normalize(vNormal);
  vec4 texColor0 = texture2D(uTexUnit0, vTexCoord0);

  float diffuse = max(0.0, dot(normal, tolight));
  fragColor = texColor0;// * diffuse;
}
