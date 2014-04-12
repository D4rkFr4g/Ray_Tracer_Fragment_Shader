#version 130

uniform vec3 uLight;
uniform sampler2D uTexUnit2;

in vec3 vNormal;
in vec3 vPosition;
in vec2 vTexCoord2;

out vec4 fragColor;

void main() {
  vec3 tolight = normalize(uLight - vPosition);
  vec3 normal = normalize(vNormal);
  vec4 texColor0 = texture2D(uTexUnit2, vTexCoord2);

  float diffuse = max(0.0, dot(normal, tolight));
  fragColor = texColor0;// * diffuse;
}
