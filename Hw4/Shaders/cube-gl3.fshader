#version 130

uniform samplerCube uTexUnit2;
uniform vec3 uLight, uColor;

in vec3 vNormal;
in vec3 vPosition;

out vec4 fragColor;

vec3 reflect(vec3 w, vec3 n)
{
    return -w + n * (dot(w,n) * 2.0);
}

void main(void)
{
    vec3 tolight = normalize(uLight - vPosition);
    vec3 normal = normalize(vNormal);
    vec3 reflected = reflect(normalize(vec3(-vPosition)), normal);

    vec4 texColor0 = textureCube(uTexUnit2, reflected) / 3.0;

    float diffuse = max(0.0, dot(normal, tolight));
    vec3 intensity = (uColor * diffuse) + vec3(texColor0);

    fragColor = vec4(intensity, 1.0);
}