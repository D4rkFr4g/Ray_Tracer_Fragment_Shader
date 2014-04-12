#version 330

uniform vec3 uLight, uLight2, uColor;

in vec3 vNormal;
in vec3 vPosition;

out vec4 fragColor;

void main() 
{
    vec3 toLight = uLight - vec3(vPosition);
    vec3 toV = -normalize(vec3(vPosition));
    toLight = normalize(toLight);
    vec3 h = normalize(toV + toLight);
    vec3 normal = normalize(vNormal);

    float specular = pow(max(0.0, dot(h, normal)), 64.0);
    float diffuse = max(0.0, dot(normal, toLight));
    vec3 intensity = vec3(0.001, 0.001, 0.001) + uColor * diffuse
        + vec3(0.6, 0.6, 0.6) *specular;

    fragColor = vec4(intensity.x, intensity.y, intensity.z, 1.0);

}
