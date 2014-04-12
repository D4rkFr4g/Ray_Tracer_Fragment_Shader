uniform vec3 uLight, uColor;

varying vec3 vNormal;
varying vec3 vPosition;
varying vec3 vTangent;

void main() 
{
    vec3 toLight = uLight - vec3(vPosition);
    vec3 toV = -normalize(vec3(vPosition));
    toLight = normalize(toLight);
    vec3 h = normalize(toV + toLight);
    vec3 normal = normalize(vNormal);

    vec3 vTangent3 = vec3(vTangent);
    vTangent3 = normalize(cross(vNormal, vTangent3));

    float nl = dot(normal, toLight);
    float dif = max(0., .75*nl + .25);
    float v = dot(vTangent3, h);
    v = pow(1.0 - v*v, 16.0);

    float r = uColor.r * dif + 0.3 * v;
    float g = uColor.g * dif + 0.3 * v;
    float b = uColor.b * dif + 0.3 * v;

    gl_FragColor = vec4(r, g, b, 1);
}
