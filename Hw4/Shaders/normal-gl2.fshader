uniform vec3 uLight, uColor;
uniform sampler2D uTexUnit1;

varying vec3 vPosition;
varying vec2 vTexCoord1;
varying vec3 vNormal;

void main() 
{
    vec3 toLight = normalize(uLight - vPosition);
    vec3 toV = -normalize(vec3(vPosition));
    vec3 h = normalize(toV + toLight);

    vec4 texNormal = texture2D(uTexUnit1, vTexCoord1);
    vec3 scaledTexNormal = normalize(2.0 * vec3(texNormal) + vNormal);

    float specular = pow(max(0.0, dot(h, scaledTexNormal)), 64.0);
    float diffuse = max(0.0, dot(vec3(scaledTexNormal), toLight));
    vec3 intensity = vec3(0.001, 0.001, 0.001) + uColor * diffuse
        + vec3(0.6, 0.6, 0.6) *specular;
    
    gl_FragColor = vec4(intensity.x, intensity.y, intensity.z, 1.0);
    //gl_FragColor = vec4(0.7, 0.7, 0.7, 1.0) * diffuse;
}