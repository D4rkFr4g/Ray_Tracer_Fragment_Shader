uniform samplerCube uTexUnit2;

varying vec3 vNormal;
varying vec3 vPosition;

vec3 reflect(vec3 w, vec3 n)
{
    return -w + n *(dot(w,n) *2.0);
}

void main(void)
{
    vec3 normal = normalize(vNormal);
    vec3 reflected = reflect(normalize(-vPosition), normal);
    gl_FragColor = textureCube(uTexUnit2, reflected);
}