uniform vec3 uLight;
uniform sampler2D uTexUnit2;
uniform int uSamples;
uniform float uSampledx, uSampledy;

varying vec3 vNormal;
varying vec3 vPosition;
varying vec2 vTexCoord2;

vec4 superSample(vec4 color)
{
	int count = 0;
	if (uSamples > 1)
	{
		int half = int(ceil(uSamples / 2.0f));
		for (int y = -half; y < half; y++)
		{
			for (int x = -uSamples; x < uSamples; x++)
			{
				color += texture2D(uTexUnit2, vTexCoord2 + vec2(uSampledx * x, uSampledy * y));
				count++;
			}
		}

		color /= count;
	}
	else
	color = texture2D(uTexUnit2, vTexCoord2);

	return color;
}

vec4 overSample(vec4 color)
{
	int count = 0;
	if (uSamples > 1)
	{
		for (int i = 0; i < uSamples; i++)
		{
			color += texture2D(uTexUnit2, vTexCoord2 + vec2(uSampledx * i, uSampledy * i));
			color += texture2D(uTexUnit2, vTexCoord2 + vec2(-uSampledx * i, uSampledy * i));
			color += texture2D(uTexUnit2, vTexCoord2 + vec2(uSampledx * i, -uSampledy * i));
			color += texture2D(uTexUnit2, vTexCoord2 + vec2(-uSampledx * i, -uSampledy * i));
			
			count += 4;
		}

		color /= count;
	}
	else
		color = texture2D(uTexUnit2, vTexCoord2);

	return color;
}

void main() {
  vec3 tolight = normalize(uLight - vPosition);
  vec3 normal = normalize(vNormal);

  vec4 texColor0 = vec4(0);  
  
  //texColor0 = superSample(texColor0);
  texColor0 = overSample(texColor0);
  
  texColor0 = texture2D(uTexUnit2, vTexCoord2);

  float diffuse = max(0.0, dot(normal, tolight));
  gl_FragColor = texColor0;// * diffuse;
}
