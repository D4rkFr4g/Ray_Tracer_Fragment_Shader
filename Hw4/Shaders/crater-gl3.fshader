#version 130

uniform vec3 uLight;
uniform sampler2D uTexUnit2;
uniform int uSamples;
uniform float uSampledx, uSampledy;

in vec3 vNormal;
in vec3 vPosition;
in vec2 vTexCoord2;

out vec4 fragColor;

vec4 superSample(vec4 color)
{
	int count = 0;
	if (uSamples > 1)
	{
		int mid = int(ceil(uSamples / 2.0f));
		for (int y = -mid; y < mid; y++)
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
  

  float diffuse = max(0.0, dot(normal, tolight));
  fragColor = texColor0;// * diffuse;
}
