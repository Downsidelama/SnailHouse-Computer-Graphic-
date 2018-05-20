#version 130

// pipeline-b�l bej�v� per-fragment attrib�tumok
in vec3 vs_out_pos;
in vec3 vs_out_normal;
in vec2 vs_out_tex0;

// kimen� �rt�k - a fragment sz�ne
out vec4 fs_out_col;

//
// uniform v�ltoz�k
//

uniform bool use_texture = false; // haszn�lunk-e text�r�t
uniform bool is_green = false; // konstans z�ld sz�nnel rajzoljunk-e

// sz�nt�r tulajdons�gok
uniform vec3 eye_pos = vec3(0, 15, 15);
uniform vec3 light_pos = vec3( 0, 5, 0 );
uniform vec4 La = vec4(0.2f, 0.2f, 0.2f, 1);

// f�nytulajdons�gok
uniform vec4 Ld = vec4(0.7f, 0.7f, 0.7f, 1);
uniform vec4 Ls = vec4(1, 1, 0.2, 1);

// anyagtulajdons�gok
uniform vec4 Ka = vec4(1, 1, 1, 1);
uniform vec4 Kd = vec4(1, 1, 1, 1);
uniform vec4 Ks = vec4(0.5, 0.5, 0.5, 1);
uniform float specular_power = 64;

uniform sampler2D texImage;

void main()
{
	
	//
	// ambiens sz�n sz�m�t�sa
	//
	vec4 ambient = La * Ka;

	//
	// diff�z sz�n sz�m�t�sa
	//

	/* seg�ts�g:
		- normaliz�l�s: http://www.opengl.org/sdk/docs/manglsl/xhtml/normalize.xml
	    - skal�ris szorzat: http://www.opengl.org/sdk/docs/manglsl/xhtml/dot.xml
	    - clamp: http://www.opengl.org/sdk/docs/manglsl/xhtml/clamp.xml
	*/

	vec3 normal = normalize(vs_out_normal);
	vec3 toLight = normalize(light_pos - vs_out_pos);
	float di = clamp( dot( toLight, normal), 0.0f, 1.0f );
	vec4 diffuse = Ld*Kd*di;

	//
	// f�nyfoltk�pz� sz�n
	//

	/* seg�ts�g:
		- reflect: http://www.opengl.org/sdk/docs/manglsl/xhtml/reflect.xml
		- power: http://www.opengl.org/sdk/docs/manglsl/xhtml/pow.xml
	*/
	vec4 specular = vec4(0);

	if ( di > 0 )
	{
		vec3 e = normalize( eye_pos - vs_out_pos );
		vec3 r = reflect( -toLight, normal );
		float si = pow( clamp( dot(e, r), 0.0f, 1.0f ), specular_power );
		specular = Ls*Ks*si;
	}
	
	// ha nem haszn�lunk text�r�t, akkor az alapsz�n a teljesen feh�r (1,1,1,1)
	if(!is_green){
		fs_out_col = (ambient + diffuse + specular) * ( use_texture ? texture(texImage, vs_out_tex0.st) : vec4(1.0f) );
	} else {
		// a norm�lvektor-p�lcik�kn�l se text�r�t, se f�nyeket nem sz�molunk: mindig ugyanolyan z�ldek
		fs_out_col = vec4(0, 1, 0, 1);
	}
}