#version 130

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_normal;
in vec2 vs_out_tex0;

// kimenõ érték - a fragment színe
out vec4 fs_out_col;

//
// uniform változók
//

uniform bool use_texture = false; // használunk-e textúrát
uniform bool is_green = false; // konstans zöld színnel rajzoljunk-e

// színtér tulajdonságok
uniform vec3 eye_pos = vec3(0, 15, 15);
uniform vec3 light_pos = vec3( 0, 5, 0 );
uniform vec4 La = vec4(0.2f, 0.2f, 0.2f, 1);

// fénytulajdonságok
uniform vec4 Ld = vec4(0.7f, 0.7f, 0.7f, 1);
uniform vec4 Ls = vec4(1, 1, 0.2, 1);

// anyagtulajdonságok
uniform vec4 Ka = vec4(1, 1, 1, 1);
uniform vec4 Kd = vec4(1, 1, 1, 1);
uniform vec4 Ks = vec4(0.5, 0.5, 0.5, 1);
uniform float specular_power = 64;

uniform sampler2D texImage;

void main()
{
	
	//
	// ambiens szín számítása
	//
	vec4 ambient = La * Ka;

	//
	// diffúz szín számítása
	//

	/* segítség:
		- normalizálás: http://www.opengl.org/sdk/docs/manglsl/xhtml/normalize.xml
	    - skaláris szorzat: http://www.opengl.org/sdk/docs/manglsl/xhtml/dot.xml
	    - clamp: http://www.opengl.org/sdk/docs/manglsl/xhtml/clamp.xml
	*/

	vec3 normal = normalize(vs_out_normal);
	vec3 toLight = normalize(light_pos - vs_out_pos);
	float di = clamp( dot( toLight, normal), 0.0f, 1.0f );
	vec4 diffuse = Ld*Kd*di;

	//
	// fényfoltképzõ szín
	//

	/* segítség:
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
	
	// ha nem használunk textúrát, akkor az alapszín a teljesen fehér (1,1,1,1)
	if(!is_green){
		fs_out_col = (ambient + diffuse + specular) * ( use_texture ? texture(texImage, vs_out_tex0.st) : vec4(1.0f) );
	} else {
		// a normálvektor-pálcikáknál se textúrát, se fényeket nem számolunk: mindig ugyanolyan zöldek
		fs_out_col = vec4(0, 1, 0, 1);
	}
}