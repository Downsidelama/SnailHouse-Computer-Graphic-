#version 130

in vec3 vs_in_pos;
in vec3 vs_in_vel;

uniform mat4 mvp;
uniform vec3 eye_pos;
uniform vec3 at;

void main()
{
	if(dot(vs_in_pos - eye_pos, at - eye_pos) > (sqrt(2) * cos(90))) {
		gl_Position = mvp * vec4(vs_in_pos, 1);
	}
	

	// így írható felül a pontprimitívek mérete (feltéve ha a kliens oldalon van egy glEnable(GL_PROGRAM_POINT_SIZE);)
	//gl_PointSize = 5;
}