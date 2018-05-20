#include "MyApp.h"
#include "GLUtils.hpp"

#include <GL/GLU.h>
#include <math.h>

#include <random>

#include "ObjParser_OGL3.h"

const float f_PI = 3.1415926f;

CMyApp::CMyApp(void)
{
	m_TextureID = 0;
}


CMyApp::~CMyApp(void)
{
}

// param�teres alak
glm::vec3 CMyApp::GetUV(float u, float v){

	glm::vec3 Fuv; // az F(u,v) param�teres f�ggv�ny�nk

	switch (shape_id){
		case 1:
			// T�RUSZ PARAM�TERES ALAK
			u *= -2 * f_PI;
			v *= 2 * f_PI;
			Fuv = glm::vec3( (2 + cosf(v))*cosf(u), sinf(v), (2 + cosf(v))*sinf(u) );
			break;
		case 2:
			// DOMBORZAT PARAM�TERES ALAK
			u = -4 * f_PI * u + 2 * f_PI;
			v =  4 * f_PI * v - 2 * f_PI;
			Fuv = glm::vec3(0.5f*u, sinf(u)*sinf(v), 0.5f*v);
			break;
		case 3:
			// V�ZA PARAM�TERES ALAK
			u = -8 * u + 4;
			v *= 2 * f_PI;
			Fuv = glm::vec3( u, cosf(u)*cosf(v), cosf(u)*sinf(v) );
			break;
		case 4:
			// TREFOIL KNOT PARAM�TERES ALAK
			u *= -2 * f_PI;
			v *=  2 * f_PI;
			Fuv = glm::vec3(
				cosf(2 * u) * cosf(v) + 3 * cosf(2 * u)*(1.5f + sinf(3 * u) / 2.0f),
				sinf(v) + cosf(3 * u),
				sinf(2 * u) * cosf(v) + 3 * sinf(2 * u)*(1.5f + sinf(3 * u) / 2.0f)
						   );
			break;
		case 5:
			// ZIZI-G�MB PARAM�TERES ALAK
			u *= -2 * f_PI;
			v *=  f_PI;
			Fuv = glm::vec3( sinf(u)*sinf(v) + 0.05f*cosf(20*v),
							 cosf(v),
							 cosf(u)*sinf(v) + 0.05f*cosf(20*u)
						   );
			break;
		case 6:
			// HULL�M PARAM�TERES ALAK
			u = -u * 50 + 25;
			v =  v * 50 - 25;
			Fuv = glm::vec3( 0.1f*u, 0.2f*cosf(glm::sqrt(u*u + v*v)), 0.1f*v );
			break;
		case 7:
			Fuv = glm::vec3((1 - expf(u))*cosf(6 * M_PI*u)*cosf(M_PI*v)*cosf(M_PI*v), (-1 + expf(u))*sinf(6 * M_PI*u)*cos(M_PI*v)*cosf(M_PI*v),1 - expf(2 * u) - sinf(2 * M_PI*v) + expf(u)*sin(2 * M_PI*v));
			Fuv = glm::vec3((1 - expf(u))*cosf(6 * M_PI*u)*cosf(M_PI*v)*cosf(M_PI*v), 1 - expf(2 * u) - sinf(2 * M_PI*v) + expf(u)*sin(2 * M_PI*v), (-1 + expf(u))*sinf(6 * M_PI*u)*cos(M_PI*v)*cosf(M_PI*v));
			break;
	}

	return Fuv;
}

// norm�lvektorok sz�mol�sa
glm::vec3 CMyApp::GetUV_normals(float u, float v){

	glm::vec3 du, dv; // u �s v szerinti parci�lis deriv�ltak
	glm::vec3 n; // a v�gs� norm�lvektor

	// az alakzat param�teres alakja az F(u,v) f�ggv�ny (l�sd: GetUV() elj�r�s)
	// Majdnem b�rmilyen (kell�en sima �s folytonos) param�teres fel�let norm�lvektor�t
	//   megkaphatjuk egy adott (u,v) pontban az al�bbi m�don:
	// n(u,v) = ( dF(u,v)/du x dF(u,v)/dv ) / | dF(u,v)/du x dF(u,v)/dv |

	// Azaz venni kell az F(u,v) vektor u �s v szerinti parci�lis deriv�ltj�t,
	// �s a kapott vektorok kereszt-szorzata (egy normaliz�l�s ut�n) lesz a norm�lvektor.


	// El�sz�r kisz�moljuk du-t �s dv-t
	switch (shape_id){
		case 1:
			// T�RUSZ PARCI�LIS DERIV�LTAK
			u *= -2 * f_PI;
			v *= 2 * f_PI;

			// A t�rusz param�teres alakja:
			// Fuv = glm::vec3((2 + cosf(v))*cosf(u), sinf(v), (2 + cosf(v))*sinf(u));
			// Sz�moljuk ki Fuv parci�lis deriv�ltjait!

			// (2 + cosf(v))*cosf(u)	deriv�ltja u szerint: (2 + cosf(v))*-sinf(u)
			// sinf(v)					deriv�ltja u szerint: 0
			// (2 + cosf(v))*sinf(u)	deriv�ltja u szerint: (2 + cosf(v))*cosf(u)
			du = glm::vec3( (2 + cosf(v))*-sinf(u), 0, (2 + cosf(v))*cosf(u) );

			// (2 + cosf(v))*cosf(u)	deriv�ltja v szerint: -sinf(v)*cosf(u)
			// sinf(v)					deriv�ltja v szerint: cosf(v)
			// (2 + cosf(v))*sinf(u)	deriv�ltja v szerint: -sinf(v)*sinf(u)
			dv = glm::vec3( -sinf(v)*cosf(u), cosf(v), -sinf(v)*sinf(u) );
			break;
		case 2:
			// DOMBORZAT PARCI�LIS DERIV�LTAK
			u = -4 * f_PI * u + 2 * f_PI;
			v = 4 * f_PI * v - 2 * f_PI;
			// Fuv = glm::vec3(0.5f*u, sinf(u)*sinf(v), 0.5f*v);
			du = glm::vec3(0.5f, cosf(u)*sinf(v), 0);
			dv = glm::vec3(0,	 sinf(u)*cosf(v), 0.5f);
			break;
		case 3:
			// V�ZA PARCI�LIS DERIV�LTAK
			u = -8 * u + 4;
			v *= 2 * f_PI;
			// Fuv = glm::vec3( u, cosf(u)*cosf(v), cosf(u)*sinf(v) );
			du = glm::vec3( 1, -sinf(u)*cosf(v), -sinf(u)*sinf(v) );
			dv = glm::vec3( 0, cosf(u)*-sinf(v),  cosf(u)*cosf(v) );
			break;
		case 4:
			// TREFOIL KNOT PARCI�LIS DERIV�LTAK
			u *= -2 * f_PI;
			v *= 2 * f_PI;
			/* Fuv = glm::vec3(
				cosf(2 * u) * cosf(v) + 3 * cosf(2 * u)*(1.5f + sinf(3 * u) / 2.0f),
				sinf(v) + cosf(3 * u),
				sinf(2 * u) * cosf(v) + 3 * sinf(2 * u)*(1.5f + sinf(3 * u) / 2.0f)
				);
			*/
			du = glm::vec3(
				-2*sinf(2*u)*cosf(v) - 6*sinf(2*u)*(sinf(3*u) / 2.0f + 1.5f) + 4.5f*cosf(2*u)*cosf(3*u),
				-3*sinf(3*u),
				 2*cosf(2*u)*cosf(v) + 6*cosf(2*u)*(sinf(3*u) / 2.0f + 1.5f) + 4.5f*sinf(2*u)*cosf(3*u)
			);
			dv = glm::vec3(
				-cosf(2*u)*sinf(v),
				cosf(v),
				-sinf(2*u)*sinf(v)
			);
			break;
		case 5:
			// ZIZI-G�MB PARCI�LIS DERIV�LTAK
			u *= -2 * f_PI;
			v *= f_PI;
			// Fuv = glm::vec3(sinf(u)*sinf(v) + 0.05f*cosf(20 * v),
			// 	cosf(v),
			// 	cosf(u)*sinf(v) + 0.05f*cosf(20 * u)
			// 	);
			du = glm::vec3( cosf(u)*sinf(v),			  0,		-sinf(u)*sinf(v) - sinf(20*u) );
			dv = glm::vec3( sinf(u)*cosf(v) - sinf(20*v), -sinf(v), cosf(u)*cosf(v) );
			break;
		case 6:
			// HULL�M PARCI�LIS DERIV�LTAK
			u = -u * 50 + 25;
			v = v * 50 - 25;
			// Fuv = glm::vec3( 0.1f*u, 0.2f*cosf(glm::sqrt(u*u + v*v)), 0.1f*v );
			// apr� "csal�s": az y �rt�k eset�n a nevez� lehet 0 (ha u = v = 0),
			// ez�rt egy 0.001f konstanst hozz�adtunk, hogy elker�lj�k a 0-val val� oszt�st
			du = glm::vec3(0.1f, -0.2f*u*sinf(glm::sqrt(u*u + v*v)) / ((glm::sqrt(u*u + v*v) + 0.001f)), 0);
			dv = glm::vec3(0,    -0.2f*v*sinf(glm::sqrt(u*u + v*v)) / ((glm::sqrt(u*u + v*v) + 0.001f)), 0.1f);
			break;
		case 7:
			du = glm::vec3(powf(-cosf(M_PI*v),2)*(expf(u)*cos(6*M_PI*u)-6*M_PI*(expf(u)-1)*sin(6*M_PI*u)), powf(cos(M_PI*v), 2)*(expf(u)*sinf(6*M_PI*u)+6*M_PI*(expf(u)-1)*cosf(6*M_PI*u)), expf(u)*(sin(2*M_PI*v)-2*expf(u)));
			dv = glm::vec3(2*M_PI*(expf(u)-1)*cosf(6*M_PI*u)*sin(M_PI*v)*cos(M_PI*v), -2*M_PI*(expf(u)-1)*sin(6*M_PI*u)*sin(M_PI*v)*cos(M_PI*v), 2*M_PI*(expf(u)-1)*cos(2*M_PI*v));


			// Helyes
			du = glm::vec3(powf(-cosf(M_PI*v), 2)*(expf(u)*cos(6 * M_PI*u) - 6 * M_PI*(expf(u) - 1)*sin(6 * M_PI*u)), expf(u)*(sin(2 * M_PI*v) - 2 * expf(u)), powf(cos(M_PI*v), 2)*(expf(u)*sinf(6 * M_PI*u) + 6 * M_PI*(expf(u) - 1)*cosf(6 * M_PI*u)));
			dv = glm::vec3(2 * M_PI*(expf(u) - 1)*cosf(6 * M_PI*u)*sin(M_PI*v)*cos(M_PI*v), 2 * M_PI*(expf(u) - 1)*cos(2 * M_PI*v), -2 * M_PI*(expf(u) - 1)*sin(6 * M_PI*u)*sin(M_PI*v)*cos(M_PI*v));

			break;
	}

	// most m�r adottak a parci�lis deriv�ltak (du �s dv)
	// a kett� kereszt-szorzata adja meg a fel�leti norm�list
	// (�gyelj�nk az ir�ny�totts�gra!)
	n = glm::cross(dv, du);
	if (shape_id == 7) {
		n = -n;
	}
	// normaliz�ljuk egys�ghossz�ra
	if (glm::length(n) > 0.0f){
		n = glm::normalize(n);
	}
	// pr�b�ljuk ki, mit l�tunk, ha mondjuk minden norm�lvektort (0, 1, 0)-ra �ll�tunk!
	//  return glm::vec3(0, 1, 0);
	// ugye milyen fura?

	return n;


}

bool CMyApp::Init()
{
	lastTickTime = SDL_GetTicks();
	// t�rl�si sz�n legyen k�kes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE);		// kapcsoljuk be a hatrafele nezo lapok eldobasat
	glEnable(GL_DEPTH_TEST);	// m�lys�gi teszt bekapcsol�sa (takar�s)
	glPolygonMode(GL_BACK, GL_LINE);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(3.0f);

	//
	// geometria letrehozasa
	//

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// TALAJ
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// attrib�tumok
	m_vb_ground.AddAttribute(0, 3);
	m_vb_ground.AddAttribute(1, 3);
	m_vb_ground.AddAttribute(2, 2);

	// 0. attrib�tum lista: poz�ci�k
	m_vb_ground.AddData(0, -10, 0, -10);
	m_vb_ground.AddData(0,  10, 0, -10);
	m_vb_ground.AddData(0, -10, 0,  10);
	m_vb_ground.AddData(0,  10, 0,  10);

	// 1. attrib�tum lista: norm�lvektorok (milyen j�, hogy az XZ s�kon fekv� talajnak trivi�lisak...)
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);

	// 2. attrib�tum: text�ra koordin�t�k
	m_vb_ground.AddData(2, 0, 0);
	m_vb_ground.AddData(2, 1, 0);
	m_vb_ground.AddData(2, 0, 1);
	m_vb_ground.AddData(2, 1, 1);

	// vertexek indexel�se
	m_vb_ground.AddIndex(1, 0, 2);
	m_vb_ground.AddIndex(1, 2, 3);

	// VAO, VBO �s IB inicializ�l�sa
	m_vb_ground.InitBuffers();

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ALAKZAT
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	m_vb_shape.AddAttribute(0, 3);
	m_vb_shape.AddAttribute(1, 3);
	m_vb_shape.AddAttribute(2, 2);

	// N*M-es r�cs bej�r�sa, a vertexek (x,y,z) koordin�t�j�t a GetUV f�ggv�nnyel kapjuk
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			float u = i / (float)N;
			float v = j / (float)M;
			// poz�ci�k: a param�teres f�ggv�ny szerint
			m_vb_shape.AddData(0, GetUV(u, v));
			// norm�lvektorok: a param�teres f�ggv�nyb�l sz�molva
			m_vb_shape.AddData(1, GetUV_normals(u, v));
			// text�ra koordin�t�k
			// itt most 1-n�l nagyobb sz�mokat is megenged�nk: ekkor a text�ra minta ism�tl�dni fog
			m_vb_shape.AddData(2, 4*(1 - u), 4*v); // a text�ra-koordin�t�kat m�s sorrendben kell bej�rni!
		}
	}
	// A norm�lvektorokat form�l� z�ld vonalak "m�sik v�ge" (�j vertexek)
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			float u = i / (float)N;
			float v = j / (float)M;
			// a poz�ci�k: az (u,v) ponthoz tartoz� vertex poz�ci�j�t eltoljuk a norm�lvektor ir�ny�ba
			// a 0.2-es szorz� csak az�rt kell, mert az egys�ghossz� p�lcik�k t�l hossz�ak lenn�nek (helyette legyenek mondjuk 0.2 hossz�ak)
			m_vb_shape.AddData(0, GetUV(u, v) + 0.2f*GetUV_normals(u, v));
			m_vb_shape.AddData(1, glm::vec3(0, 0, 0)); // ez itt mindegy - a norm�lvektorokra nem fogunk f�nyeket sz�molni
			m_vb_shape.AddData(2, glm::vec2(0, 0)); // ez is - text�r�zni sem fogjuk �ket
		}
	}


	// az indexek meghat�roz�sa minden n�gyzethez
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			// az egyik h�romsz�g a n�gyzet egyik fele
			m_vb_shape.AddIndex( (i)		+ (j)*		(N + 1),
								 (i)		+ (j + 1) *	(N + 1),
								 (i + 1)	+ (j) *		(N + 1)
							 );
			// a m�sik h�romsz�g a n�gyzet m�sik fele
			m_vb_shape.AddIndex( (i + 1)	+ (j)*		(N + 1),
								 (i)		+ (j + 1) *	(N + 1),
								 (i + 1)	+ (j + 1) *	(N + 1)
							 );
		}
	}
	// a norm�lvektorokat form�l� z�ld vonalak indexei:
	// a parametrikus fel�let�nk�n az adott (i,j) vertex indexe:
	//		(ez lesz a z�ld vonal egyik v�gpontja)
	//		P = i + j*(N + 1)
	//		
	// a hozz� tartoz�, norm�lvektor szerint eltolt vertex poz�ci� indexe:
	//		fel�letet_alkot�_vertexsz�m + P = (N + 1)*(M + 1) + i + j*(N + 1)
	//		(ez lesz a z�ld vonal m�sik v�gpontja)
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			m_vb_shape.AddIndex(i + j*(N + 1), (N + 1)*(M + 1) + i + j*(N + 1));
		}
	}

	m_vb_shape.InitBuffers();

	//
	// shaderek bet�lt�se
	//
	m_program.AttachShader(GL_VERTEX_SHADER, "dirLight.vert");
	m_program.AttachShader(GL_FRAGMENT_SHADER, "dirLight.frag");

	m_program.BindAttribLoc(0, "vs_in_pos");
	m_program.BindAttribLoc(1, "vs_in_normal");
	m_program.BindAttribLoc(2, "vs_in_tex0");

	if ( !m_program.LinkProgram() )
	{
		return false;
	}


	// r�szecsk�ket kirajzol� program
	m_particleProgram.Init({	// shaderek felsorol�sa
		{ GL_VERTEX_SHADER,		"particle.vert" },
		{ GL_FRAGMENT_SHADER,	"particle.frag" }
		},
	{	// binding-ok felsorol�sa
		{ 0, "vs_in_pos" },
		{ 1, "vs_in_vel" },
	});

	// r�szecsk�k inicializ�l�sa
	m_particlePos.reserve(m_particleCount);
	m_particleVel.reserve(m_particleCount);

	// v�letlensz�m gener�tor inicializ�l�sa
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> rnd(-20, 20);
	std::uniform_real_distribution<> rndVelocity(5, 10);

	// CPU oldali r�szecsket�mb�k felt�lt�se
	for (int i = 0; i < m_particleCount; ++i)
	{
		m_particlePos.push_back(glm::vec3(rnd(gen), 10, rnd(gen)));
		m_particleVel.push_back(glm::vec3(0, 2 * -abs(rndVelocity(gen)), 0));
	}

	// GPU-ra �tt�lteni a r�szecsk�k poz�ci�it
	m_gpuParticleBuffer.BufferData(m_particlePos);	// <=>	m_gpuParticleBuffer = m_particlePos;

													// �s v�g�l a VAO-t inicializ�lni
	m_gpuParticleVAO.Init({
		{ CreateAttribute<0, glm::vec3, 0, sizeof(glm::vec3)>, m_gpuParticleBuffer }
		});

	//
	// egy�b inicializ�l�s
	//

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);

	// text�r�k bet�lt�se
	m_TextureID = TextureFromFile("ground_texture.bmp");

	return true;
}

void CMyApp::Clean()
{
	glDeleteTextures(1, &m_TextureID);

	m_program.Clean();
}

void CMyApp::Update()
{
	// kameramozgat�shoz
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time)/1000.0f;
	m_camera.Update(delta_time);
	glm::vec3 p1 = glm::vec3(-12.0f, 20.0f, 0.0f);
	glm::vec3 p2 = glm::vec3(0.0f, 20.0f, -15.0f);
	glm::vec3 p3 = glm::vec3(12.0f, 20.0f, 0.0f);

	// Id� � [0,1]
	if (count) {
		currentTime += SDL_GetTicks() - lastTickTime;
	}
	float t = abs(sinf(currentTime / 1000.0f / 6 * M_PI));

	// Bezier g�rbe
	float x = ((1.0f - t)*((1.0f - t)*p1.x + t * p2.x) + t * ((1.0f - t) * p2.x + t * p3.x));
	float y = ((1.0f - t)*((1.0f - t)*p1.y + t * p2.y) + t * ((1.0f - t) * p2.y + t * p3.y));
	float z = ((1.0f - t)*((1.0f - t)*p1.z + t * p2.z) + t * ((1.0f - t) * p2.z + t * p3.z));

	m_camera.SetView(glm::vec3(x, y, z),		// honnan n�zz�k a sz�nteret
		glm::vec3(0, 0, 0),		// a sz�nt�r melyik pontj�t n�zz�k
		glm::vec3(0, 1, 0));

	// friss�ts�k a poz�ci�kat
	static const float energyRemaining = 1;	// t�k�letesen rugalmas �tk�z�s
	for (int i = 0; i < m_particleCount; ++i)
	{
		m_particlePos[i] += m_particleVel[i] * delta_time;

		if (m_particlePos[i].y <= 0) {
			m_particlePos[i].y = 20;
		}

	}

	// friss�ts�k a puffert
	glBindBuffer(GL_ARRAY_BUFFER, m_gpuParticleBuffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3)*m_particlePos.size(), &(m_particlePos[0][0]));
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	last_time = SDL_GetTicks();
	lastTickTime = SDL_GetTicks();
}

// talaj kirajzol�sa
void CMyApp::DrawGround(){

	glm::mat4 matWorld = glm::mat4(1.0f); // nem v�ltoztatunk rajta
	glm::mat4 matWorldIT = glm::transpose(glm::inverse(matWorld));
	glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

	m_program.SetUniform("world", matWorld);
	m_program.SetUniform("worldIT", matWorldIT);
	m_program.SetUniform("MVP", mvp);

	m_program.SetTexture("texImage", 0, m_TextureID);

	m_vb_ground.On();
	m_vb_ground.DrawIndexed(GL_TRIANGLES, 0, 6, 0);
	m_vb_ground.Off();

}

void CMyApp::DrawShape(){

	// k�tszeres�re nagy�tjuk �s feljebb toljuk, ne l�gjon a talajba
	glm::mat4 matWorld = glm::translate(glm::vec3(0.0f, 13.0f, 0.0f)) * glm::scale<float>(glm::vec3(2.0f, 2.0f, 2.0f));
	glm::mat4 matWorldIT = glm::transpose(glm::inverse(matWorld));
	//glm::mat4 mvp = m_camera.GetViewProj() *matWorld;
	glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

	m_program.SetUniform("world", matWorld);
	m_program.SetUniform("worldIT", matWorldIT);
	m_program.SetUniform("MVP", mvp);

	m_program.SetTexture("texImage", 0, m_TextureID);

	m_vb_shape.On();
	// kirajzoljuk a fel�letet: 0-t�l a k�vetkez� K = 6*N*M indexet figyelembe v�ve...
	m_vb_shape.DrawIndexed(GL_TRIANGLES, 0, 2 * 3 * N * M, 0);
	// ...�s a z�ld vonalakat is: GL_LINES a K-adik indext�l, a k�vetkez� 2*P index figyelembe v�tel�vel,
	//		ahol P = (N+1)*(M+1) a fel�letet alkot� vertexek sz�ma
	//		GL_LINES az indexek �ltal jel�lt vertexeket p�ronk�nt vonalakk� k�ti �ssze
	m_program.SetUniform("is_green", true); // sz�lunk a shadernek, hogy ne sz�moljon f�nyeket, hanem konstans z�ld sz�nnel rajzoljon
	//m_vb_shape.DrawIndexed(GL_LINES, 0, 2 * (N + 1)*(M + 1), (void*)(2 * 3 * N*M*sizeof(unsigned int)));
	m_program.SetUniform("is_green", false); // ne felejts�k el vissza�ll�tani a megjelen�t�st shaderben
	m_vb_shape.Off();

}

void CMyApp::Render()
{
	// t�r�lj�k a frampuffert (GL_COLOR_BUFFER_BIT) �s a m�lys�gi Z puffert (GL_DEPTH_BUFFER_BIT)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// megjelen�t�s m�dja
	if (is_filled){
		glPolygonMode(GL_FRONT, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT, GL_LINE);
	}

	m_program.On(); // shader program bekapcs

	// light_pos a pontszer� f�nyforr�sunk poz�ci�ja;
	// y = 10 magass�gban kering az XZ s�kon, 5 mp alatt megt�ve egy 10 sugar� k�rp�ly�t
	float r = 2 * f_PI * SDL_GetTicks() / 1000.0f / 5.0f;
	glm::vec3 light_pos = glm::vec3(10 * cosf(r), 10, 10 * sinf(r));
	
	// a shadernek �tadjuk a rajzol�shoz sz�ks�ges param�tereket
	// hol a f�ny �s a kamera
	//m_program.SetUniform("light_pos", light_pos);
	m_program.SetUniform("eye_pos", m_camera.GetEye());

	// haszn�lunk-e text�r�t
	m_program.SetUniform("use_texture", use_texture);
	
	// kirajzoljuk az alakzatokat
	//DrawGround();
	DrawShape();

	m_program.Off(); // shader program kikapcs

	// r�szecsk�k
	glEnable(GL_PROGRAM_POINT_SIZE);
	m_gpuParticleVAO.Bind();
	m_particleProgram.Use();
	m_particleProgram.SetUniform("mvp", m_camera.GetViewProj());
	m_particleProgram.SetUniform("eye_pos", m_camera.GetEye());
	m_particleProgram.SetUniform("at", m_camera.GetAt());

	glDrawArrays(GL_POINTS, 0, m_particleCount);

	glDisable(GL_PROGRAM_POINT_SIZE);

}

void CMyApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	//m_camera.KeyboardDown(key);
	switch (key.keysym.sym){
	case SDLK_q:
		//is_filled = !is_filled;
		break;
	case SDLK_t:
		//use_texture = !use_texture;
		break;
	case SDLK_SPACE:
		count = !count;
		break;
	}
}

void CMyApp::KeyboardUp(SDL_KeyboardEvent& key)
{
	//m_camera.KeyboardUp(key);
}

void CMyApp::MouseMove(SDL_MouseMotionEvent& mouse)
{
	//m_camera.MouseMove(mouse);
}

void CMyApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

// a k�t param�terbe az �j ablakm�ret sz�less�ge (_w) �s magass�ga (_h) tal�lhat�
void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_camera.Resize(_w, _h);
}