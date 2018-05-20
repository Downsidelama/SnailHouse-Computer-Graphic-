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

// paraméteres alak
glm::vec3 CMyApp::GetUV(float u, float v){

	glm::vec3 Fuv; // az F(u,v) paraméteres függvényünk

	switch (shape_id){
		case 1:
			// TÓRUSZ PARAMÉTERES ALAK
			u *= -2 * f_PI;
			v *= 2 * f_PI;
			Fuv = glm::vec3( (2 + cosf(v))*cosf(u), sinf(v), (2 + cosf(v))*sinf(u) );
			break;
		case 2:
			// DOMBORZAT PARAMÉTERES ALAK
			u = -4 * f_PI * u + 2 * f_PI;
			v =  4 * f_PI * v - 2 * f_PI;
			Fuv = glm::vec3(0.5f*u, sinf(u)*sinf(v), 0.5f*v);
			break;
		case 3:
			// VÁZA PARAMÉTERES ALAK
			u = -8 * u + 4;
			v *= 2 * f_PI;
			Fuv = glm::vec3( u, cosf(u)*cosf(v), cosf(u)*sinf(v) );
			break;
		case 4:
			// TREFOIL KNOT PARAMÉTERES ALAK
			u *= -2 * f_PI;
			v *=  2 * f_PI;
			Fuv = glm::vec3(
				cosf(2 * u) * cosf(v) + 3 * cosf(2 * u)*(1.5f + sinf(3 * u) / 2.0f),
				sinf(v) + cosf(3 * u),
				sinf(2 * u) * cosf(v) + 3 * sinf(2 * u)*(1.5f + sinf(3 * u) / 2.0f)
						   );
			break;
		case 5:
			// ZIZI-GÖMB PARAMÉTERES ALAK
			u *= -2 * f_PI;
			v *=  f_PI;
			Fuv = glm::vec3( sinf(u)*sinf(v) + 0.05f*cosf(20*v),
							 cosf(v),
							 cosf(u)*sinf(v) + 0.05f*cosf(20*u)
						   );
			break;
		case 6:
			// HULLÁM PARAMÉTERES ALAK
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

// normálvektorok számolása
glm::vec3 CMyApp::GetUV_normals(float u, float v){

	glm::vec3 du, dv; // u és v szerinti parciális deriváltak
	glm::vec3 n; // a végsõ normálvektor

	// az alakzat paraméteres alakja az F(u,v) függvény (lásd: GetUV() eljárás)
	// Majdnem bármilyen (kellõen sima és folytonos) paraméteres felület normálvektorát
	//   megkaphatjuk egy adott (u,v) pontban az alábbi módon:
	// n(u,v) = ( dF(u,v)/du x dF(u,v)/dv ) / | dF(u,v)/du x dF(u,v)/dv |

	// Azaz venni kell az F(u,v) vektor u és v szerinti parciális deriváltját,
	// és a kapott vektorok kereszt-szorzata (egy normalizálás után) lesz a normálvektor.


	// Elõször kiszámoljuk du-t és dv-t
	switch (shape_id){
		case 1:
			// TÓRUSZ PARCIÁLIS DERIVÁLTAK
			u *= -2 * f_PI;
			v *= 2 * f_PI;

			// A tórusz paraméteres alakja:
			// Fuv = glm::vec3((2 + cosf(v))*cosf(u), sinf(v), (2 + cosf(v))*sinf(u));
			// Számoljuk ki Fuv parciális deriváltjait!

			// (2 + cosf(v))*cosf(u)	deriváltja u szerint: (2 + cosf(v))*-sinf(u)
			// sinf(v)					deriváltja u szerint: 0
			// (2 + cosf(v))*sinf(u)	deriváltja u szerint: (2 + cosf(v))*cosf(u)
			du = glm::vec3( (2 + cosf(v))*-sinf(u), 0, (2 + cosf(v))*cosf(u) );

			// (2 + cosf(v))*cosf(u)	deriváltja v szerint: -sinf(v)*cosf(u)
			// sinf(v)					deriváltja v szerint: cosf(v)
			// (2 + cosf(v))*sinf(u)	deriváltja v szerint: -sinf(v)*sinf(u)
			dv = glm::vec3( -sinf(v)*cosf(u), cosf(v), -sinf(v)*sinf(u) );
			break;
		case 2:
			// DOMBORZAT PARCIÁLIS DERIVÁLTAK
			u = -4 * f_PI * u + 2 * f_PI;
			v = 4 * f_PI * v - 2 * f_PI;
			// Fuv = glm::vec3(0.5f*u, sinf(u)*sinf(v), 0.5f*v);
			du = glm::vec3(0.5f, cosf(u)*sinf(v), 0);
			dv = glm::vec3(0,	 sinf(u)*cosf(v), 0.5f);
			break;
		case 3:
			// VÁZA PARCIÁLIS DERIVÁLTAK
			u = -8 * u + 4;
			v *= 2 * f_PI;
			// Fuv = glm::vec3( u, cosf(u)*cosf(v), cosf(u)*sinf(v) );
			du = glm::vec3( 1, -sinf(u)*cosf(v), -sinf(u)*sinf(v) );
			dv = glm::vec3( 0, cosf(u)*-sinf(v),  cosf(u)*cosf(v) );
			break;
		case 4:
			// TREFOIL KNOT PARCIÁLIS DERIVÁLTAK
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
			// ZIZI-GÖMB PARCIÁLIS DERIVÁLTAK
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
			// HULLÁM PARCIÁLIS DERIVÁLTAK
			u = -u * 50 + 25;
			v = v * 50 - 25;
			// Fuv = glm::vec3( 0.1f*u, 0.2f*cosf(glm::sqrt(u*u + v*v)), 0.1f*v );
			// apró "csalás": az y érték esetén a nevezõ lehet 0 (ha u = v = 0),
			// ezért egy 0.001f konstanst hozzáadtunk, hogy elkerüljük a 0-val való osztást
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

	// most már adottak a parciális deriváltak (du és dv)
	// a kettõ kereszt-szorzata adja meg a felületi normálist
	// (Ügyeljünk az irányítottságra!)
	n = glm::cross(dv, du);
	if (shape_id == 7) {
		n = -n;
	}
	// normalizáljuk egységhosszúra
	if (glm::length(n) > 0.0f){
		n = glm::normalize(n);
	}
	// próbáljuk ki, mit látunk, ha mondjuk minden normálvektort (0, 1, 0)-ra állítunk!
	//  return glm::vec3(0, 1, 0);
	// ugye milyen fura?

	return n;


}

bool CMyApp::Init()
{
	lastTickTime = SDL_GetTicks();
	// törlési szín legyen kékes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE);		// kapcsoljuk be a hatrafele nezo lapok eldobasat
	glEnable(GL_DEPTH_TEST);	// mélységi teszt bekapcsolása (takarás)
	glPolygonMode(GL_BACK, GL_LINE);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(3.0f);

	//
	// geometria letrehozasa
	//

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// TALAJ
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// attribútumok
	m_vb_ground.AddAttribute(0, 3);
	m_vb_ground.AddAttribute(1, 3);
	m_vb_ground.AddAttribute(2, 2);

	// 0. attribútum lista: pozíciók
	m_vb_ground.AddData(0, -10, 0, -10);
	m_vb_ground.AddData(0,  10, 0, -10);
	m_vb_ground.AddData(0, -10, 0,  10);
	m_vb_ground.AddData(0,  10, 0,  10);

	// 1. attribútum lista: normálvektorok (milyen jó, hogy az XZ síkon fekvõ talajnak triviálisak...)
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);
	m_vb_ground.AddData(1, 0, 1, 0);

	// 2. attribútum: textúra koordináták
	m_vb_ground.AddData(2, 0, 0);
	m_vb_ground.AddData(2, 1, 0);
	m_vb_ground.AddData(2, 0, 1);
	m_vb_ground.AddData(2, 1, 1);

	// vertexek indexelése
	m_vb_ground.AddIndex(1, 0, 2);
	m_vb_ground.AddIndex(1, 2, 3);

	// VAO, VBO és IB inicializálása
	m_vb_ground.InitBuffers();

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ALAKZAT
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	m_vb_shape.AddAttribute(0, 3);
	m_vb_shape.AddAttribute(1, 3);
	m_vb_shape.AddAttribute(2, 2);

	// N*M-es rács bejárása, a vertexek (x,y,z) koordinátáját a GetUV függvénnyel kapjuk
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			float u = i / (float)N;
			float v = j / (float)M;
			// pozíciók: a paraméteres függvény szerint
			m_vb_shape.AddData(0, GetUV(u, v));
			// normálvektorok: a paraméteres függvénybõl számolva
			m_vb_shape.AddData(1, GetUV_normals(u, v));
			// textúra koordináták
			// itt most 1-nél nagyobb számokat is megengedünk: ekkor a textúra minta ismétlõdni fog
			m_vb_shape.AddData(2, 4*(1 - u), 4*v); // a textúra-koordinátákat más sorrendben kell bejárni!
		}
	}
	// A normálvektorokat formáló zöld vonalak "másik vége" (új vertexek)
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			float u = i / (float)N;
			float v = j / (float)M;
			// a pozíciók: az (u,v) ponthoz tartozó vertex pozícióját eltoljuk a normálvektor irányába
			// a 0.2-es szorzó csak azért kell, mert az egységhosszú pálcikák túl hosszúak lennének (helyette legyenek mondjuk 0.2 hosszúak)
			m_vb_shape.AddData(0, GetUV(u, v) + 0.2f*GetUV_normals(u, v));
			m_vb_shape.AddData(1, glm::vec3(0, 0, 0)); // ez itt mindegy - a normálvektorokra nem fogunk fényeket számolni
			m_vb_shape.AddData(2, glm::vec2(0, 0)); // ez is - textúrázni sem fogjuk õket
		}
	}


	// az indexek meghatározása minden négyzethez
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			// az egyik háromszög a négyzet egyik fele
			m_vb_shape.AddIndex( (i)		+ (j)*		(N + 1),
								 (i)		+ (j + 1) *	(N + 1),
								 (i + 1)	+ (j) *		(N + 1)
							 );
			// a másik háromszög a négyzet másik fele
			m_vb_shape.AddIndex( (i + 1)	+ (j)*		(N + 1),
								 (i)		+ (j + 1) *	(N + 1),
								 (i + 1)	+ (j + 1) *	(N + 1)
							 );
		}
	}
	// a normálvektorokat formáló zöld vonalak indexei:
	// a parametrikus felületünkön az adott (i,j) vertex indexe:
	//		(ez lesz a zöld vonal egyik végpontja)
	//		P = i + j*(N + 1)
	//		
	// a hozzá tartozó, normálvektor szerint eltolt vertex pozíció indexe:
	//		felületet_alkotó_vertexszám + P = (N + 1)*(M + 1) + i + j*(N + 1)
	//		(ez lesz a zöld vonal másik végpontja)
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= M; ++j){
			m_vb_shape.AddIndex(i + j*(N + 1), (N + 1)*(M + 1) + i + j*(N + 1));
		}
	}

	m_vb_shape.InitBuffers();

	//
	// shaderek betöltése
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


	// részecskéket kirajzoló program
	m_particleProgram.Init({	// shaderek felsorolása
		{ GL_VERTEX_SHADER,		"particle.vert" },
		{ GL_FRAGMENT_SHADER,	"particle.frag" }
		},
	{	// binding-ok felsorolása
		{ 0, "vs_in_pos" },
		{ 1, "vs_in_vel" },
	});

	// részecskék inicializálása
	m_particlePos.reserve(m_particleCount);
	m_particleVel.reserve(m_particleCount);

	// véletlenszám generátor inicializálása
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> rnd(-20, 20);
	std::uniform_real_distribution<> rndVelocity(5, 10);

	// CPU oldali részecsketömbök feltöltése
	for (int i = 0; i < m_particleCount; ++i)
	{
		m_particlePos.push_back(glm::vec3(rnd(gen), 10, rnd(gen)));
		m_particleVel.push_back(glm::vec3(0, 2 * -abs(rndVelocity(gen)), 0));
	}

	// GPU-ra áttölteni a részecskék pozícióit
	m_gpuParticleBuffer.BufferData(m_particlePos);	// <=>	m_gpuParticleBuffer = m_particlePos;

													// és végül a VAO-t inicializálni
	m_gpuParticleVAO.Init({
		{ CreateAttribute<0, glm::vec3, 0, sizeof(glm::vec3)>, m_gpuParticleBuffer }
		});

	//
	// egyéb inicializálás
	//

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);

	// textúrák betöltése
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
	// kameramozgatáshoz
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time)/1000.0f;
	m_camera.Update(delta_time);
	glm::vec3 p1 = glm::vec3(-12.0f, 20.0f, 0.0f);
	glm::vec3 p2 = glm::vec3(0.0f, 20.0f, -15.0f);
	glm::vec3 p3 = glm::vec3(12.0f, 20.0f, 0.0f);

	// Idõ € [0,1]
	if (count) {
		currentTime += SDL_GetTicks() - lastTickTime;
	}
	float t = abs(sinf(currentTime / 1000.0f / 6 * M_PI));

	// Bezier görbe
	float x = ((1.0f - t)*((1.0f - t)*p1.x + t * p2.x) + t * ((1.0f - t) * p2.x + t * p3.x));
	float y = ((1.0f - t)*((1.0f - t)*p1.y + t * p2.y) + t * ((1.0f - t) * p2.y + t * p3.y));
	float z = ((1.0f - t)*((1.0f - t)*p1.z + t * p2.z) + t * ((1.0f - t) * p2.z + t * p3.z));

	m_camera.SetView(glm::vec3(x, y, z),		// honnan nézzük a színteret
		glm::vec3(0, 0, 0),		// a színtér melyik pontját nézzük
		glm::vec3(0, 1, 0));

	// frissítsük a pozíciókat
	static const float energyRemaining = 1;	// tökéletesen rugalmas ütközés
	for (int i = 0; i < m_particleCount; ++i)
	{
		m_particlePos[i] += m_particleVel[i] * delta_time;

		if (m_particlePos[i].y <= 0) {
			m_particlePos[i].y = 20;
		}

	}

	// frissítsük a puffert
	glBindBuffer(GL_ARRAY_BUFFER, m_gpuParticleBuffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3)*m_particlePos.size(), &(m_particlePos[0][0]));
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	last_time = SDL_GetTicks();
	lastTickTime = SDL_GetTicks();
}

// talaj kirajzolása
void CMyApp::DrawGround(){

	glm::mat4 matWorld = glm::mat4(1.0f); // nem változtatunk rajta
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

	// kétszeresére nagyítjuk és feljebb toljuk, ne lógjon a talajba
	glm::mat4 matWorld = glm::translate(glm::vec3(0.0f, 13.0f, 0.0f)) * glm::scale<float>(glm::vec3(2.0f, 2.0f, 2.0f));
	glm::mat4 matWorldIT = glm::transpose(glm::inverse(matWorld));
	//glm::mat4 mvp = m_camera.GetViewProj() *matWorld;
	glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

	m_program.SetUniform("world", matWorld);
	m_program.SetUniform("worldIT", matWorldIT);
	m_program.SetUniform("MVP", mvp);

	m_program.SetTexture("texImage", 0, m_TextureID);

	m_vb_shape.On();
	// kirajzoljuk a felületet: 0-tól a következõ K = 6*N*M indexet figyelembe véve...
	m_vb_shape.DrawIndexed(GL_TRIANGLES, 0, 2 * 3 * N * M, 0);
	// ...és a zöld vonalakat is: GL_LINES a K-adik indextõl, a következõ 2*P index figyelembe vételével,
	//		ahol P = (N+1)*(M+1) a felületet alkotó vertexek száma
	//		GL_LINES az indexek által jelölt vertexeket páronként vonalakká köti össze
	m_program.SetUniform("is_green", true); // szólunk a shadernek, hogy ne számoljon fényeket, hanem konstans zöld színnel rajzoljon
	//m_vb_shape.DrawIndexed(GL_LINES, 0, 2 * (N + 1)*(M + 1), (void*)(2 * 3 * N*M*sizeof(unsigned int)));
	m_program.SetUniform("is_green", false); // ne felejtsük el visszaállítani a megjelenítést shaderben
	m_vb_shape.Off();

}

void CMyApp::Render()
{
	// töröljük a frampuffert (GL_COLOR_BUFFER_BIT) és a mélységi Z puffert (GL_DEPTH_BUFFER_BIT)
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// megjelenítés módja
	if (is_filled){
		glPolygonMode(GL_FRONT, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT, GL_LINE);
	}

	m_program.On(); // shader program bekapcs

	// light_pos a pontszerü fényforrásunk pozíciója;
	// y = 10 magasságban kering az XZ síkon, 5 mp alatt megtéve egy 10 sugarú körpályát
	float r = 2 * f_PI * SDL_GetTicks() / 1000.0f / 5.0f;
	glm::vec3 light_pos = glm::vec3(10 * cosf(r), 10, 10 * sinf(r));
	
	// a shadernek átadjuk a rajzoláshoz szükséges paramétereket
	// hol a fény és a kamera
	//m_program.SetUniform("light_pos", light_pos);
	m_program.SetUniform("eye_pos", m_camera.GetEye());

	// használunk-e textúrát
	m_program.SetUniform("use_texture", use_texture);
	
	// kirajzoljuk az alakzatokat
	//DrawGround();
	DrawShape();

	m_program.Off(); // shader program kikapcs

	// részecskék
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

// a két paraméterbe az új ablakméret szélessége (_w) és magassága (_h) található
void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_camera.Resize(_w, _h);
}