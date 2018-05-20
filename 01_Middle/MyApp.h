#pragma once

// GLEW
#include <GL/glew.h>

// SDL
#include <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>

#include "gCamera.h"
#include "gShaderProgram.h"
#include "gVertexBuffer.h"
#include "Mesh_OGL3.h"

#include "ProgramObject.h"
#include "BufferObject.h"
#include "VertexArrayObject.h"

class CMyApp
{
public:
	CMyApp(void);
	~CMyApp(void);

	bool Init();
	void Clean();

	void Update();
	void Render();

	void KeyboardDown(SDL_KeyboardEvent&);
	void KeyboardUp(SDL_KeyboardEvent&);
	void MouseMove(SDL_MouseMotionEvent&);
	void MouseDown(SDL_MouseButtonEvent&);
	void MouseUp(SDL_MouseButtonEvent&);
	void MouseWheel(SDL_MouseWheelEvent&);
	void Resize(int, int);
protected:
	// belsõ eljárások

	// paremtrikus egyenletek
	glm::vec3 GetUV(float u, float v);
	// a hozzájuk tartozó normálvektorok
	glm::vec3 GetUV_normals(float u, float v);

	// a rajzoló metódusaink
	void DrawShape();
	void DrawGround();

	// Textúra erõforrás azonosítók
	GLuint m_TextureID;

	gCamera			m_camera;
	gShaderProgram	m_program;
	glm::mat4 m_matView;

	gVertexBuffer	m_vb_ground;
	gVertexBuffer	m_vb_shape;

	int N = 75; // az alakzat felbontása
	int M = 75;

	bool is_filled = true; // kitöltött vagy drótvázas megjelenítés?
	bool use_texture = false; // textúrával vagy anélkül?

	// Melyik alakzatot szeretnénk kirajzolni?
	// 1 = TÓRUSZ
	// 2 = DOMBORZAT
	// 3 = VÁZA
	// 4 = TREFOIL CSOMÓ
	// 5 = ZIZI-GÖMB
	// 6 = HULLÁMFELÜLET
	int shape_id = 7;


	ProgramObject		m_particleProgram;
	int		m_particleCount = 1000;

	std::vector<glm::vec3>		m_particlePos{};
	std::vector<glm::vec3>		m_particleVel{};

	VertexArrayObject	m_gpuParticleVAO;
	ArrayBuffer			m_gpuParticleBuffer;

	bool count = true;
	float currentTime = 0;
	float lastTickTime = 0;

};

