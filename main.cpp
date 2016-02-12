#include <stdio.h>
#include <windows.h>
#include <stdarg.h>
#include <ctype.h>
#include "GL/glus.h"	 
#include "camera.h"
#include "LoadTexture.h" 
#include "Cube.h"
#include "AntTweakBar.h"
#include <cmath>
#include <vector>
using namespace std;

#define N 512	 

#define LENGTH	2000 	 
#define GRAVITY 981.0f 

bool g_bSkirt = true;
float AMPLITUDE = 0.5e-7f;
float WIND_SPEED = 600.0f;
GLfloat WIND_DIRECTION[2] = { 0.6f, 0.8f };
GLfloat oceanColor[4] = { 0.039f, 0.15f, 0.47f, 0.1f };
GLfloat lightDirection[3] = { 1.0f, 1.0f, 1.0f };
TwBar *bar;
Camera camera;
MeshData cube;
static const std::string cubemap_dir = "cubemap";
int nWidth = 1280;
int nHeight = 720;

GLuint cubemap_id = -1;
static GLUSprogram g_computeUpdateHtProgram; 
static GLUSprogram g_computeUpdateNormalProgram;
static GLUSprogram g_computeFftProgram;	
static GLUSprogram g_CubeProgram;
static GLint g_processColumnFftLocation;	 
static GLint g_stepsFftLocation;			 
static GLint g_totalTimeUpdateHtLocation;	

static GLuint g_textureH0;
static GLuint g_textureHt;
static GLuint g_textureDx;
static GLuint g_textureDy;
static GLuint g_textureOmg;
static GLuint g_textureIndices;
static GLuint g_textureDisplacement[6];
static GLuint g_textureNormal;

//
static bool g_bChoppy = false;
static GLint g_ChoppyLocation;
static GLUSprogram g_program;		
static GLint g_modelViewProjectionMatrixLocation;	 
static GLint g_normalMatrixLocation;		 
static GLint g_lightDirectionLocation;		 
static GLint g_colorLocation;			 
static GLint g_vertexLocation;			 
static GLint g_texCoordLocation;	   
static GLuint g_verticesVBO;  
static GLuint g_texCoordsVBO;  
static GLuint g_indicesVBO;	 
static GLuint g_vao;

//

static GLuint g_numberIndices;

static GLfloat Phillips(GLfloat waveDirection[2], GLfloat windDirection[2],GLfloat A)
{
	GLfloat k = glusVector2Dotf(waveDirection,waveDirection);
	GLfloat waveDotWind = glusVector2Dotf(waveDirection, windDirection);	

	float L = WIND_SPEED * WIND_SPEED / GRAVITY;
	float w = L / 1000;	
	float phillips = A * expf(-1.0f / (k  * L * L)) / (k * k * k ) * waveDotWind * waveDotWind;
	float dir_depend = 0.07f;

	if (waveDotWind < 0)
		phillips *= dir_depend;

	return phillips * expf(-k  * w * w) ; 
}
float Gauss()
{
	float u1 = rand() / (float)RAND_MAX;
	float u2 = rand() / (float)RAND_MAX;
	if (u1 < 1e-6f)
		u1 = 1e-6f;
	return sqrtf(-2 * logf(u1)) * cosf(2*GLUS_PI * u2);
}

void CreateH0Data()
{
	GLfloat* h0Data = (GLfloat*)malloc(N * N * 2 * sizeof(GLfloat));
	if (!h0Data)
	{
		return ;
	}
	GLfloat waveDirection[2];
	glusVector2Normalizef(WIND_DIRECTION);
	srand(0);
	for (int i = 0; i < N; i++)
	{
		// Positive N, that it matches with OpenGL z-axis.
		waveDirection[1] = ((GLfloat)-N / 2.0f + (GLfloat)i) * (2.0f * GLUS_PI / LENGTH);

		for (int j = 0; j < N; j++)
		{
			waveDirection[0] = ((GLfloat)-N / 2.0f + (GLfloat)j) * (2.0f * GLUS_PI / LENGTH);

			GLfloat phillipsSpectrum = (waveDirection[0] == 0 && waveDirection[1] == 0) ? 0 : (sqrtf(Phillips(waveDirection, WIND_DIRECTION, AMPLITUDE)));

			h0Data[i * 2 * N + j * 2 + 0] = 1 / sqrtf(2.0f) * Gauss() * phillipsSpectrum;
			h0Data[i * 2 * N + j * 2 + 1] = 1 / sqrtf(2.0f) * Gauss() * phillipsSpectrum;
		}
	}	

	glGenTextures(1, &g_textureH0);
	glBindTexture(GL_TEXTURE_2D, g_textureH0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, h0Data);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	free(h0Data);
}


void TW_CALL getFloat(void *value, void *clientData)
{
	*((float*)value) = *((float*)clientData);
}

void TW_CALL setFloat(const void *value, void *clientData)
{
	*((float*)clientData) = *((float*)value);
	CreateH0Data();
}

void TW_CALL getInt(void *value, void *clientData)
{
	*((int*)value) = *((int*)clientData);
}

void TW_CALL setInt(const void *value, void *clientData)
{
	*((int*)clientData) = *((int*)value);
	CreateH0Data();
}

void TW_CALL getBool(void *value, void *clientData)
{
	*((bool*)value) = *((bool*)clientData);
}

void TW_CALL setBool(const void *value, void *clientData)
{
	*((bool*)clientData) = *((bool*)value);
	if (g_ChoppyLocation == -1) return;
	
	glUseProgram(g_program.program);
	if (*(bool*)value > 0.5)	
		glUniform1f(g_ChoppyLocation, 1);	
	else
		glUniform1f(g_ChoppyLocation, 0);	
}
void TW_CALL setRender(const void *value, void *clientData)
{
	*((float*)clientData) = *((float*)value);
	glUseProgram(g_program.program);
	glUniform3fv(g_lightDirectionLocation, 1, lightDirection);
	glUniform4fv(g_colorLocation, 1, oceanColor);
}


GLUSboolean init(GLUSvoid)
{
	GLUStextfile computeSource;
	GLUStextfile vertexSource;
	GLUStextfile fragmentSource;
	GLUStextfile CubeSourceVS;
	GLUStextfile CubeSourceFS;
	GLUSshape gridPlane;
	GLint i;
	GLfloat matrix[16];
	
	GLint* butterflyIndices;
	GLfloat* butterflyIndicesAsFloat;
	GLfloat waveDirection[2];	
	
	camera.SetPosition(glm::vec3(0.0f, 100.0f, 0.0f));
	camera.SetLookAt(glm::vec3(0, 100, 100));
	camera.SetClipping(.1, 30000);
	camera.SetViewport(0, 0, 1280, 720);
	camera.SetFOV(45);

	//

	GLint steps = 0;
	GLint temp = N;

	while (!(temp & 0x1))
	{
		temp = temp >> 1;
		steps++;
	}

	//

	glusFileLoadText("shader/ocean_update_ht.comp.glsl", &computeSource);
	glusProgramBuildComputeFromSource(&g_computeUpdateHtProgram, (const GLchar**) &computeSource.text);
	glusFileDestroyText(&computeSource);

	glusFileLoadText("shader/ocean_fft.comp.glsl", &computeSource);
	glusProgramBuildComputeFromSource(&g_computeFftProgram, (const GLchar**) &computeSource.text);
	glusFileDestroyText(&computeSource);

	glusFileLoadText("shader/ocean_update_normal.comp.glsl", &computeSource);
	glusProgramBuildComputeFromSource(&g_computeUpdateNormalProgram, (const GLchar**) &computeSource.text);
	glusFileDestroyText(&computeSource);
	
	glusFileLoadText("shader/ocean.vert.glsl", &vertexSource);
	glusFileLoadText("shader/ocean.frag.glsl", &fragmentSource);

	GLUSboolean b1 = glusProgramBuildFromSource(&g_program, (const GLUSchar**)&vertexSource.text, 0, 0, 0, (const GLUSchar**)&fragmentSource.text);

	glusFileDestroyText(&vertexSource);
	glusFileDestroyText(&fragmentSource);

	glusFileLoadText("shader/cube_vs.glsl", &CubeSourceVS);
	glusFileLoadText("shader/cube_fs.glsl", &CubeSourceFS);
	GLUSboolean b2 = glusProgramBuildFromSource(&g_CubeProgram, (const GLchar**)&CubeSourceVS.text,0,0,0,(const GLchar**)&CubeSourceFS.text);
	glusFileDestroyText(&CubeSourceVS);
	glusFileDestroyText(&CubeSourceFS);
	
	if (!b1 || !b2)
	{
		int a = 0;
	}
	glUseProgram(g_CubeProgram.program);
	cube = CreateCube();
	cubemap_id = LoadCube(cubemap_dir);


	g_totalTimeUpdateHtLocation = glGetUniformLocation(g_computeUpdateHtProgram.program, "u_totalTime");
	g_processColumnFftLocation = glGetUniformLocation(g_computeFftProgram.program, "u_processColumn");
	g_stepsFftLocation = glGetUniformLocation(g_computeFftProgram.program, "u_steps");
	g_modelViewProjectionMatrixLocation = glGetUniformLocation(g_program.program, "u_modelViewProjectionMatrix");
	g_ChoppyLocation = glGetUniformLocation(g_program.program, "bChoppy");
	g_normalMatrixLocation = glGetUniformLocation(g_program.program, "u_normalMatrix");
	g_lightDirectionLocation = glGetUniformLocation(g_program.program, "u_lightDirection");
	g_colorLocation = glGetUniformLocation(g_program.program, "u_color");
	g_vertexLocation = glGetAttribLocation(g_program.program, "a_vertex");
	g_texCoordLocation = glGetAttribLocation(g_program.program, "a_texCoord");
		
	// Use a helper function to create a grid plane.
	glusShapeCreateRectangularGridPlanef(&gridPlane, LENGTH, LENGTH, N - 1, N - 1, GLUS_FALSE);
	g_numberIndices = gridPlane.numberIndices;
	
	// Rotate by 90 degrees, that the grid is in the x-z-plane.
	glusMatrix4x4Identityf(matrix);
	glusMatrix4x4RotateRxf(matrix, -90.0f);
	for (i = 0; i < gridPlane.numberVertices; i++)
	{
		glusMatrix4x4MultiplyPoint4f(&gridPlane.vertices[4 * i], matrix, &gridPlane.vertices[4 * i]);
	}

	//´¹Ö±È¹°ÚVertical skirt
	for (i = 0; i < N && g_bSkirt; i++)
	{
		GLUSfloat* p1 = &gridPlane.vertices[4 * i];
		GLUSfloat* p2 = &gridPlane.vertices[4 * (gridPlane.numberVertices - i)];
		*(p1 + 1) -= 100;
		*(p2 + 1) -= 100;

		if (i != 0)
		{
			GLUSfloat* p3 = &gridPlane.vertices[4 * i*N];
			*(p3 + 1) -= 100;
		}	

		if (i != 0 && i != N-1)
		{
			GLUSfloat* p4 = &gridPlane.vertices[4 * ((i+1)*N-1) ];
			*(p4 + 1) -= 100;
		}
	}


	glGenBuffers(1, &g_verticesVBO);
	glBindBuffer(GL_ARRAY_BUFFER, g_verticesVBO);
	glBufferData(GL_ARRAY_BUFFER, gridPlane.numberVertices * 4 * sizeof(GLfloat), (GLfloat*) gridPlane.vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &g_texCoordsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, g_texCoordsVBO);
	glBufferData(GL_ARRAY_BUFFER, gridPlane.numberVertices * 2 * sizeof(GLfloat), (GLfloat*) gridPlane.texCoords, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &g_indicesVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_indicesVBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, gridPlane.numberIndices * sizeof(GLuint), (GLuint*) gridPlane.indices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	glusShapeDestroyf(&gridPlane);
	
	
	// Generate H0
	CreateH0Data();

	glGenTextures(1, &g_textureHt);
	glBindTexture(GL_TEXTURE_2D, g_textureHt);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDx);
	glBindTexture(GL_TEXTURE_2D, g_textureDx);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDy);
	glBindTexture(GL_TEXTURE_2D, g_textureDy);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	glGenTextures(1, &g_textureDisplacement[0]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	glGenTextures(1, &g_textureDisplacement[1]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDisplacement[2]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[2]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDisplacement[3]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[3]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDisplacement[4]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[4]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &g_textureDisplacement[5]);
	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[5]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	glGenTextures(1, &g_textureNormal);
	glBindTexture(GL_TEXTURE_2D, g_textureNormal);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, N, N, 0, GL_RGBA, GL_FLOAT, 0);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	glBindTexture(GL_TEXTURE_2D, 0);

	//
	// Generate indices needed for inverse FFT butterfly algorithm.
	//

	butterflyIndices = (GLint*)malloc(N * sizeof(GLint));

	if (!butterflyIndices)
	{
		return GLUS_FALSE;
	}

	butterflyIndicesAsFloat = (GLfloat*)malloc(N * sizeof(GLfloat));

	if (!butterflyIndicesAsFloat)
	{
		free(butterflyIndices);

		return GLUS_FALSE;
	}

	for (i = 0; i < N; i++)
	{
		butterflyIndices[i] = i;
	}

	glusFourierButterflyShuffleFFTi(butterflyIndices, butterflyIndices, N);

	for (i = 0; i < N; i++)
	{
		butterflyIndicesAsFloat[i] = (GLfloat)butterflyIndices[i];
	}

	free(butterflyIndices);

	glGenTextures(1, &g_textureIndices);
	glBindTexture(GL_TEXTURE_1D, g_textureIndices);

	
	glTexImage1D(GL_TEXTURE_1D, 0, GL_R32F, N, 0, GL_RED, GL_FLOAT, butterflyIndicesAsFloat);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glBindTexture(GL_TEXTURE_1D, 0);
	free(butterflyIndicesAsFloat);

	//


	glUseProgram(g_computeFftProgram.program);
	glUniform1i(g_stepsFftLocation, steps);

	//

	glUseProgram(g_program.program);

	glUniform1f(g_ChoppyLocation, 0);

	glGenVertexArrays(1, &g_vao);
	glBindVertexArray(g_vao);

	glBindBuffer(GL_ARRAY_BUFFER, g_verticesVBO);
	glVertexAttribPointer(g_vertexLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(g_vertexLocation);

	glBindBuffer(GL_ARRAY_BUFFER, g_texCoordsVBO);
	glVertexAttribPointer(g_texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(g_texCoordLocation);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_indicesVBO);
	glBindVertexArray(0);

	//

	glusVector3Normalizef(lightDirection);
	glUniform3fv(g_lightDirectionLocation, 1, lightDirection);
	glUniform4fv(g_colorLocation, 1, oceanColor);
	glUseProgram(0);

	//

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	TwDefine(" GLOBAL fontscaling=1.3 ");
	TwInit(TW_OPENGL_CORE, NULL);
	bar = TwNewBar("Parameters");
	TwDefine("Parameters size='300 300'");
	TwAddVarCB(bar, "Wind speed", TW_TYPE_FLOAT, setFloat, getFloat, &WIND_SPEED, "min=100.0 max=5000.0 step=500.0 group=Spectrum");
	TwAddVarCB(bar, "Wind direction x", TW_TYPE_FLOAT, setFloat, getFloat, &WIND_DIRECTION[0], "min=0.0 max=1.0 step=0.2 group=Spectrum");
	TwAddVarCB(bar, "Wind direction y", TW_TYPE_FLOAT, setFloat, getFloat, &WIND_DIRECTION[1], "min=0.0 max=1.0 step=0.2 group=Spectrum");
	TwAddVarCB(bar, "Amplitude", TW_TYPE_FLOAT, setFloat, getFloat, &AMPLITUDE, "min=0.5e-8f max=0.5e-5f step=0.8e-7f group=Spectrum");
	TwAddVarCB(bar, "Choppy", TW_TYPE_BOOL8, setBool, getBool, &g_bChoppy, "group=Rendering");
	TwAddVarCB(bar, "OceanColor R", TW_TYPE_FLOAT, setRender, getFloat,&oceanColor[0],  "min = 0.0 max = 1.0 step = 0.2 group = Rendering");
	TwAddVarCB(bar, "OceanColor G", TW_TYPE_FLOAT, setRender, getFloat, &oceanColor[1], "min = 0.0 max = 1.0 step = 0.2 group=Rendering");
	TwAddVarCB(bar, "OceanColor B", TW_TYPE_FLOAT, setRender, getFloat, &oceanColor[2], "min = 0.0 max = 1.0 step = 0.2 group=Rendering");
	TwAddVarCB(bar, "lightDirectionX", TW_TYPE_FLOAT, setRender, getFloat, &lightDirection[0], "min = 0.0 max = 1.0 step = 0.2 group=Rendering");
	TwAddVarCB(bar, "lightDirectionY", TW_TYPE_FLOAT, setRender, getFloat, &lightDirection[1], "min = 0.0 max = 1.0 step = 0.2 group=Rendering");
	TwAddVarCB(bar, "lightDirectionZ", TW_TYPE_FLOAT, setRender, getFloat, &lightDirection[2], "min = 0.0 max = 1.0 step = 0.2 group=Rendering");
	TwAddVarCB(bar, "Skirt", TW_TYPE_BOOL8, setBool, getBool, &g_bSkirt, "group=Rendering");
	
	return GLUS_TRUE;
}

GLUSvoid reshape(GLUSint width, GLUSint height)
{
	camera.SetViewport(0, 0, width, height);
	nWidth = width;
	nHeight = height;
	TwWindowSize(width, height);
}
void drawSkybox(const glm::mat4& PV)
{
	glm::mat4 Msky = glm::scale(glm::vec3(30000.0f));
	glm::mat4 PVsky = PV;
	PVsky[3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	
	glUseProgram(g_CubeProgram.program);
	int PVM_loc = glGetUniformLocation(g_CubeProgram.program, "PVM");
	if (PVM_loc != -1)
	{
		glm::mat4 PVM = PVsky*Msky;
		glUniformMatrix4fv(PVM_loc, 1, false, glm::value_ptr(PVM));
	}

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap_id);
	glBindVertexArray(cube.mVao);
	glDrawElements(GL_TRIANGLES, cube.mNumIndices, GL_UNSIGNED_SHORT, 0);
}

GLUSboolean display(GLUSfloat time)
{
	static GLfloat totalTime = 0.0f;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	
	glm::mat4 M, V, P;
	camera.Update();
	camera.GetMatricies(P, V, M);
	glm::mat4 PV = P*V;

	glm::mat3 Normal = glm::mat3(1.0f);
	glUniformMatrix3fv(g_normalMatrixLocation, 1, GL_FALSE, glm::value_ptr(Normal));
	drawSkybox(PV);
	
	glUseProgram(g_program.program);	
	if (g_modelViewProjectionMatrixLocation != -1)
	{
		glm::mat4 PVM = PV*M;
		glUniformMatrix4fv(g_modelViewProjectionMatrixLocation, 1, false, glm::value_ptr(PVM));
	}

	//
	// Simulation part.
	//

	//
	// Update H pass.
	//

	glUseProgram(g_computeUpdateHtProgram.program);

	glBindImageTexture(0, g_textureH0, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureHt, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
	glBindImageTexture(2, g_textureDx, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
	glBindImageTexture(3, g_textureDy, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
	

	glUniform1f(g_totalTimeUpdateHtLocation, totalTime);

	// Process all vertices. No synchronization needed, so start NxN threads with local size of 1x1.
	glDispatchCompute(N, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	//  pass.
	//

	glUseProgram(g_computeFftProgram.program);

	glBindImageTexture(2, g_textureIndices, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);

	//
	// per row pass for ht.
	//

	glBindImageTexture(0, g_textureHt, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[0], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 0);

	// Process all vertices. N groups as N rows are processed. One work group is one row.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// FFT per column pass for ht.
	//

	glBindImageTexture(0, g_textureDisplacement[0], 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[1], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 1);

	// Process all vertices. N groups as N columns are processed. One work group is one column.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// per row pass for dx.
	//

	glBindImageTexture(0, g_textureDx, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[2], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 0);

	// Process all vertices. N groups as N rows are processed. One work group is one row.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// FFT per column pass for dx.
	//

	glBindImageTexture(0, g_textureDisplacement[2], 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[3], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 1);

	// Process all vertices. N groups as N columns are processed. One work group is one column.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// per row pass for dy.
	//

	glBindImageTexture(0, g_textureDy, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[4], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 0);

	// Process all vertices. N groups as N rows are processed. One work group is one row.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// FFT per column pass for dy.
	//

	glBindImageTexture(0, g_textureDisplacement[4], 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureDisplacement[5], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

	glUniform1i(g_processColumnFftLocation, 1);

	// Process all vertices. N groups as N columns are processed. One work group is one column.
	glDispatchCompute(1, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// Update normal map pass.
	//

	glUseProgram(g_computeUpdateNormalProgram.program);

	glBindImageTexture(0, g_textureDisplacement[1], 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	glBindImageTexture(1, g_textureNormal, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	// Process all vertices. No synchronization needed, so start NxN threads with local size of 1x1.
	glDispatchCompute(N, N, 1);

	// Make sure, all values are written.
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	//
	// Drawing part.
	//

	glUseProgram(g_program.program);

	glBindVertexArray(g_vao);

	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[1]);
	glActiveTexture(GL_TEXTURE1);

  	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[3]);
 	glActiveTexture(GL_TEXTURE2);
 
  	glBindTexture(GL_TEXTURE_2D, g_textureDisplacement[5]);
 	glActiveTexture(GL_TEXTURE3);	

	glBindTexture(GL_TEXTURE_2D, g_textureNormal);
	glActiveTexture(GL_TEXTURE4);
			
//	glDrawElements(GL_TRIANGLES, g_numberIndices, GL_UNSIGNED_INT, 0);
	glDrawElementsInstanced(GL_TRIANGLES, g_numberIndices, GL_UNSIGNED_INT, 0, 49);

	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindVertexArray(0);

	totalTime += time;
	TwDraw();
	return GLUS_TRUE;
}

GLUSvoid terminate(GLUSvoid)
{
	glBindTexture(GL_TEXTURE_1D, 0);
	if (g_textureIndices)
	{
		glDeleteTextures(1, &g_textureIndices);
		g_textureIndices = 0;
	}

	glBindTexture(GL_TEXTURE_2D, 0);
	if (g_textureH0)
	{
		glDeleteTextures(1, &g_textureH0);
		g_textureH0 = 0;
	}

	if (g_textureHt)
	{
		glDeleteTextures(1, &g_textureHt);
		g_textureHt = 0;
	}

	if (g_textureDisplacement[0])
	{
		glDeleteTextures(1, &g_textureDisplacement[0]);
		g_textureDisplacement[0] = 0;
	}

	if (g_textureDisplacement[1])
	{
		glDeleteTextures(1, &g_textureDisplacement[1]);
		g_textureDisplacement[1] = 0;
	}

	if (g_textureNormal)
	{
		glDeleteTextures(1, &g_textureNormal);
		g_textureNormal = 0;
	}

	//

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	if (g_verticesVBO)
	{
		glDeleteBuffers(1, &g_verticesVBO);
		g_verticesVBO = 0;
	}

	if (g_texCoordsVBO)
	{
		glDeleteBuffers(1, &g_texCoordsVBO);
		g_texCoordsVBO = 0;
	}

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	if (g_indicesVBO)
	{
		glDeleteBuffers(1, &g_indicesVBO);
		g_indicesVBO = 0;
	}

	glBindVertexArray(0);
	if (g_vao)
	{
		glDeleteVertexArrays(1, &g_vao);

		g_vao = 0;
	}

	//
	glUseProgram(0);
	glusProgramDestroy(&g_computeUpdateHtProgram);
	glusProgramDestroy(&g_computeFftProgram);
	glusProgramDestroy(&g_computeUpdateNormalProgram);
	glusProgramDestroy(&g_program);
}

void CallBackKeyboardFunc(GLUSboolean pressed, GLUSint key)
{
	if (!pressed) return;
	
	
	if (TwEventKeyboardGLUT(key, 0, 0))
	{
		FLOAT ww = WIND_SPEED;
		return;
	}
	switch (key)
	{	
	case 'w':
		camera.Move(FORWARD);
		break;
	case 'a':
		camera.Move(LEFT);
		break;
	case 's':
		camera.Move(BACK);
		break;
	case 'd':
		camera.Move(RIGHT);
		break;
	case 'q':
		camera.Move(DOWN);
		break;
	case 'e':
		camera.Move(UP);
		break;
	}
}

int oldx;
int oldy;
bool drag;
void CallBackMouseFunc(GLUSboolean pressed, GLUSint button, GLUSint x, GLUSint y)
{
	camera.SetPos(button, pressed, x, y);
	drag = false;
	if (!TwEventMouseButtonGLUT(!button, !pressed, x, y) && button == 1) {
		oldx = x;
		oldy = y;
		drag = true;
	}
}
void CallBackMotionFunc(GLUSint button, GLUSint x, GLUSint y)
{	
	if (drag) {		
		oldx = x;
		oldy = y;
		camera.Move2D(x, y);
	}
	else {
		TwMouseMotion(x, y);	
	}
}


int main(int argc, char* argv[])
{
	EGLint eglConfigAttributes[] = {
		EGL_RED_SIZE, 8,
		EGL_GREEN_SIZE, 8,
		EGL_BLUE_SIZE, 8,
		EGL_DEPTH_SIZE, 24,
		EGL_STENCIL_SIZE, 0,
		EGL_SAMPLE_BUFFERS, 1,
		EGL_SAMPLES, 4,
		EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
		EGL_NONE
	};

	EGLint eglContextAttributes[] = {
		EGL_CONTEXT_MAJOR_VERSION, 4,
		EGL_CONTEXT_MINOR_VERSION, 3,
		EGL_CONTEXT_OPENGL_FORWARD_COMPATIBLE, EGL_TRUE,
		EGL_CONTEXT_OPENGL_PROFILE_MASK, EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT,
		EGL_CONTEXT_OPENGL_DEBUG, EGL_TRUE,
		EGL_NONE
	};

	glusWindowSetInitFunc(init);
	glusWindowSetReshapeFunc(reshape);	
	glusWindowSetUpdateFunc(display);
	glusWindowSetTerminateFunc(terminate);
	glusWindowSetMouseFunc(CallBackMouseFunc);
	
	glusWindowSetMouseMoveFunc(CallBackMotionFunc);
	glusWindowSetKeyFunc(CallBackKeyboardFunc);

	if (!glusWindowCreate("Jianping Tang - Final project", nWidth, nHeight, GLUS_FALSE, GLUS_FALSE, eglConfigAttributes, eglContextAttributes, 0))
	{
		printf("Could not create window!\n");
		return -1;
	}

	glusWindowRun();

	return 0;
}
