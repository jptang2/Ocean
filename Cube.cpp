#include <glm/glm.hpp>
#include <GL/glew.h>

#include "Cube.h"



MeshData CreateCube()
{
   MeshData meshdata;

   meshdata.mScaleFactor = 0.5f;

   meshdata.mVboNormals = 0;

   using namespace glm;
   vec3 pos[8] = {vec3(-1.0f, -1.0f, -1.0f), vec3(+1.0f, -1.0f, -1.0f), vec3(+1.0f, +1.0f, -1.0f), vec3(-1.0f, +1.0f, -1.0f),
   vec3(-1.0f, -1.0f, +1.0f), vec3(+1.0f, -1.0f, +1.0f), vec3(+1.0f, +1.0f, +1.0f), vec3(-1.0f, +1.0f, +1.0f)};

   vec3 tex[8] = {vec3(-1.0f, -1.0f, -1.0f), vec3(+1.0f, -1.0f, -1.0f), vec3(+1.0f, +1.0f, -1.0f), vec3(-1.0f, +1.0f, -1.0f),
   vec3(-1.0f, -1.0f, +1.0f), vec3(+1.0f, -1.0f, +1.0f), vec3(+1.0f, +1.0f, +1.0f), vec3(-1.0f, +1.0f, +1.0f)};


   unsigned short idx[36] = { 0,2,1, 2,0,3, //bottom
                              0,5,4, 5,0,1, //front
                              1,6,5, 6,1,2, //right 
                              2,7,6, 7,2,3, //back
                              3,4,7, 4,3,0, //left
                              4,5,6, 6,7,4};//top


   //shader attrib locations
   int pos_loc = -1;
   int tex_coord_loc = -1;
   int normal_loc = -1;

   GLint program = -1;
   glGetIntegerv(GL_CURRENT_PROGRAM, &program);

   pos_loc = glGetAttribLocation(program, "pos_attrib");
   tex_coord_loc = glGetAttribLocation(program, "tex_coord_attrib");
   normal_loc = glGetAttribLocation(program, "normal_attrib");

   meshdata.mNumIndices = 36;

   glGenVertexArrays(1, &meshdata.mVao);
	glBindVertexArray(meshdata.mVao);

   //Buffer indices
   glGenBuffers(1, &meshdata.mIndexBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshdata.mIndexBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(idx), idx, GL_STATIC_DRAW);

	//Buffer vertices
	glGenBuffers(1, &meshdata.mVboVerts);
	glBindBuffer(GL_ARRAY_BUFFER, meshdata.mVboVerts);
	glBufferData(GL_ARRAY_BUFFER, sizeof(pos), pos, GL_STATIC_DRAW);
	glEnableVertexAttribArray(pos_loc);
	glVertexAttribPointer(pos_loc, 3, GL_FLOAT, 0, 0, 0);
		

   // buffer for vertex texture coordinates
	glGenBuffers(1, &meshdata.mVboTexCoords);
	glBindBuffer(GL_ARRAY_BUFFER, meshdata.mVboTexCoords);
	glBufferData(GL_ARRAY_BUFFER, sizeof(tex), tex, GL_STATIC_DRAW);
	glEnableVertexAttribArray(tex_coord_loc);
	glVertexAttribPointer(tex_coord_loc, 3, GL_FLOAT, 0, 0, 0);

   // no normals
   meshdata.mVboNormals = 0;
   glDisableVertexAttribArray(normal_loc);
      
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);

   return meshdata;
}
