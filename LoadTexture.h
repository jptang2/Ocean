#ifndef __LOADTEXTURE_H__
#define __LOADTEXTURE_H__

#include <windows.h>
#include <string>
#include "GL/glew.h"
#include "GL/gl.h"


GLuint LoadTexture(const std::string& fname);
GLuint LoadCube(const std::string& fname);


#endif