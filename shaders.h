/* Adapted from code originally written by pknowles */

/*
NOTE: make sure to call glewInit before loading shaders

use getShader() to load, compile shaders and return a program
use glUseProgram(program) to activate it
use glUseProgram(0) to return to fixed pipeline rendering
use glDeleteProgram() to free resources
*/

#ifndef SHADERS_H
#define SHADERS_H

#if __cplusplus
extern "C" {
#endif


#define CHECK_GL_ERROR oglError(__LINE__, __FILE__)
int oglError(int line, const char* file);
unsigned int getShader(const char* vertexFile, const char* fragmentFile);


#if __cplusplus
}
#endif


#endif
