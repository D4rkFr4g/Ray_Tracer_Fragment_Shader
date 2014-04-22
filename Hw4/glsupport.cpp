#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "glsupport.h"

using namespace std;

void checkGlErrors() {
  const GLenum errCode = glGetError();

  if (errCode != GL_NO_ERROR) {
    string error("GL Error: ");
    error += errCode;
    cerr << error << endl;
    throw runtime_error(error);
  }
}

// Dump text file into a character vector, throws exception on error
static void readTextFile(const char *fn, vector<char>& data) {
  // Sets ios::binary bit to prevent end of line translation, so that the
  // number of bytes we read equals file size
  ifstream ifs(fn, ios::binary);
  if (!ifs)
    throw runtime_error(string("Cannot open file ") + fn);

  // Sets bits to report IO error using exception
  ifs.exceptions(ios::eofbit | ios::failbit | ios::badbit);
  ifs.seekg(0, ios::end);
  size_t len = ifs.tellg();

  data.resize(len);

  ifs.seekg(0, ios::beg);
  ifs.read(&data[0], len);
}

// Print info regarding an GL object
static void printInfoLog(GLuint obj, const string& filename) {
  GLint infologLength = 0;
  GLint charsWritten  = 0;
  glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infologLength);
  if (infologLength > 0) {
    string infoLog(infologLength, ' ');
    glGetInfoLogARB(obj, infologLength, &charsWritten, &infoLog[0]);
    std::cerr << "##### Log [" << filename << "]:\n" << infoLog << endl;
  }
}

const char* modifyFS(const char* fs, int numLights, int numGeom)
{
   /*	PURPOSE:		Find and replace array sizes
      RECEIVES:	fs - Fragment Shader
      RETURNS:
      REMARKS:		
      */

   string frag(fs);

   string one = frag.substr(0,32);
   string two = frag.substr(34,48);
   string three = frag.substr(84);

   stringstream ss;
   ss << numLights;
   string str = ss.str();
   
   one += "" + str;

   ss << numGeom;
   str = ss.str();
   two += "" + str;

   frag = one + two + three;

   // Convert back to const char*
   char* S = new char[frag.length() + 1];
   std::strcpy(S,frag.c_str());

   return (const char*) S;

}

void readAndCompileSingleShader(GLuint shaderHandle, const char *fn) {
  vector<char> source;
  readTextFile(fn, source);

  const char *ptrs[] = {&source[0]};
  
  const GLint lens[] = {(GLint)source.size()};
  glShaderSource(shaderHandle, 1, ptrs, lens);   // load the shader sources

  glCompileShader(shaderHandle);

  printInfoLog(shaderHandle, fn);

  GLint compiled = 0;
  glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &compiled);
  if (!compiled)
    throw runtime_error("fails to compile GL shader");
}

void readAndCompileSingleShader(GLuint shaderHandle, const char*fn, int numLights, int numGeom)
{
  vector<char> source;
  readTextFile(fn, source);

  const char *ptrs[] = {&source[0]};

  const char *filename = "./Shaders/ray-tracer-gl3.fshader";
  string f1 = filename;
  string f2 = fn;
  if (f1.compare(f2) == 0)
   const char* mfn = modifyFS(*ptrs, numLights, numGeom);

  readAndCompileSingleShader(shaderHandle, fn);
}
void linkShader(GLuint programHandle, GLuint vs, GLuint fs) {
  glAttachShader(programHandle, vs);
  glAttachShader(programHandle, fs);

  glLinkProgram(programHandle);

  glDetachShader(programHandle, vs);
  glDetachShader(programHandle, fs);

  GLint linked = 0;
  glGetProgramiv(programHandle, GL_LINK_STATUS, &linked);
  printInfoLog(programHandle, "linking");

  if (!linked)
    throw runtime_error("fails to link shaders");
}


void readAndCompileShader(GLuint programHandle, const char * vertexShaderFileName, const char * fragmentShaderFileName) {
  GlShader vs(GL_VERTEX_SHADER);
  GlShader fs(GL_FRAGMENT_SHADER);

  readAndCompileSingleShader(vs, vertexShaderFileName);
  readAndCompileSingleShader(fs, fragmentShaderFileName);

  linkShader(programHandle, vs, fs);
}

void readAndCompileShader(GLuint programHandle, const char * vertexShaderFileName, const char * fragmentShaderFileName, int numLights, int numGeom)
{
   GlShader vs(GL_VERTEX_SHADER);
   GlShader fs(GL_FRAGMENT_SHADER);

   readAndCompileSingleShader(vs, vertexShaderFileName);
   readAndCompileSingleShader(fs, fragmentShaderFileName, numLights, numGeom);

   linkShader(programHandle, vs, fs);
}