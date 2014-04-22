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

string writeTextFile(const char *fn, string s)
{
  string filename(fn);
  filename += "b";

  ofstream myfile;
  myfile.open (filename);
  myfile << s;
  myfile.close();

  return filename;
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

string modifyFS(const char* fn, const char* fs, int numLights, int numGeom)
{
   /*	PURPOSE:		Find and replace array sizes
      RECEIVES:	fs - Fragment Shader
      RETURNS:
      REMARKS:		
      */

   int lightStride = 3;
   int geomStride = 6;

   string frag(fs);

   string one = frag.substr(0,32);
   string two = frag.substr(34,48);
   string three = frag.substr(84,136);
   string four = frag.substr(221,67);
   string five = frag.substr(289);

   stringstream ss;
   ss.str(std::string());
   
   ss << numLights / lightStride;
   string str = ss.str();
   one += "" + str;

   ss.str(std::string());
   ss << numGeom / geomStride;
   str = ss.str();
   two += "" + str;

   ss.str(std::string());
   ss << (numLights == 0 ? 1 : numLights); // GLSL doesn't like [0]
   str = ss.str(); 
   three += "" + str;

   ss.str(std::string());
   ss << numGeom;
   str = ss.str();
   four += "" + str;

   // Remove junk at end of file
   bool isRemoving = true;
   while (isRemoving)
   {
      string test = five.substr(five.length() - 1);
      if (five.substr(five.length() - 1).compare("}") != 0)
         five = five.substr(0, five.length() - 1);
      else
         isRemoving = false;
   }

   frag = one + two + three + four + five;

   return writeTextFile(fn, frag);
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
  {
     string s = modifyFS(fn, *ptrs, numLights, numGeom);
     const char* mfn = s.c_str();
     readAndCompileSingleShader(shaderHandle, mfn);
  }
  else
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