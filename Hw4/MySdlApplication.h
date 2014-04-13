#ifndef _MYSDLAPPLICATION_H_
#define _MYSDLAPPLICATION_H_

#ifdef WIN32
   #include<windows.h> //for windows
#endif

#include <SDL.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "glsupport.h"

using namespace std;

class MySdlApplication
{
   private:
      // Globals
      bool g_running;
      SDL_Window* g_display;

      // Functions
      void keyboard(const char * key);
      void mouse(SDL_MouseButtonEvent button);
      void motion(const int x, const int y);

   public:
      // Functions
      MySdlApplication();
      int onExecute();
      bool onInit();
      void onEvent(SDL_Event* Event);
      void keyboard();
      void onLoop();
      void onRender();
      void onCleanup();
      void initScene();
      void initScene2();
      void loadScene();
};

#endif