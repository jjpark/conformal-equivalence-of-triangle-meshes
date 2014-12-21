#include "main.hpp"
#include "RgbImage.h"
GLfloat mouseX, mouseY = -3;
GLfloat cameraAngleX = 0;
GLfloat cameraAngleY = 194;
int mesh_or_par;
bool mouseLeftDown = false;
GLuint  texture[2];     // Storage For One Texture ( NEW )
int drawMode;
char* filename = "";
int tile=4677; //781
double t=0;
float zoom=-4;
float translateY=-2;
float translateX=2.51;
float rotateZ=0;
bool mouseRightDown;
bool mouseMiddleDown;
/*2.51 -2 0 194 -4
* Read a texture map from a BMP bitmap file.
*/

void loadTextureFromBmp(char *filename)
{   
   //glClearColor (0.0, 0.0, 0.0, 0.0);
   glEnable(GL_DEPTH_TEST);

   RgbImage theTexMap( filename );

   // Pixel alignment: each row is word aligned (aligned to a 4 byte boundary)
   //    Therefore, no need to call glPixelStore( GL_UNPACK_ALIGNMENT, ... );

  
   glGenTextures(1, &texture[0]);               // Create The Texture
      glBindTexture(GL_TEXTURE_2D, texture[0]);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      // Typical Texture Generation Using Data From The Bitmap
   
      glTexImage2D(GL_TEXTURE_2D, 0, 3, theTexMap.GetNumCols(), theTexMap.GetNumRows(), 0, GL_RGB, GL_UNSIGNED_BYTE, theTexMap.ImageData() );

}
GLfloat *vertices1;
GLfloat *normals1;
GLfloat *tex2;
GLfloat *tex;

#define checkImageWidth 64
#define checkImageHeight 64

static GLubyte checkImage[checkImageHeight][checkImageWidth][4];
int num_tri;

void set_vertex (GLfloat *v, int n, GLfloat *norm, GLfloat *texture,int d) {
  mesh_or_par = d;
  num_tri = n;
  vertices1 = v;
  normals1 = norm;
  tex = texture;
  tex2 = new GLfloat[num_tri * 3 * 2];
int z;
for (z = 0; z < num_tri; z++) {
  tex2[z*6+0] = 0;
  tex2[z*6+1] = 0;
  tex2[z*6+2] = 0.5;
  tex2[z*6+3] = 0;
  tex2[z*6+4] = 0.5;
  tex2[z*6+5] = 0.5;
}

}

void makeCheckImage(void)
{
   int i, j, c;
    
   for (i = 0; i < checkImageHeight; i++) 
   {
      for (j = 0; j < checkImageWidth; j++) 
      {
   c = ( ( ((i&0x4)==0) ^ ((j&0x4)==0) ) )*1;
         checkImage[i][j][0] = (GLubyte) 255-255*c;
         checkImage[i][j][1] = (GLubyte) 174-50*c;
         checkImage[i][j][2] = (GLubyte) 185;
         checkImage[i][j][3] = (GLubyte) 255;
      }
   }
}

void loadTextureFromFile(char *filename)
{   
   glClearColor (1.0, 1.0, 1.0, 0.0);
   glShadeModel(GL_SMOOTH);
   glEnable(GL_DEPTH_TEST);
   makeCheckImage();

   // Pixel alignment: each row is word aligned (aligned to a 4 byte boundary)
   //    Therefore, no need to call glPixelStore( GL_UNPACK_ALIGNMENT, ... );

  
  glGenTextures(1, &texture[0]);          // Create The Texture
    glBindTexture(GL_TEXTURE_2D, texture[0]);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    // Typical Texture Generation Using Data From The Bitmap
  
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, checkImageWidth, checkImageHeight, 
                0, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);
 

}

/*
* Draw the texture in the OpenGL graphics window
*/
void drawScene(void)
{
 
   //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   //glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, texture[0]);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

   glVertexPointer(2,GL_FLOAT, 0, vertices1);
   glNormalPointer(GL_FLOAT, 0, normals1);
   glTexCoordPointer(2,GL_FLOAT,0,tex);
   glLoadIdentity();
   glTranslatef(0.0,0.0,-5.000001);
   glTranslatef(translateX,translateY,zoom);
   //glRotatef(160,0,0,1);
   //glRotatef(160,0,0,1);

   //glRotatef(150,1,0,0);
   //glRotatef()
   glRotatef(rotateZ,0,0,1);
   glRotatef(cameraAngleX, 1, 0, 0);
   glRotatef(cameraAngleY, 0, 1, 0);
   glDrawArrays(GL_TRIANGLES, 0, num_tri * 3);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
   glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  
}

void drawScene2(void)
{
 
   //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, texture[0]);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

   glVertexPointer(3,GL_FLOAT, 0, vertices1);
   glNormalPointer(GL_FLOAT, 0, normals1);
   glTexCoordPointer(2,GL_FLOAT,0,tex);
   glLoadIdentity();
   glTranslatef(0.0,0.0,-5.000001);
   glTranslatef(translateX,translateY,zoom);
   //glRotatef(160,0,0,1);
   //glRotatef(160,0,0,1);

   //glRotatef(150,1,0,0);
   //glRotatef()
   glRotatef(rotateZ,0,0,1);
   glRotatef(cameraAngleX, 1, 0, 0);
   glRotatef(cameraAngleY, 0, 1, 0);
   glDrawArrays(GL_TRIANGLES, 0, num_tri * 3);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
   glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  
}


void resizeWindow(int x, int y)
{
 if (y == 0 || x == 0) return;  //Nothing is visible then, so return
    //Set a new projection matrix
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity();
    //Angle of view:40 degrees
    //Near clipping plane distance: 0.5
    //Far clipping plane distance: 20.0
     
    gluPerspective(40.0,(GLdouble)x/(GLdouble)y,0.5,35.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);  //Use the whole window for rendering
}

void keyboard (unsigned char key, int x, int y)
{
   switch (key) { 
      case 27:
         exit(0);
         break;
      default:
         break;
   }
}

void drawScenes(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  if (!mesh_or_par)
    drawScene();
  else
    drawScene2();
  glShadeModel(GL_SMOOTH);
  glutSwapBuffers();
}

void init(int argc, char **argv) {
   drawMode=0;
    mouseLeftDown = mouseRightDown = mouseMiddleDown = false;

   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize(750, 750);
   glutInitWindowPosition(100, 100);
   glutCreateWindow(argv[0]);
   glutReshapeFunc(resizeWindow);
   glutKeyboardFunc(keyboardCB);
   glutMouseFunc(mouseCB);
   glutMotionFunc(mouseMotionCB);
   glutDisplayFunc(drawScenes);
   glutTimerFunc(33, timerCB, 33); 

   
}

void init_material(void) {
   GLfloat mat_specular[] = { 1,1,1};
   GLfloat mat_shininess[] = { 5 };
   GLfloat mat_ambient[] = {0.2, 0.2,0.2};
   glShadeModel(GL_SMOOTH);

   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);

   glEnable(GL_DEPTH_TEST);

}

void start(int argc, char **argv) {
  init(argc, argv);
  initLights();
  init_material();

  loadTextureFromFile( filename );
   glutMainLoop();
}

void initLights() {
    glEnable(GL_LIGHTING);
    
    GLfloat amb[] = { 1,1,1, 0.0 };
    GLfloat diff[]= { 1,1,1};
    GLfloat spec[]= { 1,1,1};
    GLfloat lightpos[]= { -2,4,6 };

    glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
    //glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);
    glLightfv(GL_LIGHT0, GL_SPECULAR, spec);
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
    glEnable(GL_LIGHT0);
  


}

void timerCB(int millisec)
{
    glutTimerFunc(millisec, timerCB, millisec);
    glutPostRedisplay();
}
void mouseCB(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }

    else if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseMiddleDown = true;
        }
        else if(state == GLUT_UP)
            mouseMiddleDown = false;
    }
    else if(button == 3){
        if (state == GLUT_DOWN) {
        zoom += 0.1f;

        }
        else if(state == GLUT_UP) {
        }
    }
    else if(button == 4) {
        if (state == GLUT_DOWN) {
            zoom -= 0.1f;
        }
        else if(state == GLUT_UP) {
        }
    }

}


void mouseMotionCB(int x, int y)
{
    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        translateY -= (y - mouseY) * 0.005f;
        mouseY = y;
        translateX += (x - mouseX) * 0.005f;
        mouseX = x;

    }
}


void keyboardCB(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 27: // ESCAPE
        exit(0);
        break;

    case 'd':
        drawMode = ++drawMode % 3;
        if(drawMode == 0)        // fill mode
        {   glEnable(GL_LIGHT0);
            glDisable(GL_LIGHT1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        else if(drawMode == 1)  // wireframe mode
        {   GLfloat amb[] = { 0,0,0, 0.0 };
            GLfloat diff[]= { 0,0,0};
            GLfloat spec[]= { 0,0,0};
            GLfloat lightpos[]= { 0,0,0 };

            glLightfv(GL_LIGHT1, GL_AMBIENT, amb);
            glLightfv(GL_LIGHT1, GL_SPECULAR, spec);
            glLightfv(GL_LIGHT1, GL_POSITION, lightpos);
            glEnable(GL_LIGHT1);
            glDisable(GL_LIGHT0);

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);
        }
        else                    // point mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);
        }
        break;

    default:
        ;
    }
}