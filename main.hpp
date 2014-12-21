#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <stdint.h>
#include <FlexLexer.h>
#include <string>
#include <list>
#include <stack>
#include "Eigen/Dense"
#include "data.h"
#include <vector>
#include <stdlib.h>
#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"
#include <unordered_map>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;



void initLights();

void set_vertex (GLfloat *, int, GLfloat *, GLfloat *,int);

void makeCheckImage(void);

void loadTextureFromFile(char *filename);

void drawScene(void);

void drawScenes();

void drawScene2();

void init(int argc, char **argv);

void start(int argc, char **argv);

void mouseCB(int button, int state, int x, int y);

void mouseMotionCB(int x, int y);

void timerCB(int millisec);

void init_material(void);

void keyboardCB(unsigned char key, int x, int y);

void loadTextureFromFile(char *filename);

//   ---------helper.cpp----------
void vectorXdToArray (VectorXd v, GLfloat* g);

void vectorToArray (vector<float> v, GLfloat*g);

VectorXd ArrayToVectorXd (GLfloat *a, int n);

VectorXd vectorToVectorXd (vector<float> vec);

// ----------structures----------
struct object {
  vector<vector<float> > nums;
  vector<vector<int> > facets;
};

struct triangle {
  int tri_num;
  vector<int> vertex_indices;
  vector<vector<float> > vertices_coords;
  VectorXd normal;
  
  void calculate_normal() {

  	Vector3d v1 = vectorToVectorXd(vertices_coords[0]);
    Vector3d v2 = vectorToVectorXd(vertices_coords[1]);
  	Vector3d v3 = vectorToVectorXd(vertices_coords[2]);

  	Vector3d edge1 = v2 - v1;
  	Vector3d edge2 = v3 - v1;

  	normal = (edge1.cross(edge2)).normalized();

  }
};


#endif