#ifndef RAYTRACING_H_1
#define RAYTRACING_H_1
#include <vector>
#include <cfloat>
#include "mesh.h"

//Welcome to your MAIN PROJECT...
//THIS IS THE MOST relevant code for you!
//this is an important file, raytracing.cpp is what you need to fill out
//In principle, you can do the entire project ONLY by working in these two files

extern Mesh MyMesh; //Main mesh
extern std::vector<Vec3Df> MyLightPositions;
extern Vec3Df MyCameraPosition; //currCamera
extern unsigned int WindowSize_X;//window resolution width
extern unsigned int WindowSize_Y;//window resolution height
extern unsigned int RayTracingResolutionX;  // largeur fenetre
extern unsigned int RayTracingResolutionY;  // largeur fenetre

extern unsigned int recursionLevel;

//use this function for any preprocessing of the mesh.
void init();

// Intersection structure
struct Intersection {
	float distance = FLT_MAX;
	int index = -1;
	Vec3Df intersect;
	Vec3Df normal;
	Vec3Df dest;
	Vec3Df origin;
};

//you can use this function to transform a click to an origin and destination
//the last two values will be changed. There is no need to define this function.
//it is defined elsewhere
void produceRay(int x_I, int y_I, Vec3Df & origin, Vec3Df & dest);

//Intersection with sphere
bool intersectionWithSphere(const Vec3Df &rayOrigin, const Vec3Df &rayDir, const Vec3Df &sphereOrigin, const float &radius);

//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int recursionLevel);

// calculates intersection between a ray from the origin to the dest and the mesh MyMesh
Intersection calcIntersection(Vec3Df origin, Vec3Df dest);

// Shading function(for both hard and soft)
Vec3Df shade(Intersection intersection, int level);

//your heatmap function
float performRayTracingHeatmap(const Vec3Df & origin, const Vec3Df & dest);

//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

// Get the material based on triangle index
Material getMat(int index);

// Return the ReflectionLevel
int getReflectionLevel();

//want keyboard interaction? Here it is...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);

#endif