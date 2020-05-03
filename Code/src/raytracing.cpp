#include <stdio.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#endif
#ifdef WIN32
#include <windows.h>
#include <algorithm>
#endif
#include "raytracing.h"
#include <time.h>
#include "Box.h"


//temporary variables
//these are only used to illustrate 
//a simple debug drawing. A ray 
Vec3Df testRayOrigin;
Vec3Df testRayDestination;
Vec3Df testRayColor;
std::pair<Vec3Df, Vec3Df> bounds;
Box box;

int reflectionLevel = 0;

std::vector<Triangle> triangles;
std::vector<Vertex> vertices;

const unsigned int maxReflectionLevel = recursionLevel;

const Intersection empty = { // empty intersection, since we want to stop the recursion.
	FLT_MAX,
	-1,
	Vec3Df(),
	Vec3Df(),
	Vec3Df(),
	Vec3Df()
};

// For the diffuse, the light we are gonna use 
bool lightMode;

int getReflectionLevel() {
	return reflectionLevel;
}

/* determines whether or not a given ray (origin, dir) intersects with the general
** bounding box.
** Source: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
*/
bool inBoundingBox(Vec3Df origin, Vec3Df dir) {
	Vec3Df min = bounds.first;
	Vec3Df max = bounds.second;

	float xmin = min[0];
	float ymin = min[1];
	float zmin = min[2];

	float xmax = max[0];
	float ymax = max[1];
	float zmax = max[2];

	float txmin = (xmin - origin[0]) / dir[0];
	float txmax = (xmax - origin[0]) / dir[0];

	if (txmin > txmax) {
		std::swap(txmin, txmax);
	}

	float tymin = (ymin - origin[1]) / dir[1];
	float tymax = (ymax - origin[1]) / dir[1];

	if (tymin > tymax) {
		std::swap(tymin, tymax);
	}

	if ((txmin > tymax) || (tymin > txmax)) {
		return false;
	}

	if (tymin > txmin) {
		txmin = tymin;
	}

	if (tymax < txmax) {
		txmax = tymax;
	}

	float tzmin = (zmin - origin[2]) / dir[2];
	float tzmax = (zmax - origin[2]) / dir[2];

	if (tzmin > tzmax) {
		std::swap(tzmin, tzmax);
	}

	if ((txmin > tzmax) || (tzmin > txmax)) {
		return false;
	}

	return true;
}

std::pair<Vec3Df, Vec3Df> boundingBox() {
	std::vector<Vertex> vertices = MyMesh.vertices;
	Vec3Df min = vertices[0].p;
	Vec3Df max = vertices[0].p;

	for (int i = 1; i < vertices.size(); i++) {
		Vec3Df current = vertices[i].p;
		for (int j = 0; j < 3; j++) {
			if (current[j] < min[j]) {
				min[j] = current[j];
				
			}
		}
		
		for (int j = 0; j < 3; j++) {
			if (current[j] > max[j]) {
				max[j] = current[j];
				
			}
		}
	}
	std::cout << "Bounding box min value: " << min << std::endl;
	std::cout << "Bounding box max value: " << max << std::endl;
	return std::pair<Vec3Df, Vec3Df>(min, max);
}

float intPow(float base, int exp) {
	float ret = 1;
	for (int i = 0; i < exp; i++) {
		ret *= base;
	}
	return ret;
}

/*
** http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
** https://graphics.stanford.edu/~mdfisher/Code/Engine/Plane.cpp.html
** https://stackoverflow.com/questions/7168484/3d-line-segment-and-plane-intersection
** Intresection with a plane
*/
Vec3Df intersectionWithPlane(Vec3Df &dir, Vec3Df &origin, const Vec3Df &planeNormal, Vec3Df &planePoint)
{
	dir.normalize();
	origin.normalize();
	planePoint.normalize();

	float t = Vec3Df::dotProduct((planePoint - origin), planeNormal) / Vec3Df::dotProduct(dir, planeNormal);
	Vec3Df res = origin + t * dir;
	return res;
}


/*
** http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
** https://graphics.stanford.edu/~mdfisher/Code/Engine/Plane.cpp.html
** https://stackoverflow.com/questions/7168484/3d-line-segment-and-plane-intersection
** http://www.cplusplus.com/forum/general/96882/
** Check wether we have intersection with a sphere or not.
*/
bool intersectionWithSphere(const Vec3Df &rayOrigin, const Vec3Df &rayDir, const Vec3Df &sphereOrigin, const float &radius)
{
	float a = rayOrigin[0] - sphereOrigin[0];
	float b = rayOrigin[1] - sphereOrigin[1];
	float c = rayOrigin[2] - sphereOrigin[2];

	float A = (rayDir[0] * rayDir[0]) + (rayDir[1] * rayDir[1]) + (rayDir[2] * rayDir[2]);
	float B = 2 * ((rayDir[0] * a) + (rayDir[1] * b) + (rayDir[2] * c));
	float C = (a*a) + (b*b) + (c*c) - (radius*radius);

	float discriminant = (B*B) - (4 * A*C);
	if (discriminant >= 0) {
		return true;
	}
	else {
		return false;
	}

	/*
	//If the ray has one intersection point the two intersection points are the same
	//If the ray has two intersection points the two intersection points are different
	//If the ray has no intersection points the intersection points are null;
	*/
}

// Get the material based on triangle index
Material getMat(int index) {
	// index = (triangles.size() - index + 1) % (triangles.size() - 1);
	int matIndex = MyMesh.triangleMaterials[index];
	return MyMesh.materials[matIndex];
}

/*
** Calculate the diffuse and return it with it's color
** Source: CG slides
*/
Vec3Df diffuse(Intersection intersection, Vec3Df lightPos) {
	Intersection lightIntersection = calcIntersection(intersection.intersect, lightPos);
	if (lightIntersection.index != -1) {
		return Vec3Df(0, 0, 0);
	}

	Vec3Df lightRay = lightPos - intersection.intersect;
	lightRay.normalize();

	float diffuseTerm = Vec3Df::dotProduct(intersection.normal, lightRay);
	if (diffuseTerm < 0) {
		return Vec3Df(0, 0, 0);
	}
	return getMat(intersection.index).Kd() * diffuseTerm;
}

Vec3Df blinnPhong(Intersection intersection, const Vec3Df & lightPos, const Vec3Df & cameraPos) {
	Intersection lightIntersection = calcIntersection(intersection.intersect, lightPos);
	if (lightIntersection.index != -1) {
		return Vec3Df(0, 0, 0);
	}
	Material mat = getMat(intersection.index);
	
	Vec3Df camRay = cameraPos - intersection.intersect;
	Vec3Df lightRay = lightPos - intersection.intersect;
	camRay.normalize();
	lightRay.normalize();
	Vec3Df halfVector = lightRay + camRay; // no need to devide by 2 due to normalization later on
	
	float specularTerm = Vec3Df::dotProduct(intersection.normal, halfVector.unit());
	
	if (Vec3Df::dotProduct(intersection.normal, lightRay.unit()) < 0) {
		return Vec3Df(0, 0, 0);
	}
	float powSpecular = pow(specularTerm, mat.Ns());
	return mat.Ks() * powSpecular;
}

/*
** This function is optional but it helps a lot
** Function to add offsets towards the direction of the reflection so the ray doesnt hit itself
*/
void Offset(Vec3Df* hit, Vec3Df* dest) {
	Vec3Df vector = (*dest) - (*hit);
	vector.normalize();
	vector *= 0.01f;
	*hit += vector;
}

/*
** https://cs.nyu.edu/~perlin/courses/fall2005ugrad/phong.html
** http://ray-tracing-conept.blogspot.com/2015/01/hard-and-soft-shadows.html
** Implement the hard shadows calculater
*/
Vec3Df shade(Intersection intersection, int level) {
	Vec3Df result = Vec3Df(0, 0, 0);

	Material material = getMat(intersection.index);

	for (unsigned int i = 0; i < MyLightPositions.size(); i++) {
		Vec3Df L = MyLightPositions[i];

		/* Start of shading block */
		if (material.has_Ka()) {
			result += material.Ka();
		}

		if (material.has_Kd()) {
			result += diffuse(intersection, L);
		}

		if (material.has_Ks()) {
			result += blinnPhong(intersection, L, MyCameraPosition);
		}
		/* End of shading block */
	}

	if (result.p[0] > 1) {
		result.p[0] = 1;
	}
	if (result.p[1] > 1) {
		result.p[1] = 1;
	}
	if (result.p[2] > 1) {
		result.p[2] = 1;
	}
	if (result.p[0] < 0) {
		result.p[0] = 0;
	}
	if (result.p[1] < 0) {
		result.p[1] = 0;
	}
	if (result.p[2] < 0) {
		result.p[2] = 0;
	}
	return result;
}

/*
** Need to calculate how much gets reflected and how much gets reflected
** source: https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
*/
float reflectionFraction(Intersection intersection) {
	float muRatio = getMat(intersection.index).Ni();

	float mu1 = 1;			// still need to figure out how much this is from the material
	float mu2 = 0.75;

	Vec3Df incident = intersection.origin - intersection.intersect;
	Vec3Df normal = intersection.normal;

	bool TIR = false;		// ???

	if (TIR) {
		return 1;
	}

	float cosThetaI = -Vec3Df::dotProduct(incident, normal);
	float sinThetaISquared = intPow(muRatio, 2) * (1 - intPow(muRatio, 2));

	float R0 = intPow((mu1 - mu2) / (mu1 + mu2), 2);

	if (mu1 <= mu2) {
		return R0 + (1 - R0) * intPow(1 - cosThetaI, 5);
	}
	else {
		float cosThetaT = std::sqrt(1 - sinThetaISquared);
		return R0 + (1 - R0) * intPow(1 - cosThetaT, 5);
	}
}

float calcMuFraction(Intersection intersection) {
	float ret = getMat(intersection.index).Ni();
	if (Vec3Df::dotProduct(intersection.normal, intersection.dest - intersection.origin) < 0) {
		return ret;
	}
	else {
		return 1 / ret;
	}
}

/*
** calculates reflection at intersection
** source: https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
*/
Vec3Df calcReflection(Intersection intersection, int recursionLevel) {
	Vec3Df dir = intersection.intersect - intersection.origin;
	dir.normalize();

	//calculate the reflection vector
	Vec3Df refDir = dir - (2 * Vec3Df::dotProduct(intersection.normal, dir)*intersection.normal);
	//hit point on mirror(Or another reflective surface)
	Vec3Df hit = intersection.intersect;
	//destination point
	Vec3Df dest = hit + refDir;

	return performRayTracing(hit, dest, recursionLevel);
}

/*
** calculates refraction at intersection
** source: https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
*/
Vec3Df calcRefraction(Intersection intersection, int recursionLevel) {
	Vec3Df incident = intersection.intersect - intersection.origin;
	Vec3Df normal = intersection.normal;

	float muRatio = calcMuFraction(intersection);

	float cosTheta = -Vec3Df::dotProduct(incident, normal);
	float sinThetaSquared = intPow(muRatio, 2) * (1 - intPow(cosTheta, 2));

	if (sinThetaSquared > 1) {
		return Vec3Df(0, 0, 0);
	}

	Vec3Df refractDir = muRatio * incident + (muRatio * cosTheta - std::sqrt(1 - sinThetaSquared)) * normal;

	return performRayTracing(intersection.intersect, intersection.intersect + refractDir, recursionLevel);
}

//function to get max vec
Vec3Df computeMaxVec() {
	std::vector<Triangle> triangles = MyMesh.triangles;
	std::vector<Vertex> vertices = MyMesh.vertices;
	Vec3Df max;

	max[0] = vertices.at(triangles.at(0).v[0]).p[0];
	max[1] = vertices.at(triangles.at(0).v[0]).p[1];
	max[2] = vertices.at(triangles.at(0).v[0]).p[2];

	for (int i = 0; i < triangles.size(); i++) {
		Triangle curTriangle = triangles.at(i);

		for (int vernumber = 0; vernumber < 3; vernumber++) {
			if (max[0]<vertices.at(curTriangle.v[vernumber]).p[0]) { max[0] = vertices.at(curTriangle.v[vernumber]).p[0]; }
			if (max[1]<vertices.at(curTriangle.v[vernumber]).p[1]) { max[1] = vertices.at(curTriangle.v[vernumber]).p[1]; }
			if (max[2]<vertices.at(curTriangle.v[vernumber]).p[2]) { max[2] = vertices.at(curTriangle.v[vernumber]).p[2]; }
		}
	}
	return max;
}

//function to get min vec
Vec3Df computeMinVec() {
	std::vector<Triangle> triangles = MyMesh.triangles;
	std::vector<Vertex> vertices = MyMesh.vertices;
	Vec3Df min;

	min[0] = vertices.at(triangles.at(0).v[0]).p[0];
	min[1] = vertices.at(triangles.at(0).v[0]).p[1];
	min[2] = vertices.at(triangles.at(0).v[0]).p[2];

	for (int i = 0; i < triangles.size(); i++) {
		Triangle curTriangle = triangles.at(i);

		for (int vernumber = 0; vernumber < 3; vernumber++) {
			if (min[0]>vertices.at(curTriangle.v[vernumber]).p[0]) { min[0] = vertices.at(curTriangle.v[vernumber]).p[0]; }
			if (min[1]>vertices.at(curTriangle.v[vernumber]).p[1]) { min[1] = vertices.at(curTriangle.v[vernumber]).p[1]; }
			if (min[2]>vertices.at(curTriangle.v[vernumber]).p[2]) { min[2] = vertices.at(curTriangle.v[vernumber]).p[2]; }
		}
	}
	return min;
}

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should, if they 
	//are exported as WavefrontOBJ.
	//PLEASE ADAPT THE LINE BELOW TO THE FULL PATH OF THE dodgeColorTest.obj
	//model, e.g., "C:/temp/myData/GraphicsIsFun/dodgeColorTest.obj", 
	//otherwise the application will not load properly
	//MyMesh.loadMesh("dodgeColorTest.obj", true);
    //MyMesh.loadMesh("scene3.obj", true);
	//MyMesh.loadMesh("cube.obj", true);
	//MyMesh.loadMesh("refraction.obj", true);
	MyMesh.loadMesh("scene.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(Vec3Df(-1.75357, 2.82853, 3.55423));
	MyLightPositions.push_back(Vec3Df(4.04358, 1.82305, 0.34666));
	MyLightPositions.push_back(Vec3Df(-1.45883, 3.66948, -0.0327641));


	bounds = boundingBox();

	triangles = MyMesh.triangles;
	vertices = MyMesh.vertices;

	std::vector<int> triangleIndices(triangles.size());
	for (int i = 0; i < triangles.size(); i++) {
		triangleIndices[i] = i;
	}


	box = Box(triangleIndices, computeMinVec(), computeMaxVec());
}

float performRayTracingHeatmap(const Vec3Df & origin, const Vec3Df & dest) {
	float start = clock();
	performRayTracing(origin, dest, maxReflectionLevel);
	float end = clock();
	return start - end;
}

Intersection calcIntersection(Vec3Df origin, Vec3Df dest) {
	Vec3Df dir = dest - origin;
	Vec3Df intersectionPoint;
	Vec3Df intersectionNormal;

    /////////////
    std::vector<int> TrianglesFromTheKDtree;
    if(box.hit(box, origin, dir)){
        TrianglesFromTheKDtree = box.trace(origin, dir);
    }
    /////////////
    
    
	if (!inBoundingBox(origin, dir)) {
		return empty;
	}

	float smallestT = INFINITY;
	int relevantIndex = -1;

    for (int i = 0; i < TrianglesFromTheKDtree.size(); i++) {
        Triangle curTriangle = triangles[TrianglesFromTheKDtree.at(i)];

        Vec3Df v0 = vertices.at(curTriangle.v[0]).p;
        Vec3Df v1 = vertices.at(curTriangle.v[1]).p;
        Vec3Df v2 = vertices.at(curTriangle.v[2]).p;
		
		Vec3Df aEdge = v1 - v0;
        Vec3Df bEdge = v2 - v0;

        Vec3Df normal = Vec3Df::crossProduct(aEdge, bEdge).unit();
        float D = Vec3Df::dotProduct(normal, v0);

		float t = (D - Vec3Df::dotProduct(origin, normal)) / Vec3Df::dotProduct(dir, normal);

		if (Vec3Df::dotProduct(dir, normal) < 0 && t > 0 && t < smallestT) {
			Vec3Df p = origin + t * dir;

			Vec3Df pv0 = p - v0;

			// using Cramers rule, see: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates

			float d00 = Vec3Df::dotProduct(aEdge, aEdge);
			float d01 = Vec3Df::dotProduct(aEdge, bEdge);
			float d11 = Vec3Df::dotProduct(bEdge, bEdge);
			float d20 = Vec3Df::dotProduct(pv0, aEdge);
			float d21 = Vec3Df::dotProduct(pv0, bEdge);

			float denom = d00 * d11 - d01 * d01;

			float a = (d11 * d20 - d01 * d21) / denom;
			float b = (d00 * d21 - d01 * d20) / denom;

			if (a >= 0 && b >= 0 && a + b <= 1) {
				smallestT = t;
				relevantIndex = TrianglesFromTheKDtree[i];

				intersectionPoint = p;
				intersectionNormal = normal;
			}

		}
	}

	if (smallestT == INFINITY) {
		return empty;
	}
	else {
		return {
			smallestT,
			relevantIndex,
			intersectionPoint,
			intersectionNormal,
			dest,
			origin
		};
	}
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int recursionLevel)
{
	reflectionLevel = 0;
	recursionLevel--;
	
	if (recursionLevel < 1) {
		return Vec3Df(0, 0, 0);
	}

	Intersection intersection = calcIntersection(origin, dest);

	if (intersection.index != -1) {
		Vec3Df ret = Vec3Df(0, 0, 0);
		Material mat = getMat(intersection.index);

		// decide what to calculate based on the illum modes
		// source: http://paulbourke.net/dataformats/mtl/
		switch (mat.illum()) {
			case 0:
				ret = mat.Kd();
				break;
		
			case 1:
				for (int i = 0; i < MyLightPositions.size(); i++) {
					ret += diffuse(intersection, MyLightPositions[i]);
				}
				break;

			case 2:
				ret = shade(intersection, recursionLevel);
				break;

			case 3:
			case 4:
			case 5:
				ret = shade(intersection, recursionLevel) + 
					mat.Ks() * calcReflection(intersection, recursionLevel);
				break;

			case 6:
				ret = shade(intersection, recursionLevel) +
					mat.Ks() * calcReflection(intersection, recursionLevel) +
					(Vec3Df(1, 1, 1) - mat.Ks()) * mat.Tr() * calcRefraction(intersection, recursionLevel);
				break;
		}
		
		return ret;
	}
	else {
		return Vec3Df(0, 0, 0);
	}
}

void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame

	//let's draw the mesh
	MyMesh.draw();
	
	//let's draw the lights in the scene as points
	glPushAttrib(GL_ALL_ATTRIB_BITS); //store all GL attributes
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i=0;i<MyLightPositions.size();++i)
		glVertex3fv(MyLightPositions[i].pointer());
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state. 
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous 
	//state after the pop.


	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(testRayColor[0], testRayColor[1], testRayColor[2]);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glColor3f(testRayColor[0], testRayColor[1], testRayColor[2]);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();
	
	//draw whatever else you want...
	////glutSolidSphere(1,10,10);
	////allows you to draw a sphere at the origin.
	////using a glTranslate, it can be shifted to whereever you want
	////if you produce a sphere renderer, this 
	////triangulated sphere is nice for the preview
}

//yourKeyboardFunc is used to deal with keyboard input.
//t is the character that was pressed
//x,y is the mouse position in pixels
//rayOrigin, rayDestination is the ray that is going in the view direction UNDERNEATH your mouse position.
//
//A few keys are already reserved: 
//'L' adds a light positioned at the camera location to the MyLightPositions vector
//'l' modifies the last added light to the current 
//    camera position (by default, there is only one light, so move it with l)
//    ATTENTION These lights do NOT affect the real-time rendering. 
//    You should use them for the raytracing.
//'r' calls the function performRaytracing on EVERY pixel, using the correct associated ray. 
//    It then stores the result in an image "result.bmp".
//    Initially, this function is fast (performRaytracing simply returns 
//    the target of the ray - see the code above), but once you replaced 
//    this function and raytracing is in place, it might take a 
//    while to complete...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination)
{

	Intersection intersection = calcIntersection(rayOrigin, rayDestination);

	testRayOrigin=MyCameraPosition;	
	testRayDestination=intersection.intersect;
	testRayColor = performRayTracing(testRayOrigin, testRayDestination, maxReflectionLevel);
	
	
	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;	
}
