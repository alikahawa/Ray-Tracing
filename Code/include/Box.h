//
//  Box.h
//  Assignment_5_Project
//
//  Created by Mohamad Awab Alkhiami & Ghiyath Alaswad on 21/06/2018.
//

#ifndef Box_h
#define Box_h
#include<algorithm>
class Box {
private:
    std::vector<int> triangleIndices;
    std::vector<Box> Boxes;
    Vec3Df min;
    Vec3Df max;
    
public:
    Box() {
        
    }
    
    Box(std::vector<int> _triangles, Vec3Df _min, Vec3Df _max) {
        this->min = _min;
        this->max = _max;
        this->triangleIndices = _triangles;
        
        if (triangleIndices.size() < 50) {
            //it is fine because the triangles vector has less than 50 triangles
        }
        else {
            // Longest axis
            char _theLongestAxes = findLongestAxis();
            
            Vec3Df leftMAx;
            Vec3Df rightMin;
            
            computeTheLeftMaxAndRightMin(min, max, leftMAx, rightMin, _theLongestAxes);
            
            
            std::vector<int> _LTriangles;
            std::vector<int> _RTriangles;
            
            split(&this->triangleIndices, &_LTriangles, &_RTriangles, this->min, leftMAx, rightMin, this->max);
            
            
            Box _L = Box(_LTriangles, this->min, leftMAx);
            Box _R = Box(_RTriangles, rightMin, this->max);
            Boxes.push_back(_L);
            Boxes.push_back(_R);
            
        }
    }
    
    
    
    void split(std::vector<int> *triangles, std::vector<int> *_LTriangles, std::vector<int> *_RTriangles, Vec3Df Lmin, Vec3Df Lmax, Vec3Df Rmin, Vec3Df Rmax) {
        
        std::vector<int> _Temp;
        
        for (const int& _tri : *triangles) {
            if (hasVertexInBox(_tri, Lmin, Lmax) && hasVertexInBox(_tri, Rmin, Rmax)) {//in both boxes
                _Temp.push_back(_tri);
            }
            else if (hasVertexInBox(_tri, Lmin, Lmax)) {//in left box
                _LTriangles->push_back(_tri);
            }
            else if (hasVertexInBox(_tri, Rmin, Rmax)) {//in right box
                _RTriangles->push_back(_tri);
            }
        }
        *triangles = _Temp;
    }
    
    
    //function that takes thie min max vec3Df and decide that we gonna split our box according to X,Y or Z
    char findLongestAxis() {
        Vec3Df min = this->min;
        Vec3Df max = this->max;
        
        float _X = max[0] - min[0];
        float _Y = max[1] - min[1];
        float _Z = max[2] - min[2];
        if (_X == std::max(_X, _Y) && _X == std::max(_X, _Z)) {
            return 'x';
        }
        else if (_Y == std::max(_X, _Y) && _Y == std::max(_Y, _Z)) {
            return 'y';
        }
        else if (_Z == std::max(_X, _Z) && _X == std::max(_Y, _Z)) {
            return 'z';
        }
        
        return 'x';
    }
    
    
    //function to tell wether the ray is hitting the nox or not
    bool hit(Box box, Vec3Df origin, Vec3Df dir) {
        Vec3Df min = box.min;
        Vec3Df max = box.max;
        
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
    
    void computeTheLeftMaxAndRightMin(Vec3Df &finalmin, Vec3Df &finalmax, Vec3Df &leftMAx, Vec3Df &rightMin, char _theLongestAxes) {
        
        switch (_theLongestAxes) {
            case 'x':
                leftMAx[0] = finalmin[0] + ((finalmax[0] - finalmin[0]) / 2);
                leftMAx[1] = max[1];
                leftMAx[2] = max[2];
                
                rightMin[0] = min[0] + ((max[0] - min[0]) / 2);
                rightMin[1] = min[1];
                rightMin[2] = min[2];
                break;
                
            case 'y':
                leftMAx[0] = max[0];
                leftMAx[1] = min[1] + ((max[1] - min[1]) / 2);
                leftMAx[2] = max[2];
                
                rightMin[0] = min[0];
                rightMin[1] = min[1] + ((max[1] - min[1]) / 2);
                rightMin[2] = min[2];
                break;
                
            case 'z':
                leftMAx[0] = max[0];
                leftMAx[1] = max[1];
                leftMAx[2] = min[2] + ((max[2] - min[2]) / 2);
                
                rightMin[0] = min[0];
                rightMin[1] = min[1];
                rightMin[2] = min[2] + ((max[2] - min[2]) / 2);
                break;
        }
    }
    
    //function to determine wether a vertix or more of a triangle is in a box
    bool hasVertexInBox(int triangle, Vec3Df _min, Vec3Df _max) {
        for (int i = 0; i < 3; i++) {
            if (isVertexInBox(MyMesh.vertices[MyMesh.triangles[triangle].v[i]], _min, _max)) {
                return true;
            }
        }
        return false;
    }
    
    bool isVertexInBox(Vertex vertex, Vec3Df _min, Vec3Df _max) {
        Vec3Df pos = vertex.p;
        return pos[0] >= _min[0] && pos[0] <= _max[0] &&
        pos[1] >= _min[1] && pos[1] <= _max[1] &&
        pos[2] >= _min[2] && pos[2] <= _max[2];
    }
    
    
    std::vector<int> trace(const Vec3Df & origin, const Vec3Df & dir) {
        std::vector<int> returnable;
        
        //        if (hit(this->Boxes[0], origin,dir)  && hit(this->Boxes[1], origin,dir)) {
        //            //box has items or box hit ???
        //            if(hit(this->Boxes[0], origin,dir)) returnable.push_back(this->Boxes[0].trace(origin,dir));
        //            if(hit(this->Boxes[1], origin,dir)) returnable.push_back(this->Boxes[1].trace(origin,dir));
        //        }
        
        
        
            if(this->triangleIndices.size()>0){
                returnable.insert(returnable.end(), this->triangleIndices.begin(),this->triangleIndices.end());
            }
        
        if(this->Boxes.size()>0){
            if(hit(this->Boxes[0], origin,dir)){
                std::vector<int> _TEMPVECTOR = this->Boxes[0].trace(origin,dir);
                returnable.insert(returnable.end(),_TEMPVECTOR.begin(),_TEMPVECTOR.end());
            }
            
            if (hit(this->Boxes[1], origin,dir)){
                std::vector<int> _TEMPVECTOR = this->Boxes[1].trace(origin,dir);
                returnable.insert(returnable.end(),_TEMPVECTOR.begin(),_TEMPVECTOR.end());
            }
        }
        
        return returnable;
    }
    
};

#endif /* Box_h */
