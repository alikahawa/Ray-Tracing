

//PLEASE IGNORE THIS FILE !!!! 
//IT DEALS WITH THE TRACKBALL
//IT IS COMPLEX and has French comments ;) 
//good times...



/** \file mouse.h
 
Tools for using a trackball/mouse in an OpenGL window
 
*/
#ifndef MOUSE_H
#define MOUSE_H
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#endif
#ifdef WIN32
#include <windows.h>
#endif
#include <math.h>
#include "matrix.h"
#include "stdio.h"
#include "Vec3D.h"
static const float speedfact = 0.2;

/** Display function */
void display();

/** View transformation */
GLdouble tb_matrix[16] =   { 1,0,0,0,
                             0,1,0,0,
                             0,0,1,0,
                             0,0,0,1  };
GLdouble tb_inverse[16] =  { 1,0,0,0,
                             0,1,0,0,
                             0,0,1,0,
                             0,0,0,1  };

/** Mouse variables */
int tb_oldX, tb_oldY, tb_rotateXY=0, tb_translateXY=0, tb_moveZ=0;


/** Initialize the current transformation matrix with the initial viewpoint */
void tbInitTransform()
{
    glGetDoublev( GL_MODELVIEW_MATRIX, tb_matrix );
    inverse( tb_matrix, tb_inverse );
}

/** Apply view transformation */
void tbVisuTransform()
{
    glMultMatrixd( tb_matrix );
};

/** Print help */
void tbHelp()
{
    printf("Left button: turn in XY,\n");
    printf("Right button: translate in XY,\n");
    printf("Middle button: move along Z.\n");
}

/** Process mouse buttons */
void tbMouseFunc( int button, int state, int x, int y )
{
    /* Left button press */
    if( button==GLUT_LEFT_BUTTON && state==GLUT_DOWN )
    {
        tb_rotateXY = 1;
        tb_oldX = x;
        tb_oldY = y;
    }
    /* Left button release */
    else if( button==GLUT_LEFT_BUTTON && state==GLUT_UP )
    {
        tb_rotateXY = 0;
    }
    /* Middle button press */
    if( button==GLUT_MIDDLE_BUTTON && state==GLUT_DOWN )
    {
        tb_moveZ = 1;
        tb_oldX = x;
        tb_oldY = y;
    }
    /* Middle button release */
    else if( button==GLUT_MIDDLE_BUTTON && state==GLUT_UP )
    {
        tb_moveZ = 0;
    }
    /* Right button press */
    else if( button==GLUT_RIGHT_BUTTON && state==GLUT_DOWN )
    {
        tb_translateXY = 1;
        tb_oldX = x;
        tb_oldY = y;
    }
    /* Right button release */
    else if( button==GLUT_RIGHT_BUTTON && state==GLUT_UP )
    {
        tb_translateXY = 0;
    }
}


/** Process mouse motion for translation */
void tbMotionFunc( int x, int y )
{
    double dx,dy,nrm, tx,ty,tz;

    if( tb_rotateXY || tb_translateXY || tb_moveZ )
    {
        /* motion amount */
        dx = x - tb_oldX;
        dy = tb_oldY - y; /* Vertical axis direction is inverted */

        if( tb_rotateXY )
        {
            tx = tb_matrix[12];
            tb_matrix[12]=0;
            ty = tb_matrix[13];
            tb_matrix[13]=0;
            tz = tb_matrix[14];
            tb_matrix[14]=0;

            nrm = sqrt( dx*dx+dy*dy+dx*dx+dy*dy )*speedfact;
            glLoadIdentity();
            glRotatef( nrm, -dy, dx, 0 );/* Axis perpendicular to motion */
            glMultMatrixd( tb_matrix );
            glGetDoublev( GL_MODELVIEW_MATRIX, tb_matrix );

            tb_matrix[12] = tx;
            tb_matrix[13] = ty;
            tb_matrix[14] = tz;
        }
        else if( tb_translateXY )
        {
            tb_matrix[12] += dx/100.0*speedfact;
            tb_matrix[13] += dy/100.0*speedfact;
        }
        else if( fabs(dx)>fabs(dy) )
        { // rotation z
            tx = tb_matrix[12];
            tb_matrix[12]=0;
            ty = tb_matrix[13];
            tb_matrix[13]=0;
            tz = tb_matrix[14];
            tb_matrix[14]=0;

            glLoadIdentity();
            glRotatef( dx, 0,0,-1 );/* Axis perpendicular to the screen */
            glMultMatrixd( tb_matrix );
            glGetDoublev( GL_MODELVIEW_MATRIX, tb_matrix );

            tb_matrix[12] = tx;
            tb_matrix[13] = ty;
            tb_matrix[14] = tz;
        }
        else if( fabs(dy)>fabs(dx) )
        {
            tb_matrix[14] -= dy/100.0*speedfact;
        }
        tb_oldX = x;
        tb_oldY = y;
        inverse( tb_matrix, tb_inverse );
        glutPostRedisplay();
    }
}

/** Process mouse motion for rotation */
void tbRotate( double angle, double x, double y, double z )
{
    double tx,ty,tz;

    tx = tb_matrix[12];
    tb_matrix[12]=0;
    ty = tb_matrix[13];
    tb_matrix[13]=0;
    tz = tb_matrix[14];
    tb_matrix[14]=0;

    glLoadIdentity();
    glRotatef( angle, x, y, z );
    glMultMatrixd( tb_matrix );
    glGetDoublev( GL_MODELVIEW_MATRIX, tb_matrix );

    tb_matrix[12] = tx;
    tb_matrix[13] = ty;
    tb_matrix[14] = tz;

    inverse( tb_matrix, tb_inverse );
    glutPostRedisplay();
}

/// Projection to world space
void tbProject( const GLdouble *m, const GLdouble* p, GLdouble* q )
{
    double pp[4];
    //cout<<"tb, matrix: "; printMatrix(tb_matrix); cout<<endl;
    //cout<<"tb, inverse: "; printMatrix(tb_inverse); cout<<endl;
    project( m, p, pp );
    //cout<<"proj: "<<pp[0]<<", "<<pp[1]<<", "<<pp[2]<<", "<<pp[3]<<endl;
    project( tb_inverse, pp, q );
    //cout<<"projRep: "<<q[0]<<", "<<q[1]<<", "<<q[2]<<", "<<q[3]<<endl;
}

void tbProject( const GLdouble* p, GLdouble* q )
{
    //cout<<"proj: "<<pp[0]<<", "<<pp[1]<<", "<<pp[2]<<", "<<pp[3]<<endl;
    project( tb_inverse, p, q );
    //cout<<"projRep: "<<q[0]<<", "<<q[1]<<", "<<q[2]<<", "<<q[3]<<endl;
}
Vec3Df getCameraPosition()
{
	const GLdouble p[]={0,0,0,1};
	GLdouble LightP[4];
	tbProject(p, LightP);
	Vec3Df LightPos;
	LightPos[0]=LightP[0];
	LightPos[1]=LightP[1];
	LightPos[2]=LightP[2];
	return LightPos;
}
Vec3Df getWorldPositionOfPixel(unsigned int px, unsigned int py)
{

	double mv[16];
	double pr[16];
	int vp[4];
	glGetDoublev(GL_MODELVIEW_MATRIX,mv);
	glGetDoublev(GL_PROJECTION_MATRIX,pr);
	glGetIntegerv(GL_VIEWPORT,vp);

	double x,y,z;
	gluUnProject(double(px),double(py),0,mv,pr,vp,&x,&y,&z);



	return Vec3Df(x,y,z);
}
#endif
