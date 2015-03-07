/* CS6555 Computer Animation, Fall 2014
 * Lab 1£ºBasic Key - Framing Motion Control System
 * Due by Sept 25, 2014
 * Edited by Fei Yan (Email:hcifaye@gwu.edu)
 * reference: Lab 0's code
*/

// window
#include "stdafx.h"

// standard
#include <assert.h>
#include <math.h>

// glut
#include <GL/glut.h>

//================================
// global variables
//================================
// screen size
int g_screenWidth = 0;
int g_screenHeight = 0;

// Number of Points for Spline
static int points = 0; //points index started from 0
static int number = 7; //The total number of points 

// time variables for generate Q(t)
static GLfloat t = 0;

// The final M matrix for glLoadMatrixf()
static GLfloat M[16] = { 0 };

// The Catmul-Rom Spline M Marix
static GLfloat CRSplineM[16] = { -0.5f, 1.0f, -0.5f, 0.0f,	  // Column 1
								1.5f, -2.5f, 0.0f, 1.0f,      // Column 2
								-1.5f, 2.0f, 0.5f, 0.0f,      // Column 3
								0.5f, -0.5f, 0.0f, 0.0f };    // Column 4

// The B Spline M Marix
static GLfloat BSplineM[16]= { -1.0f/6.0f, 3.0f/6.0f, -3.0f/6.0f, 1.0f/6.0f,  // Column 1
								3.0f/6.0f, -6.0f/6.0f, 0.0f/6.0f, 4.0f/6.0f,  // Column 2
								-3.0f/6.0f, 3.0f/6.0f, 3.0f/6.0f, 1.0f/6.0f,  // Column 3
								1.0f/6.0f, 0.0f/6.0f, 0.0f /6.0f, 0.0f/6.0f };// Column 4

// 7 Pionts in Quternion£º the first 4 numbers represent w, x, y, z in quaternion, and the rest 3 numbers represent the position x,y,z in world Cartisian System
static GLfloat point_quaternion[7][7] = { { 1, 0, 0, 0, -5, 0, -5 },   //point 1
										  { 0, 1, 0, 0, -3, 3, -10 },  //point 2
										  { 0, 0, 1, 0, -1, 1, -15 },   //point 3
										  { 0, 0, 0, 1, 0, -5, -20 },  //point 4
										  { 0, 0, 1, 0, 1, 1, -15 },   //point 5
										  { 0, 1, 0, 0, 3, 3, -10 },  //point 6
										  { 1, 0, 0, 0, 5, 0, -5 } }; //point 7


// 7 Pionts in Euler Angle£º the first 3 numbers represent x_angle, y_angle, z_angle in Euler angle, and the rest 3 numbers represent the position x,y,z in world Cartisian System
static GLfloat point_euler[7][6] = { { 90, 0, 45, -5, 0, -5 },		//point 1
									 { 70, 20, 65, -3, 3, -10 },	//point 2
									 { 50, 40, 85, -1, 1, -15 },	//point 3
									 { 30, 60, 105, 0, -5, -20 },	//point 4
									 { 50, 40, 85, 1, 1, -15 },		//point 5
									 { 70, 20, 65, 3, 3, -10 },		//point 6
									 { 90, 0, 45, 5, 0, -5 }, };	//point 7


//==================================================================================================
// Blending Function : Q(t) = T*M*G, finding out the vector position of time t
//==================================================================================================
GLfloat blend(GLfloat T[4], GLfloat MS[16], GLfloat G[4])
{
	// B[4] is the result of T*M
	GLfloat B[4] = { 0 };
	B[0] = T[0] * MS[0] + T[1] * MS[1] + T[2] * MS[2] + T[3] * MS[3];	 //column 1
	B[1] = T[0] * MS[4] + T[1] * MS[5] + T[2] * MS[6] + T[3] * MS[7];	 //column 2
	B[2] = T[0] * MS[8] + T[1] * MS[9] + T[2] * MS[10] + T[3] * MS[11];  //column 3
	B[3] = T[0] * MS[12] + T[1] * MS[13] + T[2] * MS[14] + T[3] * MS[15];//column 4

	// Generate the result of T*M*G
	GLfloat Qt = B[0] * G[0] + B[1] * G[1] + B[2] * G[2] + B[3] * G[3];

	return Qt;
}

//==================================================================================================
// Unit Quaternion : generate unit quaternion with given quaternion array
//==================================================================================================
void Normalization(GLfloat N_tempM[7]) 
{
	GLfloat squa_quaterion = N_tempM[0] * N_tempM[0] + N_tempM[1] * N_tempM[1] + N_tempM[2] * N_tempM[2] + N_tempM[3] * N_tempM[3];
		if (squa_quaterion != 0) // avoid being divided by 0
		{
			GLfloat base_quaternion = sqrt(squa_quaterion);
			N_tempM[0] = N_tempM[0] / base_quaternion;
			N_tempM[1] = N_tempM[1] / base_quaternion;
			N_tempM[2] = N_tempM[2] / base_quaternion;
			N_tempM[3] = N_tempM[3] / base_quaternion;
		}
}


//==================================================================================================
// Quaternion to Rotation Matrix : generate Rotation Matrix with Given Quaternion and Position
//==================================================================================================
void QuaternionRoatationM(GLfloat Q_tempM[7], GLfloat R[16])
{
	GLfloat w = Q_tempM[0];
	GLfloat x = Q_tempM[1];
	GLfloat y = Q_tempM[2];
	GLfloat z = Q_tempM[3];
	R[0] = 1.0f - 2.0f*y*y - 2.0f*z*z; //column1 row1
	R[1] = 2.0f*x*y + 2.0f*w*z;        //....... row2
	R[2] = 2.0f*x*z - 2.0f*w*y;		   //....... row3
	R[3] = 0.0f;					   //....... row4
	R[4] = 2.0f*x*y - 2.0f*w*z;		   //column2 row1
	R[5] = 1.0f - 2.0f*x*x - 2.0f*z*z; //....... row2
	R[6] = 2.0f*y*z + 2.0f*w*x;		   //....... row3
	R[7] = 0.0f;					   //....... row4
	R[8] = 2.0f*x*z + 2.0f*w*y;		   //column3 row1
	R[9] = 2.0f*y*z - 2.0f*w*x;		   //....... row2
	R[10] = 1.0f - 2.0f*x*x - 2.0f*y*y;//....... row3
	R[11] = 0.0f;					   //....... row4
	R[12] = Q_tempM[4];				   //column4 row1
	R[13] = Q_tempM[5];			       //....... row2
	R[14] = Q_tempM[6];			       //....... row3
	R[15] = 1.0f;					   //....... row4
}

//==================================================================================================
// Euler angle to Quaternion : generate quaternion with Given Euler Angle
//==================================================================================================
void EulerToQuaternion(GLfloat E_tempM[7])
{
	GLfloat a = E_tempM[0] / 2;
	GLfloat b = E_tempM[1] / 2;
	GLfloat c = E_tempM[2] / 2;

	// // put each value 1 position afterward to make room for the change of 3 Euler Angle turned to 4 Quaternion
	E_tempM[6] = E_tempM[5];
	E_tempM[5] = E_tempM[4];
	E_tempM[4] = E_tempM[3];
	E_tempM[0] = cos(c)*cos(b)*cos(c) + sin(c)*sin(b)*sin(a); //w
	E_tempM[1] = sin(c)*cos(b)*cos(c) - cos(c)*sin(b)*sin(a); //x
	E_tempM[2] = cos(c)*sin(b)*cos(c) + sin(c)*cos(b)*sin(a); //y
	E_tempM[3] = cos(c)*cos(b)*sin(c) - sin(c)*sin(b)*cos(a); //z
}


//=========================================================================================================
// Quaternion Interpolating Function : generate interpolation with given quaterions, postions and spline styles 
//=========================================================================================================
void q_interpolate(GLfloat p_quaternion[6][7], GLfloat SplineM[16])
{
	// Set up T matrix T = {t*t*t,t*t,t,1}
	GLfloat TMatrix_q[4] = { t*t*t, t*t, t, 1 };

	// Set up temporate matrix to store the interpolation track (pose and position)
	GLfloat tempM[7];

	// Loop to generate the interpolation track based on 4 points every time
	// i = 0, get w's changing along the Q(t) curve to tempM[0];
	// i = 1, get x's changing along the Q(t) curve to tempM[1];
	// i = 2, get y's changing along the Q(t) curve to tempM[2];
	// i = 3, get z's changing along the Q(t) curve to tempM[3];
	// i = 4, get x's changing along the Q(t) curve to tempM[4];
	// i = 5, get y's changing along the Q(t) curve to tempM[5];
	// i = 6, get z's changing along the Q(t) curve to tempM[6];

	for (int i = 0; i < 7; i++)
	{
		// the value of points would be changed by timer function in the following
		GLfloat GMatrix_q[4] = { p_quaternion[points][i],
								 p_quaternion[(points + 1)][i],
								 p_quaternion[(points + 2)][i],
								 p_quaternion[(points + 3)][i]};

		tempM[i] = blend(TMatrix_q, SplineM, GMatrix_q);
	}

	Normalization(tempM);
	QuaternionRoatationM(tempM,M);
}

//=========================================================================================================
// Euler Angle Interpolating Function : generate interpolation with given Euler angles, postions and spline styles 
//=========================================================================================================
void e_interpolate(GLfloat p_euler[7][6], GLfloat SplineM[16])
{
// Set up T matrix T = {t*t*t,t*t,t,1}
GLfloat TMatrix_e[4] = { t*t*t, t*t, t, 1 };

// Set up temporate matrix to store the interpolation track (pose and position)
GLfloat tempM[7];

// Loop to generate the interpolation track based on 4 points every time
// i = 0, get w's changing along the Q(t) curve to tempM[0];
// i = 1, get x's changing along the Q(t) curve to tempM[1];
// i = 2, get y's changing along the Q(t) curve to tempM[2];
// i = 3, get z's changing along the Q(t) curve to tempM[3];
// i = 4, get x's changing along the Q(t) curve to tempM[4];
// i = 5, get y's changing along the Q(t) curve to tempM[5];
// i = 6, get z's changing along the Q(t) curve to tempM[6];

	for (int i = 0; i < 7; i++)
	{
		// the value of points would be changed by timer function in the following
		GLfloat GMatrix_e[4] = { p_euler[points][i],
								 p_euler[(points + 1)][i],
								 p_euler[(points + 2)][i],
								 p_euler[(points + 3)][i] };

		tempM[i] = blend(TMatrix_e, SplineM, GMatrix_e);
	}

	EulerToQuaternion(tempM);
	Normalization(tempM);
	QuaternionRoatationM(tempM,M);
}


//=========================================================================================
// teapot animation : please change the notation position to get different animation effect
//=========================================================================================
void teapotAnimation()
{
	q_interpolate(point_quaternion,CRSplineM);
	//q_interpolate(point_quaternion,BSplineM);
	//e_interpolate(point_euler, CRSplineM);
	//e_interpolate(point_euler, BSplineM);
	glLoadMatrixf(M);
	
	// render objects
	glutSolidTeapot(1.0);

}
//================================
// timer : triggered every 16ms ( about 60 frames per second )
//================================
void timer(int value) {
	// render
	glutPostRedisplay();

	// Set time increase by 0.01, changing the value of points from 0 to 2
	t = t + 0.01;
	if (t >= 1)
	{
		t = 0;
		if (points < number - 4)
		{
			points++;
		}
		else
		{
			points = 0;
		}
	}
	// reset timer
	glutTimerFunc(16, timer, 0);
}


//================================
// render
//================================
void render(void) {
	// clear buffer
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render state
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// light source attributes
	GLfloat LightAmbient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat LightDiffuse[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat LightSpecular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat LightPosition[] = { 5.0f, 5.0f, 5.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

	// surface material attributes
	GLfloat material_Ka[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	GLfloat material_Kd[] = { 0.43f, 0.47f, 0.54f, 1.0f };
	GLfloat material_Ks[] = { 0.33f, 0.33f, 0.52f, 1.0f };
	GLfloat material_Ke[] = { 0.1f, 0.0f, 0.1f, 1.0f };
	GLfloat material_Se = 10;

	glMaterialfv(GL_FRONT, GL_AMBIENT, material_Ka);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_Kd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_Ks);
	glMaterialfv(GL_FRONT, GL_EMISSION, material_Ke);
	glMaterialf(GL_FRONT, GL_SHININESS, material_Se);


	// modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	// animation
	teapotAnimation();	

	// disable lighting
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

	// swap back and front buffers
	glutSwapBuffers();
}

//================================
// keyboard input
//================================
void keyboard(unsigned char key, int x, int y) {}

//================================
// reshape : update viewport and projection matrix when the window is resized
//================================
void reshape(int w, int h) {
	// screen size
	g_screenWidth = w;
	g_screenHeight = h;

	// viewport
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	// projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w / (GLfloat)h, 1.0, 2000.0);
}

//================================
// main
//================================
int main(int argc, char** argv) {
	// create opengL window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1000, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Lab 1 - Computer Animation - Fei Yan");

	// set callback functions
	glutDisplayFunc(render);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(16, timer, 0);

	// main loop
	glutMainLoop();

	return 0;
}