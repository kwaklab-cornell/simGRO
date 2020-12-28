#include <GLUT/glut.h>
#include <cstdlib>
#include <sstream>
#include "simGRO2.h"

using namespace std;

#define WIDTH 1200
#define HEIGHT 600
#define MAXOBJ 10000

// The time in milliseconds between timer ticks
#define TIMERMSECS 20

// Global variables for measuring time (in milli-seconds)
int startTime;
int prevTime;
int elapsedTime;
int timeSincePrevFrame;
float playSpeed = 5;
int saveMode = 0;
int saveFrame = 10000000;
GLubyte saveBuffer[2*WIDTH*HEIGHT*3];
ofstream saveFilelist;

stringstream title("Default");

// ----- Function Prototypes -----
static void init();
static void init1();
static void init2();
static void reshape(GLsizei w, GLsizei h);
static void animate(int value);
static void render();
static void render1();
static void render2();
static void renderall();
static void key(unsigned char k, int x, int y);
static void createmenus();
static void processmenuevents(int option);
void drawpol2(int age);
void drawIC(int age);
void drawPF(int age);
void drawDNA(int width);

// Global variables for defining objects
GLuint dl_pol2, dl_pf, dl_ic, dl_dna;
int viewDistance = 2400;

// Global variables for molecular objects
extern DNAstat DNArep;
static int objCount;
static int objType[MAXOBJ];
static int objAge[MAXOBJ];
static int objPosition[MAXOBJ];

// Windows
int mainWindow, subWindow1, subWindow2;

// ---- Function Implementations -----
int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// Set up the GLUT window
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (WIDTH, HEIGHT);
	mainWindow = glutCreateWindow ("GROgu v0.5");
	init();

	// Set up the callbacks
	glutKeyboardFunc(key);
	glutReshapeFunc(reshape);
	glutDisplayFunc(render);
	glutPostRedisplay();
	glutIdleFunc(renderall);

	// Start the timer
	glutTimerFunc(TIMERMSECS, animate, 0);

	// Initialize the time variables
	startTime = glutGet(GLUT_ELAPSED_TIME);
	prevTime = startTime;

	// Sub-window1 for animation
	subWindow1 = glutCreateSubWindow(mainWindow,0,0,WIDTH,HEIGHT/3);
	glutDisplayFunc(render1);
	glutKeyboardFunc(key);
	init1();

	// Sub-window2 for pol2 profile
	subWindow2 = glutCreateSubWindow(mainWindow,0,HEIGHT/3,WIDTH,HEIGHT*2/3);
	glutDisplayFunc(render2);
	glutKeyboardFunc(key);
	init2();

	// Start the main loop
	reshape(WIDTH,HEIGHT);
	glutMainLoop();
	return 0;
}

static void init()
{
	InitSimGRO();
	createmenus();
//	initmypng(WIDTH, HEIGHT);
}

static void init1()
{

	glClearColor(0, 0, 0, 0);
	glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	

	// Create Pol2 display lists
	dl_pol2 = glGenLists(8);
	for(int i=0;i<8;i++)
	{
		glNewList(dl_pol2+i, GL_COMPILE);
		drawpol2(i);
		glEndList();
	}
	
	// Pausing factors
	dl_pf = glGenLists(8);
	for(int i=0;i<8;i++)
	{
		glNewList(dl_pf+i, GL_COMPILE);
		drawPF(i);
		glEndList();
	}
	// Initiation complex
	dl_ic = glGenLists(8);
	for(int i=0;i<8;i++)
	{
		glNewList(dl_ic+i, GL_COMPILE);
		drawIC(i);
		glEndList();
	}
	// DNAs
	dl_dna = glGenLists(3);
	glNewList(dl_dna, GL_COMPILE);
	drawDNA(1);
    glEndList();
	glNewList(dl_dna+1, GL_COMPILE);
	drawDNA(2);
    glEndList();
	glNewList(dl_dna+2, GL_COMPILE);
	drawDNA(3);
    glEndList();
	createmenus();
}

static void init2()
{
	glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glClearColor(1, 1, 1, 0);
	createmenus();
}

static void reshape (GLsizei w, GLsizei h)
{
	// Reshape sub-window1
	glutSetWindow(subWindow1);
	glutPositionWindow(0,0);
	glutReshapeWindow(w,h/3);

	glViewport (0, 0, (GLsizei) w, (GLsizei) h/3);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	gluPerspective (5.5,(GLfloat)w/(GLfloat)h*3, 1, 10000);
//	glFrustum(0, (GLdouble) w, -(GLdouble) h/8, (GLdouble) h/8);
	glMatrixMode(GL_MODELVIEW);

	// Reshape sub-window2
	glutSetWindow(subWindow2);
	glutPositionWindow(0,h/3);
	glutReshapeWindow(w,h*2/3);
	
	glViewport (0, -(GLsizei) h*3/8, (GLsizei) w, (GLsizei) h*3/4);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	gluPerspective (5.5,(GLfloat)w/(GLfloat)h*3, 1, 10000);
//	gluOrtho2D(0, (GLdouble) w, -(GLdouble) h*3/8, (GLdouble) h*3/8);
	glMatrixMode(GL_MODELVIEW);

}

void saveScreen()
{
	stringstream savefn;
	char sfn[256];
	savefn.str("");
	savefn<<"screencap/sc"<<saveFrame<<".png";
	savefn>>sfn;
//	writemypng(sfn,WIDTH,HEIGHT,saveBuffer);
	saveFilelist<<"sc"<<saveFrame<<".png"<<endl;
}

static void animate(int value)
{
	// Set up the next timer tick (do this first)
	glutTimerFunc(TIMERMSECS, animate, 0);
	// Measure the elapsed time
	int currTime = glutGet(GLUT_ELAPSED_TIME);
	timeSincePrevFrame = currTime - prevTime;
	elapsedTime = currTime - startTime;
	if(saveMode) // if at the saveMode
	{
		UpdateSimGRO(TIMERMSECS*playSpeed);
		saveScreen();
		saveFrame++;
	}
	else if(timeSincePrevFrame*playSpeed>1000) UpdateSimGRO(999);
	else UpdateSimGRO(timeSincePrevFrame*playSpeed);

	Print(elapsedTime,timeSincePrevFrame, saveFrame);

	// Force a redisplay to render the new image
	glutPostRedisplay();

	prevTime = currTime;
}

static void render()
{
	glutSetWindow(mainWindow);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glutSwapBuffers();
}

void rendertext(stringstream &ss, float x, float y,float aspect=1,float red=0, float green=0, float blue=0)
{
	glPushMatrix();
	glRasterPos2f(x/100*viewDistance,y/100*viewDistance);
	int ss_size = ss.str().size();
	glColor3f(red,green,blue);
	for (int i=0;i<ss_size;i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, ss.str().at(i));
	glPopMatrix();
}

void gridtest()
{
    stringstream ss;
    int ss_size;
    for(int i=-50;i<50;i++)
        for(int j=-50;j<50;j++)
        {
            glRasterPos2f((float)i/100*viewDistance,(float)j/100*viewDistance);
            ss.str("");
            ss<<i<<":"<<j;
            ss_size=ss.str().size();
            for(int k=0;k<ss_size;k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ss.str().at(k));
        }
}


static void render1()
{

	// Clear the screen
	glutSetWindow(subWindow1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Viewing transformation
	glLoadIdentity();	
	gluLookAt(viewDistance*0.23,0,viewDistance,viewDistance*0.23,viewDistance*0.01,0,0,1,0);

	// Labels
	glColor3f(1,0.5,0);
//    gridtest();
    glColor3f(1,1,1);
	rendertext(title,3,5,2,1,1,1);
	stringstream ss;
	ss<<"TF";
	rendertext(ss,-3,4,2,0.4,0.3,1);
	ss.str("");
	ss<<"on once per "<<IC_on_hl<<"s";
	rendertext(ss,-1,4,2,1,1,1);
	ss.str("");
	ss<<"TF";
	rendertext(ss,-3,3,2,0.4,0.3,1);
	ss.str("");
	ss<<"off after "<<IC_off_hl<<"s on average";
	rendertext(ss,-1,3,2,1,1,1);
	ss.str("");
	ss<<"Pol II";
	rendertext(ss,-3,-2,2,1,0.3,0.4);
	ss.str("");
	ss<<"escape once per "<<PF_off_hl<<"s";
	rendertext(ss,0,-2,2,1,1,1);
	
	
	// Draw Objects
	DNArep.ObjectList(objCount, objType, objAge, objPosition);
	if(objCount>MAXOBJ) objCount = MAXOBJ;

	// Draw DNA
	glPushMatrix();
	if(viewDistance>2000) glCallList(dl_dna);
	else if(viewDistance>1000) glCallList(dl_dna+1);
	else glCallList(dl_dna+2);
	glPopMatrix();
	
	for(int i=0;i<objCount;i++)
	{
		switch(objType[i])
		{
			case IC:
				glPushMatrix();
				glTranslatef(objPosition[i],0,0);
				if(objAge[i]>0) glCallList(dl_ic+objAge[i]-1);
				else glCallList(dl_ic-objAge[i]-1);
				glPopMatrix();
				break;
			case POL2:
				// Draw Pol2 body
				glPushMatrix();
				glTranslatef(objPosition[i],0,0);
				if(objAge[i]>0) glCallList(dl_pol2+objAge[i]-1);
				else glCallList(dl_pol2-objAge[i]-1);
				glPopMatrix();
				break;
			case POL2PF:
				glPushMatrix();
				glTranslatef(objPosition[i]+10,10,0);
				if(objAge[i]>0) glCallList(dl_pf+objAge[i]-1);
				else glCallList(dl_pf-objAge[i]-1);
				glPopMatrix();
				break;
		}
	}
	


	// Swap the buffers (if double-buffered) to show the rendered image
	glutSwapBuffers();
	// Save screen
	if(saveMode)
	{
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, saveBuffer+WIDTH*3*HEIGHT*2/3);
	}
	glFlush();
}

int scale_f=1;
float scale_e=1;

	

void renderlabel(float axisval, int position)
{
	stringstream ss;
    glPushMatrix();
	glRasterPos2f(-0.03*viewDistance,(float)position/100*viewDistance*0.07);
    ss<<axisval<<"%";
	int ss_size=ss.str().size();
	for (int i=0;i<ss_size;i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, ss.str().at(i));
	glPopMatrix();
}

static void render2()
{
	// Render the screen
	glutSetWindow(subWindow2);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Viewing transformation
	glLoadIdentity();	
	gluLookAt(viewDistance*0.23,viewDistance*0.033,viewDistance,viewDistance*0.23,viewDistance*0.033,0,0,1,0);

	// Draw horizontal gridlines
	glColor3f(0.7,0.7,0.7);
	glLineWidth(1);
	glBegin(GL_LINES);
	glVertex2f(0,viewDistance*0.07);
	glVertex2f(DNALEN,viewDistance*0.07);
	if(scale_f==1)
	{
		glVertex2f(0,viewDistance*0.07*0.5);
		glVertex2f(DNALEN,viewDistance*0.07*0.5);
		glVertex2f(0,viewDistance*0.07*0.2);
		glVertex2f(DNALEN,viewDistance*0.07*0.2);
	}
	else if(scale_f==5)
	{
		glVertex2f(0,viewDistance*0.07*0.4);
		glVertex2f(DNALEN,viewDistance*0.07*0.4);
		glVertex2f(0,viewDistance*0.07*0.2);
		glVertex2f(DNALEN,viewDistance*0.07*0.2);
	}
	else
	{
		glVertex2f(0,viewDistance*0.07*0.5);
		glVertex2f(DNALEN,viewDistance*0.07*0.5);
		glVertex2f(0,viewDistance*0.07*0.25);
		glVertex2f(DNALEN,viewDistance*0.07*0.25);
	}

	// Draw vertical gridline
	glVertex2f(DNALEN,0);
	glVertex2f(DNALEN,viewDistance*0.07);
	glEnd();

	// Draw TSS and PAS site grids
	glLineWidth(2);
	glBegin(GL_LINES);
	glVertex2f(TSSPOS,0);
	glVertex2f(TSSPOS,viewDistance*0.07);
	glVertex2f(PASPOS,0);
	glVertex2f(PASPOS,viewDistance*0.07);
	glEnd();

	
	
	//Draw profile
	glColor3f(.9,0,0.1);
	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	for(int i=0;i<DNALEN;i++) glVertex2f(i,Metagene[i]/(scale_e*scale_f)*viewDistance*0.07);
	glEnd();

	// Draw Axis
	glColor3f(0,0,0);
	glLineWidth(3);
	glBegin(GL_LINES);
	glVertex2f(0,0);
	glVertex2f(DNALEN,0);
	glEnd();

	//Draw axis
	glColor3f(0,0,0);
	glLineWidth(3);
	glBegin(GL_LINES);
	glVertex2f(0,viewDistance*0.07);
	glVertex2f(0,0);
	glEnd();
	
	// Labels
	float axisval=scale_e*scale_f*100;
	if(scale_f==1)
	{
		renderlabel(axisval,100);
		renderlabel(axisval*0.5,50);
		renderlabel(axisval*0.2,20);
	}
	else if(scale_f==5)
	{
		renderlabel(axisval,100);
		renderlabel(axisval*0.4,40);
		renderlabel(axisval*0.2,20);
	}
	else
	{
		renderlabel(axisval,100);
		renderlabel(axisval*0.5,50);
		renderlabel(axisval*0.25,25);
	}

	glLineWidth(2);
	stringstream ss;
	ss.str("TSS");
	rendertext(ss,(float)(TSSPOS-20)/viewDistance*100,-0.5);
	ss.str("3' end (Poly-A site)");
	rendertext(ss,(float)(PASPOS-20)/viewDistance*100,-0.5);
	ss.str("Real-time Pol II profile");
	rendertext(ss,17,-1);
	switch(Assay_res_type)
	{
		case 0:
			ss.str("2 bp resolution");
			rendertext(ss,45,-1);
			break;
		case 1:
			ss.str("10 bp resolution");
			rendertext(ss,45,-1);
			break;
		case 2:
			ss.str("50 bp resolution");
			rendertext(ss,45,-1);
			break;
	}
			
	// Swap the buffers (if double-buffered) to show the rendered image
	glutSwapBuffers();
	if(saveMode)
	{
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, saveBuffer);
	}	
	glFlush();
}

static void renderall()
{
	render1();
	render2();
}

int pref_n;
char **pref_filename;
char **pref_description;

static void key(unsigned char k, int x, int y)
{
	if(k>='1' && k-'1'<pref_n)
	{
		ResetSimGRO(pref_filename[k-'1']);
		title.str(pref_description[k-'1']);
		return;
	}
	else switch (k) {
		case 27:  /* Escape */
			exit(0); // quit
			break;
		case 'h':
			if(saveMode==0)
			{
				saveMode=1;
				saveFilelist.open("screencap/scList.txt");
			}
			break;
		case 'j':
			if(saveMode==1)
			{
				saveMode=0;
				saveFrame=10000000;
				saveFilelist.close();
			}
			break;
		case 'a':
			viewDistance/=1.05;
			if(viewDistance<DNArep.TSS*4) viewDistance*=1.05;
			break;
		case 's':
			viewDistance*=1.05;
			if(viewDistance>DNArep.Length*2) viewDistance/=1.05;
			break;	
		case 'q':
			viewDistance=DNArep.TSS*4;
			break;
		case 'w':
			viewDistance=2400;
			break;	
		case 'e':
			viewDistance=DNArep.Length*2;
			break;
		case 'f':
			DNArep.Promoter_status[IC]=1;
			break;
		case 'g':
			DNArep.Promoter_status[IC]=-8;
			break;
		case 'z':
			if(scale_f==1)
			{
				scale_f=5;
				scale_e/=10;
			}
			else if(scale_f==5) scale_f=2;
			else scale_f=1;
			break;
		case 'x':
			if(scale_f==5)
			{
				scale_f=1;
				scale_e*=10;
			}
			else if(scale_f==2) scale_f=5;
			else scale_f=2;
			break;
		case 'c':
			if(playSpeed == 1) playSpeed = 10;
			else playSpeed=1;
			cout<<"Play speed = "<<playSpeed<<endl;
			break;
		case 'v':
			for(list<Pol2stat>::iterator it = DNArep.Pol2list.begin();it!=DNArep.Pol2list.end();it++)
			{
				if(it->Status[PAUSED])
				{
					  it->Status[PFBOUND]=-8;
					  it->Status[PAUSED]=0;
					  it->Status[ELONGATING]=1;
				}
			}
			break;
		case 'r':
			if(POP_N==500) POP_N=1000;
			else if(POP_N==1000) POP_N=2000;
			else if(POP_N==2000) POP_N=5000;
			else if(POP_N==5000) POP_N=10000;
			else if(POP_N==10000) POP_N=20000;
			cout<<"Simulation depth = "<<POP_N<<" copies of DNA templates"<<endl;
			break;
		case 't':
			if(POP_N==1000) POP_N=500;
			else if(POP_N==2000) POP_N=1000;
			else if(POP_N==5000) POP_N=2000;
			else if(POP_N==10000) POP_N=5000;
			else if(POP_N==20000) POP_N=10000;
			cout<<"Simulation depth = "<<POP_N<<" copies of DNA templates"<<endl;
			break;
		default:
			return;
	}
	// Force a redraw of the screen in order to update the display
	glutPostRedisplay();
}


static void createmenus()
{
	ifstream preflist("gu/pref.txt");
	preflist>>pref_n;
	preflist.ignore(1024,'\n');
	pref_filename = new char*[pref_n];
	pref_description = new char*[pref_n];
	for(int i=0;i<pref_n;i++)
	{
		pref_filename[i] = new char[1024];
		pref_description[i] = new char[1024];
		preflist.getline(pref_description[i],1024);
	}
	for(int i=0;i<pref_n;i++) preflist>>pref_filename[i];
	preflist.close();

 	int menu,submenu1, submenu2, submenu3, submenu4, submenu5;

	submenu1 = glutCreateMenu(processmenuevents);
	glutAddMenuEntry("Default",0);
	for(int i=0;i<pref_n;i++) glutAddMenuEntry(pref_description[i],(i+1)*100);

	submenu2 = glutCreateMenu(processmenuevents);
	glutAddMenuEntry("Near",1);
	glutAddMenuEntry("Mid",2);
	glutAddMenuEntry("Far",3);

	submenu3 = glutCreateMenu(processmenuevents);
	glutAddMenuEntry("x 1",11);
	glutAddMenuEntry("x 2",12);
	glutAddMenuEntry("x 5",13);
	glutAddMenuEntry("x 10",14);
	
	submenu4 = glutCreateMenu(processmenuevents);
	glutAddMenuEntry("2 bp",21);
	glutAddMenuEntry("10 bp",22);
	glutAddMenuEntry("50 bp",23);
	
	submenu5 = glutCreateMenu(processmenuevents);
	glutAddMenuEntry("500",31);
	glutAddMenuEntry("1000",32);
	glutAddMenuEntry("2000",33);
	glutAddMenuEntry("5000",34);
	glutAddMenuEntry("10000",35);
	glutAddMenuEntry("20000",36);

	menu = glutCreateMenu(processmenuevents);
	glutAddSubMenu("Load parameters",submenu1);
	glutAddSubMenu("Zoom",submenu2);
	glutAddSubMenu("Play speed",submenu3);
	glutAddSubMenu("Assay resolution",submenu4);
	glutAddSubMenu("Simulation depth",submenu5);
	glutAddMenuEntry("Increase scale",5);
	glutAddMenuEntry("Decrease scale",6);
	
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

static void processmenuevents(int option)
{
	switch(option)
	{
		case 0:
			ResetSimGRO();
			break;
		case 1:
			viewDistance = DNArep.TSS*4;
			break;
		case 2:
			viewDistance = 2400;
			break;
		case 3:
			viewDistance = DNArep.Length*2;
			break;
		case 4:
			PrintStatus(elapsedTime,timeSincePrevFrame);
			break;
		case 5:
			if(scale_f==5)
			{
				scale_f=1;
				scale_e*=10;
			}
			else if(scale_f==2) scale_f=5;
			else scale_f=2;
			break;
		case 6:
			if(scale_f==1)
			{
				scale_f=5;
				scale_e/=10;
			}
			else if(scale_f==5) scale_f=2;
			else scale_f=1;
			break;
		case 11:
			playSpeed = 1;
			break;
		case 12:
			playSpeed = 2;
			break;
		case 13:
			playSpeed = 5;
			break;
		case 14:
			playSpeed = 10;
			break;
		case 21:
			Assay_res_type=0;
			break;
		case 22:
			Assay_res_type=1;
			break;
		case 23:
			Assay_res_type=2;
			break;
		case 31:
			POP_N = 500;
			break;
		case 32:
			POP_N = 1000;
			break;
		case 33:
			POP_N = 2000;
			break;
		case 34:
			POP_N = 5000;
			break;
		case 35:
			POP_N = 10000;
			break;
		case 36:
			POP_N = 20000;
			break;
	}	
	if(option>=100)
	{
		ResetSimGRO(pref_filename[option/100-1]);
		title.str(pref_description[option/100-1]);
	}
}

void drawpol2(int age)
{
	glColor4f(.9,.9,.9,(float)age/7);
	glBegin(GL_POLYGON);
    glVertex2f(-20,-11);
    glVertex2f(-19.2,-12.8);
    glVertex2f(-17.8,-14.2);
    glVertex2f(-16,-15);
    glVertex2f(-1,-15);
    glVertex2f(0.8,-14.2);
    glVertex2f(12.2,-2.8);
    glVertex2f(13,-1);
    glVertex2f(13,1);
    glVertex2f(12.2,2.8);
    glVertex2f(0.8,14.2);
    glVertex2f(-1,15);
    glVertex2f(-16,15);
    glVertex2f(-17.8,14.2);
    glVertex2f(-19.2,12.8);
    glVertex2f(-20,11);
//    glVertex2f(-20,11);
    glEnd();
	glColor4f(0.9,0,0.1,(float)age/7);
	glBegin(GL_POLYGON);
    glVertex2f(-20+2,11);
    glVertex2f(-19.2+1.4,12.8-1.4);
    glVertex2f(-17.8+1.4,14.2-1.4);
    glVertex2f(-16,15-2);
    glVertex2f(-1,15-2);
    glVertex2f(0.8-1.4,14.2-1.4);
    glVertex2f(12.2-1.4,2.8-1.4);
    glVertex2f(13-2,1);
    glVertex2f(13-2,-1);
    glVertex2f(12.2-1.4,-2.8+1.4);
    glVertex2f(0.8-1.4,-14.2+1.4);
    glVertex2f(-1,-15+2);
    glVertex2f(-16,-15+2);
    glVertex2f(-17.8+1.4,-14.2+1.4);
    glVertex2f(-19.2+1.4,-12.8+1.4);
    glVertex2f(-20+2,-11);
//    glVertex2f(-20,11);
    glEnd();
}

void drawIC(int age)
{
	glColor4f(0.9,0.9,0.9,(float)age/7);
	glBegin(GL_POLYGON);
    glVertex2f(30,16);
    glVertex2f(29.2,17.8);
    glVertex2f(27.8,19.2);
    glVertex2f(26,20);
    glVertex2f(-26,20);
    glVertex2f(-27.8,19.2);
    glVertex2f(-29.2,17.8);
    glVertex2f(-30,16);
    glVertex2f(-30,-1);
    glVertex2f(-29.2,-2.8);
    glVertex2f(-27.8,-4.2);
    glVertex2f(-26,-5);
    glVertex2f(26,-5);
    glVertex2f(27.8,-4.2);
    glVertex2f(29.2,-2.8);
    glVertex2f(30,-1);
//    glVertex2f(-30,16);
    glEnd();
	glColor4f(0.1,0,0.9,(float)age/7);
	glBegin(GL_POLYGON);
	glVertex2f(-30+2,16);
    glVertex2f(-29.2+1.4,17.8-1.4);
    glVertex2f(-27.8+1.4,19.2-1.4);
    glVertex2f(-26,20-2);
    glVertex2f(26,20-2);
    glVertex2f(27.8-1.4,19.2-1.4);
    glVertex2f(29.2-1.4,17.8-1.4);
    glVertex2f(30-2,16);
    glVertex2f(30-2,-1);
    glVertex2f(29.2-1.4,-2.8+1.4);
    glVertex2f(27.8-1.4,-4.2+1.4);
    glVertex2f(26,-5+2);
    glVertex2f(-26,-5+2);
    glVertex2f(-27.8+1.4,-4.2+1.4);
    glVertex2f(-29.2+1.4,-2.8+1.4);
    glVertex2f(-30+2,-1);
//    glVertex2f(-30,16);
    glEnd();
}

void drawPF(int age)
{
	glColor4f(0.9,0.9,0.9,(float)age/7);
	glBegin(GL_POLYGON);
    for(int i=0;i<32;i++)
    {
        glVertex2f(12*cos((float)i/32*2*3.141593),12*sin((float)i/32*2*3.141593));
    }
	glEnd();
    glColor4f(0.9,0.9,0,(float)age/7);
	glBegin(GL_POLYGON);
	for(int i=0;i<32;i++)
    {
        glVertex2f(10*cos((float)i/32*2*3.141593),10*sin((float)i/32*2*3.141593));
    }
    glEnd();
}

void drawDNA(int width)
{
	glColor3f(0.7,0.7,0.7);
    glLineWidth(width);
	glBegin(GL_LINE_STRIP);
	for(int i=0;i<DNArep.Length;i++)
	{
		glVertex2f(i,3*sin((float)i*2*3.141593/10));
	}
    glEnd();
	glBegin(GL_LINE_STRIP);
	for(int i=0;i<DNArep.Length;i++)
	{
		glVertex2f(i,3*cos((float)i*2*3.141593/10));
	}
    glEnd();
	glColor3f(1,1,1);
    glLineWidth(width+1);
    glBegin(GL_LINES);
    glVertex2f(TSSPOS,10);
    glVertex2f(TSSPOS,30);
    glVertex2f(TSSPOS,30);
    glVertex2f(TSSPOS+20,30);
    glVertex2f(TSSPOS+10,40);
    glVertex2f(TSSPOS+20,30);
    glVertex2f(TSSPOS+10,20);
    glVertex2f(TSSPOS+20,30);
	glVertex2f(PASPOS-2,10);
	glVertex2f(PASPOS-2,30);
	glVertex2f(PASPOS+2,10);
	glVertex2f(PASPOS+2,30);
	glEnd();
}

