#include<iostream>//Lighting implemented
#include<GL/gl.h>
#include<GL/glut.h>
#include<GL/glu.h>
#include <stdlib.h>
#include<string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include<cmath>
#include "ifs.h"

#define UINTSIZE sizeof(unsigned int)
#define FLOATSIZE sizeof(float)
using namespace std;
IFS_DATA* MODEL;//global.Make changes to this struct
float dx=0,dy=0,dz=0;//translation steps variables.
float totx=0,toty=0,totz=0;//used in rotation about object's own axis.Pass these to all functions so that totx,toty,totz get updated
float scale=1;//scaling state variable
int thetax=0,thetay=0,thetaz=0;//rotation state variables in degrees along the 3 world axes
float eyex=0,eyey=0,eyez=5;
float targx=0,targy=0,targz=0;
float upx=0,upy=1,upz=0;
float alphax=0,alphay=0,alphaz=0;//rotation along 3 object axes
float delta1=0;//rotation angle about current look vector
float delta2=0;
float delta3=0;
const double PI=acos(-1.0);
float forward[3], side[3], up[3];
    
    
ssize_t read_uint32(int infd, unsigned int* uint_star) {
    if (read(infd, uint_star, UINTSIZE) == UINTSIZE) {
	return UINTSIZE;
    } else {
	fprintf(stderr, "Error reading size of a uint32 field\n");
	exit(-1);
    }
}

ssize_t read_float32(int infd, float* float_star) {
    if (read(infd, float_star, FLOATSIZE) == FLOATSIZE) {
        return FLOATSIZE;
    } else {
        fprintf(stderr, "Error reading size of a float32 field\n");
        exit(-1);
    }
}

ssize_t read_string32 (int infd, char** buf) {
    unsigned int mem_len = 0;
    read_uint32(infd, &mem_len);
    void *ptr;
    ptr= realloc(*buf, mem_len);
    *buf=(char*)ptr;
    if (mem_len == read(infd, *buf, mem_len)) {
	return mem_len;
    }
    else {
	fprintf(stderr, "Error reading a string32 field\n");
	exit(-1);
    }
}

IFS_DATA* load_ifs_file (const char* filename) {
    int infd;
    IFS_DATA* ifs_data = NULL;
    float version;
    char* ifstag = NULL;
    unsigned int i;
    unsigned int nVertices = 0;
    unsigned int nTriangles = 0;
    unsigned int tmp_Index = 0;
    
    if ((infd = open(filename, O_RDONLY)) < 2) {
	fprintf(stderr, "Error opening an input IFS file\n");
	exit(-1);
    }

    ifs_data = (IFS_DATA*) malloc(sizeof(IFS_DATA));
    ifs_data->modelName = NULL;
    ifs_data->numVertices = 0;
    ifs_data->vertices = NULL;
    ifs_data->numTriangles = 0;
    ifs_data->triangles = NULL;
    
    read_string32(infd, &ifstag);
    if (strcmp(ifstag, "IFS") != 0) {
	fprintf(stderr, "Not IFS filetype\n");
	exit(-1);
    }
    free(ifstag); ifstag = NULL;

    read_float32(infd, &version);
    if (version != 1.0) {
	fprintf(stderr, "Invalid version number: %f\n", version);
	exit(-1);
    }

    read_string32(infd, &(ifs_data->modelName));
    
    read_string32(infd, &ifstag);
    if (strcmp(ifstag, "VERTICES") != 0) {
	fprintf(stderr, "Not IFS filetype\n");
	exit(-1);
    }
    free(ifstag); ifstag = NULL;

    read_uint32(infd, &nVertices);
    ifs_data->numVertices = nVertices;

    ifs_data->vertices = (VERTEX*) malloc(nVertices * sizeof(VERTEX));
    for (i =0; i < ifs_data->numVertices; ++i) {
	ifs_data->vertices[i].id = i;
	read_float32(infd, &((ifs_data->vertices)[i].x));
	read_float32(infd, &((ifs_data->vertices)[i].y));
	read_float32(infd, &((ifs_data->vertices)[i].z));
    }
    
    read_string32(infd, &ifstag);
    if (strcmp(ifstag, "TRIANGLES") != 0) {
	fprintf(stderr, "Not IFS filetype\n");
	exit(-1);
    }
    free(ifstag); ifstag = NULL;
	
    read_uint32(infd, &nTriangles);
    ifs_data->numTriangles = nTriangles;
    
    ifs_data->triangles = (TRIANGLE*) malloc(nTriangles * sizeof(TRIANGLE));
    for (i =0; i < ifs_data->numTriangles; ++i) {
	read_uint32(infd, &tmp_Index);
	if (tmp_Index >= nVertices) {
	    fprintf(stderr, "Invalid Vertex index\n");
	    exit(-1);
	}
	ifs_data->triangles[i].a = &((ifs_data->vertices)[tmp_Index]);
	read_uint32(infd, &tmp_Index);
	if (tmp_Index >= nVertices) {
	    fprintf(stderr, "Invalid Vertex index\n");
	    exit(-1);
	}
	ifs_data->triangles[i].b = &((ifs_data->vertices)[tmp_Index]);
	read_uint32(infd, &tmp_Index);
	if (tmp_Index >= nVertices) {
	    fprintf(stderr, "Invalid Vertex index\n");
	    exit(-1);
	}
	ifs_data->triangles[i].c = &((ifs_data->vertices)[tmp_Index]);
    }

    if (close(infd) == -1) {
	fprintf(stderr, "Error closing an input IFS file\n");
	exit(-1);
    }
	
    return ifs_data;
}

void free_ifs_data(IFS_DATA** ifs_data) {
    if (ifs_data) {
	free((*ifs_data)->modelName);
	free((*ifs_data)->vertices);
	free((*ifs_data)->triangles);
    }
    free(*ifs_data);
    *ifs_data = NULL;
}

void print_ifs_summary(FILE* target, IFS_DATA* ifs_data) {
    unsigned int i;
    fprintf(target, "=====  IFS  SUMMARY  =====\n");
    fprintf(target, " Model name          : %s\n", ifs_data->modelName);
    fprintf(target, " Number of vertices  : %d\n", ifs_data->numVertices);
    for (i=0; i<ifs_data->numVertices; ++i) {
	fprintf(target, " v_%06d : (%8f, %8f, %8f)\n",
		ifs_data->vertices[i].id,
		ifs_data->vertices[i].x, 
		ifs_data->vertices[i].y,
		ifs_data->vertices[i].z);
    }
    fprintf(target, " Number of triangles : %d\n", ifs_data->numTriangles);
    for (i=0; i<ifs_data->numTriangles; ++i) {
	fprintf(target, " t_%06d : (v_%06d, v_%06d, v_%06d)\n", i,
		(ifs_data->triangles[i].a)->id, 
		(ifs_data->triangles[i].b)->id,
		(ifs_data->triangles[i].c)->id);
    }
    fprintf(target, "===== END OF SUMMARY =====\n");
}
//-------------------------------------------------Vector class-----------------
class Vector4{
public:
double dim[4];

Vector4(){
for(int i=0;i<4;i++)dim[i]=0;
}

Vector4(float a,float b,float c,float d)
{dim[0]=a;
 dim[1]=b;
 dim[2]=c;
 dim[3]=d;}

void homogenize(){
 if(dim[3]!=1)dim[0]/=dim[3],dim[1]/=dim[3],dim[2]/=dim[3],dim[3]=1;
 }
 
Vector4& operator= (const Vector4& param)
{ for(int i=0;i<4;i++)
 	dim[i]=param.dim[i];
  return *this;
}
void print()
{cout<<dim[0]<<" "<<dim[1]<<" "<<dim[2]<<" "<<dim[3]<<endl;}

};
//---------------------------------------------------------Matrix Class-------------------------------
class Matrix4{
public:
float a[4][4];
Matrix4(){//initialized as identity matrix
for(int i=0;i<4;i++)
	{for(int j=0;j<4;j++)
		{	a[i][j]=0;	
		}
	}
	}
void LoadIdentity()
{for(int i=0;i<4;i++)
	{for(int j=0;j<4;j++)
		{if(i==j)
			a[i][j]=1;
		 else
		 	a[i][j]=0;	
		}
	}
}
void transpose()
{
 for(int i=0;i<4;i++)
 	{for(int j=i+1;j<4;j++)
 		{float temp=a[i][j];
 		 a[i][j]=a[j][i];
 		 a[j][i]=temp;}
 	}
}

Matrix4 operator *(double& scalar)
{Matrix4 res;
 for(int i=0;i<4;i++)
 	{for(int j=0;j<4;j++)
 		res.a[i][j]=scalar*a[i][j];
 	}
 return res;
}

Vector4 operator *(Vector4& vec)//Matrix4*Vector4
{Vector4 res;
 for(int i=0;i<4;i++)
 	{res.dim[i]=0;
 	 for(int j=0;j<4;j++)
 		res.dim[i]+=a[i][j]*vec.dim[j];
 	}
 res.homogenize();
 return res;}

Matrix4 operator *(Matrix4& param)//Matrix*Matrix
{Matrix4 res;
 for(int i=0;i<4;i++)//rows
 {for(int j=0;j<4;j++)
 	{
 	 for(int k=0;k<4;k++)
 	 	res.a[i][j]+=a[i][k]*param.a[k][j];
 	}
 }
return res;}

Matrix4& operator =(const Matrix4& param)   
{for(int i=0;i<4;i++)
 	{for(int j=0;j<4;j++)
 		a[i][j]=param.a[i][j];
 		}
return *this;
}

void print()
{for(int i=0;i<4;i++)
	{for(int j=0;j<4;j++)
		printf("%f ",a[i][j]);
 	 printf("\n");}
}

}CameraRotMatrix,CameraTransMatrix;//used in mylookat
//----------------------------------------------------------------------------------------------------
void Translate()//in world system
{totx+=dx;
 toty+=dy;
 totz+=dz;
  
for(int i=0;i<MODEL->numVertices;i++){
 	MODEL->vertices[i].x+=dx;
 	MODEL->vertices[i].y+=dy;
 	MODEL->vertices[i].z+=dz;
 	}
}
void Scale()
{totx*=scale;
 toty*=scale;
 totz*=scale;
 for(int i=0;i<MODEL->numVertices;i++){
 	MODEL->vertices[i].x*=scale;
        MODEL->vertices[i].y*=scale;
        MODEL->vertices[i].z*=scale;}
}
void RotateX()//Rotate about world x-axis
{float temp=toty;
 toty=temp*cosf(thetax*PI/180.0)-totz*sinf(thetax*PI/180.0);
 totz=temp*sinf(thetax*PI/180.0)+totz*cosf(thetax*PI/180.0);
for(int i=0;i<MODEL->numVertices;i++){
 	float t=MODEL->vertices[i].y;
        MODEL->vertices[i].y=t*cosf(thetax*PI/180.0)-MODEL->vertices[i].z*sinf(thetax*PI/180.0);
        MODEL->vertices[i].z=t*sinf(thetax*PI/180.0)+MODEL->vertices[i].z*cosf(thetax*PI/180.0);
        }
}
void RotateY()//Rotate about world y-axis
{float temp=totx;
totx=temp*cosf(thetay*PI/180.0)-totz*sinf(thetay*PI/180.0);
totz=temp*sinf(thetay*PI/180.0)+totz*cosf(thetay*PI/180.0);
for(int i=0;i<MODEL->numVertices;i++){
 	 float t=MODEL->vertices[i].x;
 MODEL->vertices[i].x=t*cosf(thetay*PI/180.0)-MODEL->vertices[i].z*sinf(thetay*PI/180.0);
 MODEL->vertices[i].z=t*sinf(thetay*PI/180.0)+MODEL->vertices[i].z*cosf(thetay*PI/180.0);
      }

}
void RotateZ()//Rotate about world-z axis
{float temp=totx;
totx=temp*cosf(thetaz*PI/180.0)+toty*sinf(thetaz*PI/180.0);
toty=-temp*sinf(thetaz*PI/180.0)+toty*cosf(thetaz*PI/180.0);
 for(int i=0;i<MODEL->numVertices;i++){
 	float t=MODEL->vertices[i].x;
 MODEL->vertices[i].x=t*cosf(thetaz*PI/180.0)+MODEL->vertices[i].y*sinf(thetaz*PI/180.0);
 MODEL->vertices[i].y=-t*sinf(thetaz*PI/180.0)+MODEL->vertices[i].y*cosf(thetaz*PI/180.0);
 	} 
 }

void Rotate(float angle,float x,float y,float z)//Used in rotate object about camera axis.Actually rotate about arbitrary world axis 
{
//Do adjustment for totx,toty,totz
float sinAngle, cosAngle;
 float mag = sqrtf(x * x + y * y + z * z);
 
 sinAngle = sinf ( angle * PI / 180.0f );
 cosAngle = cosf ( angle * PI / 180.0f );
 float xx, yy, zz, xy, yz, zx, xs, ys, zs;
      float oneMinusCos;
      //ESMatrix rotMat;
 
      x /= mag;
      y /= mag;
      z /= mag;
 
      xx = x * x;
      yy = y * y;
      zz = z * z;
      xy = x * y;
      yz = y * z;
      zx = z * x;
      xs = x * sinAngle;
      ys = y * sinAngle;
      zs = z * sinAngle;
      oneMinusCos = 1.0f - cosAngle;
      
 float tempx,tempy,tempz;
 tempx=totx*((oneMinusCos * xx) + cosAngle)+toty*((oneMinusCos * xy) - zs)+totz*((oneMinusCos * zx) + ys);
 tempy=totx*((oneMinusCos * xy) + zs)+toty*((oneMinusCos * yy) + cosAngle)+totz*((oneMinusCos * yz) - xs);
 tempz=totx*((oneMinusCos * zx) - ys)+toty*((oneMinusCos * yz) + xs)+totz*((oneMinusCos * zz) + cosAngle);
 
 for(int i=0;i<MODEL->numVertices;i++){
  if ( mag > 0.0f )
   {
      float resx,resy,resz;
      resx=MODEL->vertices[i].x*((oneMinusCos * xx) + cosAngle)+MODEL->vertices[i].y*((oneMinusCos * xy) - zs)+MODEL->vertices[i].z*((oneMinusCos * zx) + ys);
      resy=MODEL->vertices[i].x*((oneMinusCos * xy) + zs)+MODEL->vertices[i].y*((oneMinusCos * yy) + cosAngle)+MODEL->vertices[i].z*((oneMinusCos * yz) - xs);
      resz=MODEL->vertices[i].x*((oneMinusCos * zx) - ys)+MODEL->vertices[i].y*((oneMinusCos * yz) + xs)+MODEL->vertices[i].z*((oneMinusCos * zz) + cosAngle);
      MODEL->vertices[i].x=resx;
      MODEL->vertices[i].y=resy;
      MODEL->vertices[i].z=resz; 
      } 
 }
}

void rotateX()//rotate about object's x.here we won't change our totx,toty,totz
{thetax=alphax;
 dx=-totx;
 dy=-toty;
 dz=-totz;
 float tempx=totx;
 float tempy=toty;
 float tempz=totz;
 Translate();
 RotateX();
 dx=-dx;
 dy=-dy;
 dz=-dz;
 Translate();
 //totx=tempx;//totx,toty,totz is restored
 //toty=tempy;
 //totz=tempz;
 }
 
void rotateY()
{thetay=alphay;
 dx=-totx,dy=-toty,dz=-totz;//brought the figure about origin
 float tempx=totx;
 float tempy=toty;
 float tempz=totz;
 Translate();
 RotateY();
 dx=-dx;
 dy=-dy;
 dz=-dz;
 Translate();
 //totx=tempx;//totx,toty,totz is restored
 //toty=tempy;
 //totz=tempz;
}

void rotateZ()
{thetaz=alphaz;
 dx=-totx,dy=-toty,dz=-totz;//brought the figure about origin
 float tempx=totx;
 float tempy=toty;
 float tempz=totz;
 Translate();
 RotateZ();
 dx=-dx;
 dy=-dy;
 dz=-dz;
 Translate();
 //totx=tempx;//totx,toty,totz is restored
 //toty=tempy;
 //totz=tempz;
 }
void normalize(GLfloat x[3])
{GLdouble mag=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
 mag=sqrtf(mag);
 x[0]/=mag;
 x[1]/=mag;
 x[2]/=mag;
} 
void cross(GLfloat a[3],GLfloat b[3],GLfloat c[3])
{c[0]=a[1]*b[2]-a[2]*b[1];
 c[1]=a[2]*b[0]-a[0]*b[2];
 c[2]=a[0]*b[1]-a[1]*b[0];
} 
void myLookAt()
{
    int i;
    forward[0] = targx - eyex;
    forward[1] = targy - eyey;
    forward[2] = targz - eyez;
 
    up[0] = upx;
    up[1] = upy;
    up[2] = upz;
 
    normalize(forward);
 
    /* Side = forward x up which is u vector*/
    cross(forward, up, side);
    normalize(side);
 float m[4][4];
    
    /* Recompute up as: up = side x forward */
    cross(side, forward, up);
    normalize(up);
    CameraRotMatrix.LoadIdentity();
    CameraRotMatrix.a[0][0] = side[0];
    CameraRotMatrix.a[1][0] = side[1];
    CameraRotMatrix.a[2][0] = side[2];

    CameraRotMatrix.a[0][1] = up[0];
    CameraRotMatrix.a[1][1] = up[1];
    CameraRotMatrix.a[2][1] = up[2];

    CameraRotMatrix.a[0][2] = -forward[0];
    CameraRotMatrix.a[1][2] = -forward[1];
    CameraRotMatrix.a[2][2] = -forward[2];
    CameraRotMatrix.transpose();
    CameraTransMatrix.LoadIdentity();
    CameraTransMatrix.a[0][3]=-eyex;
    CameraTransMatrix.a[1][3]=-eyey;
    CameraTransMatrix.a[2][3]=-eyez;
    }  
void CameraRot()//Handles rotation about camera axes.Rotating world 
{
dx=-eyex;
dy=-eyey;
dz=-eyez;
Translate();
Rotate(delta1,forward[0],forward[1],forward[2]);
Rotate(delta2,up[0],up[1],up[2]);
Rotate(delta3,side[0],side[1],side[2]);
dx=eyex;
dy=eyey;
dz=eyez;
Translate();
}
void display(void)
{glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
 glColor3f(1.0,1.0,1.0);
 GLfloat model[16];
 //printf("camerarot\n");
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  myLookAt();//Computes Camera Matrix
 
 glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
 float light_position[]={eyex,eyey,eyez,0};//place the light at eye position.Camera is stationary here.
 glLightfv(GL_LIGHT0, GL_POSITION, light_position);
 glEnable(GL_LIGHTING);
 glEnable(GL_LIGHT0);
   
 glBegin(GL_TRIANGLES);
 for(int i=0;i<MODEL->numTriangles;i++){
 	TRIANGLE tri=MODEL->triangles[i];
 	Vector4 A(tri.a->x,tri.a->y,tri.a->z,1);
 	Vector4 B(tri.b->x,tri.b->y,tri.b->z,1);
 	Vector4 C(tri.c->x,tri.c->y,tri.c->z,1);
 	A=CameraRotMatrix*CameraTransMatrix*A;
 	B=CameraRotMatrix*CameraTransMatrix*B;
 	C=CameraRotMatrix*CameraTransMatrix*C;
 	glVertex3f(A.dim[0],A.dim[1],A.dim[2]);
 	glVertex3f(B.dim[0],B.dim[1],B.dim[2]);
 	glVertex3f(C.dim[0],C.dim[1],C.dim[2]);
         
 }
 glEnd();
 glFlush();
}

void reshape(int w,int h)
{ glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity ();
   gluPerspective(30,1,2,20);//PerspectiveProjection matrix.Never set znear to zero.Defines the view volume
 
}

void keyboard(unsigned char key,int x,int y)
{switch(key){
 case 'a'://move along the direction of the camera towards the target
 eyex+=.1*(targx-eyex);
 eyey+=.1*(targy-eyey);
 eyez+=.1*(targz-eyez);
 break;
 case 'A':
 eyex-=.1*(targx-eyex);
 eyey-=.1*(targy-eyey);
 eyez-=.1*(targz-eyez);
 break;
 case 'n':
 upx+=.1;
 break;
 case 'N':
 upx-=.1;
 break;
 case 'm':
 upy+=.1;
 break;
 case 'M':
 upy-=.1;
 break;
 case 'l':
 upz+=.1;
 break;
 case 'L':
 upz-=.1;
 break;
 case 'x':
 //move camera forward along x-axis
 eyex+=.1;
 break;
 case 'X':
 //move camera backward along x
 eyex-=.1;
 break;
 case 'y':
 eyey+=.1;
 //move camera forward along y
 break;
 case 'Y':
 eyey-=.1;
 //move camera backward along y
 break;
 case 'z':
 //move camera forward along -ve z
 eyez-=.1;
 break;
 case 'Z':
 eyez+=.1;
 //move camera backward along z
 break;
 case 'i':
 //rotate clockwise  about the look vector
 delta1=-1;
 delta2=0;
 delta3=0;
 CameraRot();
 break;
 case 'I':
 delta1=1;
 delta2=0;
 delta3=0;
 CameraRot();
 //Setting coordinates for target
 break;
 case 'j'://rotating about up.....To make changes
 delta1=0;
 delta2=-1;//moves camera anticlockwise to up
 delta3=0;
 CameraRot();
 break;
 case 'J':
 delta1=0;
 delta2=1;
 delta3=0;
 CameraRot();
 break;
 case 'k'://rotate anti clokwise about side
 delta1=0;
 delta2=0;
 delta3=1;
 CameraRot();
 break;
 case 'K':
 delta1=0;
 delta2=0;
 delta3=-1;
 CameraRot();
 break;
 case 'e':
 targx+=.1;
 break;
 case 'E':
 targx-=.1;
 break;
 case 'r':
 targy+=.1;
 break;
 case 'R':
 targy-=.1;
 break;
 case 't':
 targz+=.1;
 break;
 case 'T':
 targz-=.1;
 break;
 //Rotation of object about world axis
 case 'u':
 thetax=1;
 RotateX();
 break;
 case 'U':
 thetax=-1;
 RotateX();
 break;
 case 'v':
 thetay=1;
 RotateY();	
 break;
 case 'V':
 thetay=-1;
 RotateY();
 break;
 case 'w'://rotate clockwise about Z
 thetaz=1;
 RotateZ();
 break;
 case 'W':
 thetaz=-1;
 RotateZ();
 break;
 case '+'://scale
 scale=1.2;
 Scale();
 break;
 case '-'://scale
 scale=.9;
 Scale();
 break;
 }
glutPostRedisplay();//calls redisplay immediately
}
void specialkeys(int key,int x,int y)
{switch(key)
 {//translation controls
  case GLUT_KEY_LEFT:
  dx=-0.1;
  dy=0;
  dz=0;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  case GLUT_KEY_RIGHT:
  dx=.1;
  dy=0;
  dz=0;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  case GLUT_KEY_UP:
  dx=0;
  dy=.1;
  dz=0;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  case GLUT_KEY_DOWN:
  dx=0;
  dy=-.1;
  dz=0;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  case GLUT_KEY_HOME:
  dx=0;
  dy=0;
  dz=.1;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  case GLUT_KEY_END:
  dx=0;
  dy=0;
  dz=-.1;
  //totx+=dx;
 //toty+=dy;
 //totz+=dz;
  Translate();
  break;
  //Rotation about object axis controls
  case GLUT_KEY_F1:
  alphax=1;
  rotateX();
  break;
  case GLUT_KEY_F2:
  alphax=-1;
  rotateX();
  break;
  case GLUT_KEY_F3:
  alphay=1;
  rotateY();
  break;
  case GLUT_KEY_F4:
  alphay=-1;
  rotateY();
  break;
  case GLUT_KEY_F5:
  alphaz=1;
  rotateZ();
  break;
  case GLUT_KEY_F6:
  alphaz=-1;
  rotateZ();
 }
 
glutPostRedisplay();//calls redisplay immediately
 
}
void init (void) 
{   glClearColor (0.0, 0.0, 0.0, 0.0);
    GLfloat mat_specular[]={1.0,1.0,1.0,1.0};
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 0.0, 0.0, 1.0, 0.0 };
    glShadeModel(GL_SMOOTH);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
    }

int main(int argc,char** argv)
{MODEL=load_ifs_file(argv[1]);//Load Model
 
    glutInit(&argc, argv);
    
    glutInitDisplayMode (GLUT_DEPTH|GLUT_SINGLE | GLUT_RGB);
    
    glutInitWindowSize (800,600); 
    
    glutInitWindowPosition (100,100);
    
    glutCreateWindow (argv[0]);
    
    init ();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialkeys);
    glutMainLoop(); 
 return 0;
}
