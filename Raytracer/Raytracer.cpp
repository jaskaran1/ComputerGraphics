/*
Author:Jaskaran Singh
1)Primary ray(Completed)
2)Intersections(Completed)
3)Illumination(Complete)
4)Reflection(Incomplete)
*/

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
char *filename=0;
//different display modes
#define MODE_DISPLAY 1
int mode=MODE_DISPLAY;
const double PI=acos(-1.0);
//you may want to make these smaller for debugging purposes
#define WIDTH 640.0/2
#define HEIGHT 480.0/2

//the field of view of the camera
#define fov 60.0
#define MAX_RAY_DEPTH 3
#define INF 1e9

//unsigned char buffer[HEIGHT][WIDTH][3];
using namespace std;
double determinant(double a,double b,double c,double d)
{return a*d-b*c;
}

struct Vector3{
double dim[3];
Vector3(){
for(int i=0;i<3;i++)dim[i]=0;
}

Vector3(float a,float b,float c)
{dim[0]=a;
 dim[1]=b;
 dim[2]=c;}

double area(Vector3 &b)//returns magnitude of area sweeped by vector b and *this
{Vector3 result;
 result.dim[0]=determinant(dim[1],dim[2],b.dim[1],b.dim[2]);
 result.dim[1]=-determinant(dim[0],dim[2],b.dim[0],b.dim[2]);
 result.dim[2]=determinant(dim[0],dim[1],b.dim[0],b.dim[1]);
 double mag=sqrt(result.dim[0]*result.dim[0]+result.dim[1]*result.dim[1]+result.dim[2]*result.dim[2]);//absintheta
 return mag/2.;
}

Vector3 cross(Vector3 &b)
{Vector3 result;
 result.dim[0]=determinant(dim[1],dim[2],b.dim[1],b.dim[2]);
 result.dim[1]=-determinant(dim[0],dim[2],b.dim[0],b.dim[2]);
 result.dim[2]=determinant(dim[0],dim[1],b.dim[0],b.dim[1]);
 return result;
}
void normalize()
{double len=dim[0]*dim[0]+dim[1]*dim[1]+dim[2]*dim[2];
  len=sqrt(len);
  dim[0]/=len;
  dim[1]/=len;
  dim[2]/=len;
 
}
};
struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Ray
{double origin[3];
 double direction[3];
 void normalize()
 {double len=direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2];
  len=sqrt(len);
  direction[0]/=len;
  direction[1]/=len;
  direction[2]/=len;
 }
};


struct Triangle
{
  struct Vertex v[3];
  
  bool IntersectPlane(Ray ray,double &t)//Checks for intersection with the plane
  {//Calculate plane normal by cross product
   //Intersect the ray with the plane of the triangle
   Vector3 v01(v[1].position[0]-v[0].position[0] , v[1].position[1]-v[0].position[1] , v[1].position[2]-v[0].position[2]);
   Vector3 v02(v[2].position[0]-v[0].position[0] , v[2].position[1]-v[0].position[1] , v[2].position[2]-v[0].position[2]);
   Vector3 Normal;
   Normal=v01.cross(v02);
   Normal.normalize();
   double Nx=Normal.dim[0];
   double Ny=Normal.dim[1];
   double Nz=Normal.dim[2];
  // printf("%lf,%lf,%lf\n",Nx,Ny,Nz);
      //The plane passes through the point v[0] and has normal Nx,Ny,Nz
   double D=Nx*v[0].position[0]+Ny*v[0].position[1]+Nz*v[0].position[2];//Distance of plane from (0,0,0) with sign
   double Vd=Nx*ray.direction[0]+Ny*ray.direction[1]+Nz*ray.direction[2];
   double Vo=-(Nx*ray.origin[0]+Ny*ray.origin[1]+Nz*ray.origin[2]-D);
   
   if(Vd==0)//Ray parallel to plane and if Vo is 0 then ray is contained in the plane.This ray will hit triangle not on its face.Thus
   	{//if(Vo==0)
   	//	{//t=-1;//Ray lies in the plane of the triangle
   	//	 return true;}
    	 return false;//Ray is parallel to the plane of triangle	
   	}
   t=Vo/Vd;
   if(t>1e-11)
   	return true;
   return false;		
  }
  
  bool IntersectTriangle(Ray ray,double xi,double yi,double zi,double &Nx,double &Ny,double &Nz)
  {//Find dominant component of normal of plane.Project triangle onto this plane
  //Normal has to be calculated using interpolation
   double alpha,beta,DET;
   Vertex v0=v[0];
   Vertex v1=v[1];
   Vertex v2=v[2];
   double x0=v0.position[0];
   double x1=v1.position[0];
   double x2=v2.position[0];
   double y0=v0.position[1];
   double y1=v1.position[1];
   double y2=v2.position[1];
   double z0=v0.position[2];
   double z1=v1.position[2];
   double z2=v2.position[2];
   Vector3 v01,v02,vi0,vi1,vi2;	 
 	 v01.dim[0]=v1.position[0]-v0.position[0];
         v01.dim[1]=v1.position[1]-v0.position[1];
         v01.dim[2]=v1.position[2]-v0.position[2];
 
 	 v02.dim[0]=v2.position[0]-v0.position[0];
         v02.dim[1]=v2.position[1]-v0.position[1];
         v02.dim[2]=v2.position[2]-v0.position[2];
 
 	 vi0.dim[0]=v0.position[0]-xi;//vector from intersection point to V0
 	 vi0.dim[1]=v0.position[1]-yi;
 	 vi0.dim[2]=v0.position[2]-zi;
 
 	 vi1.dim[0]=v1.position[0]-xi;//vector from intersection point to V1
 	 vi1.dim[1]=v1.position[1]-yi;
 	 vi1.dim[2]=v1.position[2]-zi;
 
 	 vi2.dim[0]=v2.position[0]-xi;//vector from intersection point to V2
 	 vi2.dim[1]=v2.position[1]-yi;
 	 vi2.dim[2]=v2.position[2]-zi;	 	 
 	 
 double AREA=v01.area(v02);//cross product v01 X v02
 double area0=vi1.area(vi2);//area of triangle opposite to vertex v0
 double area1=vi0.area(vi2);//area of triangle opposite to vertex v1
 double area2=vi0.area(vi1);//area of triangle opposite to vertex v2
 Nx=area0*v[0].normal[0]+area1*v[1].normal[0]+area2*v[2].normal[0];
 Ny=area0*v[0].normal[1]+area1*v[1].normal[1]+area2*v[2].normal[1];
 Nz=area0*v[0].normal[2]+area1*v[1].normal[2]+area2*v[2].normal[2];
 
   
   if(abs(Nx)>=abs(Ny) && abs(Nx)>=abs(Nz))
   	{//magnitude of x component of normal is largest
   	 //we remove equation of x
   	 DET=determinant(y1-y0,y2-y0,z1-z0,z2-z0);
   	 alpha=determinant(yi-y0,y2-y0,zi-z0,z2-z0)/DET;
   	 beta=determinant(y1-y0,yi-y0,z1-z0,zi-z0)/DET;
   	}
   else if(abs(Ny)>=abs(Nx) && abs(Ny)>=abs(Nz))
   	{//y component largest
   	 //remove y equation
   	 DET=determinant(x1-x0,x2-x0,z1-z0,z2-z0);
   	 alpha=determinant(xi-x0,x2-x0,zi-z0,z2-z0)/DET;
   	 beta=determinant(x1-x0,xi-x0,z1-z0,zi-z0)/DET;
   	}
   else
   	{//z component largest
   	 DET=determinant(x1-x0,x2-x0,y1-y0,y2-y0);
   	 alpha=determinant(xi-x0,x2-x0,yi-y0,y2-y0)/DET;
   	 beta=determinant(x1-x0,xi-x0,y1-y0,yi-y0)/DET;
        }
   return alpha>=0 && beta>=0 && alpha+beta<=1;} 	
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
  
  bool IntersectSphere(Ray ray,double &tnear)
  {
  
  double B=2*(ray.direction[0]*(ray.origin[0]-position[0])+ray.direction[1]*(ray.origin[1]-position[1])+ray.direction[2]*(ray.origin[2]-position[2]));
  
  double C=(ray.origin[0]-position[0])*(ray.origin[0]-position[0])+(ray.origin[1]-position[1])*(ray.origin[1]-position[1])+(ray.origin[2]-position[2])*(ray.origin[2]-position[2])-radius*radius;
 //printf("A=%lfB=%lfC=%lf\n",A,B,C);
  if(B*B-4*C<=0)
  	return false;//didn't intersect 
  double D=sqrt(B*B-4*C);
  double t1=(-B-D)*0.5;
  double t2=(-B+D)*0.5;
  if(t1>=0&&t2>0)
  	{tnear=t1;
  	 //printf("%lf\n",tnear);
  	 return true;
  	}
  if(t1<0&&t2>1e-11)
  	{tnear=t2;
  	//printf("%lf\n",tnear);
  	return true;
  	}
  return false;		
  		
  }
};

struct Light
{
  double position[3];
  double color[3];
};

struct pixel
{double r,g,b;
};

vector<Triangle>triangles;//vector of triangles
vector<Sphere>spheres;//array of spheres
vector<Light>lights;//array of lights
double ambient_light[3];//constant through  a scene

Sphere SphereIntersection(Ray& curr,double& xi,double& yi,double& zi,double& Nx,double& Ny,double& Nz,int& flag,double& nearest)
{//printf("Curr ray=[%lf,%lf,%lf]\n",curr.direction[0],curr.direction[1],curr.direction[2]);
 Sphere s;
 //s.radius=-1;//initialize value
 for(vector<Sphere>::iterator it=spheres.begin();it!=spheres.end();it++)
 {double tnear;
  if(it->IntersectSphere(curr,tnear))//Intersection is confirmed only when t0 is positive.else the ray originates inside or it intersects backward
    	{//printf("tnear=%lf\n",tnear);
    	 if(tnear<nearest)	
   		{flag=1;
   		 
   		 nearest=tnear;
    		 
    		 xi=curr.origin[0]+tnear*curr.direction[0];//update points
 	 	 yi=curr.origin[1]+tnear*curr.direction[1];
 	 	 zi=curr.origin[2]+tnear*curr.direction[2];
 	 	 //printf("zi=%lf\n",zi);
    		 Nx=(xi-it->position[0])/it->radius;//update normals
    		 Ny=(yi-it->position[1])/it->radius;
    		 Nz=(zi-it->position[2])/it->radius;
    		 //printf("Nz=%lf\n",Nz);
    		 double len=sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
    		 Nx/=len;
    		 Ny/=len;
    		 Nz/=len;
    		 //printf("%lf,%lf,%lf\n",xi,yi,zi);
    		 s=*it;
    		}  
  	}
 
 }
 return s;//Returns the intersected sphere
}

Triangle TriangleIntersection(Ray& curr,double& xi,double& yi,double& zi,double& Nx,double& Ny,double& Nz,int& flag,double& nearest)
{
 Triangle t;
 for(vector<Triangle>::iterator it=triangles.begin();it!=triangles.end();it++)
  {double tnear;
   if(it->IntersectPlane(curr,tnear))//Plane intersection working correctly
  	{double ai,bi,ci;//They're the point of intersection on the plane whereas xi,yi,zi are the points of intersection on the triangle
  	 ai=curr.origin[0]+tnear*curr.direction[0];
 	 bi=curr.origin[1]+tnear*curr.direction[1];
 	 ci=curr.origin[2]+tnear*curr.direction[2];
    	 double Na,Nb,Nc;	 
  	 if(it->IntersectTriangle(curr,ai,bi,ci,Na,Nb,Nc))//Only if this is true then intersection is done
  		{double dist=sqrt((ai-curr.origin[0])*(ai-curr.origin[0])+(bi-curr.origin[1])*(bi-curr.origin[1])+(ci-curr.origin[2])*(ci-curr.origin[2]));//calculate distance
  	 	 if(dist<nearest)//then we update xi,yi,zi and nearest
  	 		{xi=ai;
  	 		 yi=bi;
  	 		 zi=ci;
  	 		 Nx=Na;
  	 		 Ny=Nb;
  	 		 Nz=Nc;
  	 		 flag=2; 
  	 		 nearest=dist;
  	 	 	 t=*it;
  	 		}
  		}
  	}
  }
return t;
}
Ray ComputeReflectedRay(Ray Lm,double Nx,double Ny,double Nz,double xi,double yi,double zi)
{	Lm.normalize();
        double magnorm=sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
    	Nx/=magnorm;
    	Ny/=magnorm;
    	Nz/=magnorm;
    	//printf("Nx=%lf,Ny=%lf,Nz=%lf",Nx,Ny,Nz);
	double costheta=Nx*Lm.direction[0]+Ny*Lm.direction[1]+Nz*Lm.direction[2];
	double k=2*costheta;
	Ray reflect;
	reflect.origin[0]=xi;
	reflect.origin[1]=yi;
	reflect.origin[2]=zi;
	reflect.direction[0]=k*Nx-Lm.direction[0];
	reflect.direction[1]=k*Ny-Lm.direction[1];
	reflect.direction[2]=k*Nz-Lm.direction[2];
        reflect.normalize();
        return reflect;}
    
	
void Interpolate(double xi,double yi,double zi,Triangle T,double color_diffuse[],double color_specular[],double& shininess)
{Vertex v0=T.v[0];
 Vertex v1=T.v[1];
 Vertex v2=T.v[2];
 Vector3 v01,v02,vi0,vi1,vi2;	 
 	 v01.dim[0]=v1.position[0]-v0.position[0];
         v01.dim[1]=v1.position[1]-v0.position[1];
         v01.dim[2]=v1.position[2]-v0.position[2];
 
 	 v02.dim[0]=v2.position[0]-v0.position[0];
         v02.dim[1]=v2.position[1]-v0.position[1];
         v02.dim[2]=v2.position[2]-v0.position[2];
 
 	 vi0.dim[0]=v0.position[0]-xi;//vector from intersection point to V0
 	 vi0.dim[1]=v0.position[1]-yi;
 	 vi0.dim[2]=v0.position[2]-zi;
 
 	 vi1.dim[0]=v1.position[0]-xi;//vector from intersection point to V1
 	 vi1.dim[1]=v1.position[1]-yi;
 	 vi1.dim[2]=v1.position[2]-zi;
 
 	 vi2.dim[0]=v2.position[0]-xi;//vector from intersection point to V2
 	 vi2.dim[1]=v2.position[1]-yi;
 	 vi2.dim[2]=v2.position[2]-zi;	 	 
 	 
 double AREA=v01.area(v02);//cross product v01 X v02
 double area0=vi1.area(vi2);//area of triangle opposite to vertex v0
 double area1=vi0.area(vi2);//area of triangle opposite to vertex v1
 double area2=vi0.area(vi1);//area of triangle opposite to vertex v2
 
 color_diffuse[0]=(area0*v0.color_diffuse[0]+area1*v1.color_diffuse[0]+area2*v2.color_diffuse[0])/AREA;//
 color_diffuse[1]=(area0*v0.color_diffuse[1]+area1*v1.color_diffuse[1]+area2*v2.color_diffuse[1])/AREA;//
 color_diffuse[2]=(area0*v0.color_diffuse[2]+area1*v1.color_diffuse[2]+area2*v2.color_diffuse[2])/AREA;
 
 color_specular[0]=(area0*v0.color_specular[0]+area1*v1.color_specular[0]+area2*v2.color_specular[0])/AREA;
 color_specular[1]=(area0*v0.color_specular[1]+area1*v1.color_specular[1]+area2*v2.color_specular[1])/AREA;
 color_specular[2]=(area0*v0.color_specular[2]+area1*v1.color_specular[2]+area2*v2.color_specular[2])/AREA;
 
 shininess=(area0*v0.shininess+area1*v1.shininess+area2*v2.shininess)/AREA;
 
}
double CalcSpecularComponent(Ray reflected,Ray view,double shine)
{
view.direction[0]=-view.direction[0];
view.direction[1]=-view.direction[1];
view.direction[2]=-view.direction[2];
double ans=reflected.direction[0]*(view.direction[0]) + reflected.direction[1]*(view.direction[1]) + reflected.direction[2]*(view.direction[2]);
if(ans<=0)
	{//printf("Spec is negative\n");
	 return 0;}
ans=pow(ans,shine);
	{
	return ans;}
}

double CalcDiffuseComponent(Ray shadowray,double Nx,double Ny,double Nz)
{double ans=Nx*(shadowray.direction[0])+Ny*(shadowray.direction[1])+Nz*(shadowray.direction[2]);
 ans/=(sqrt((Nx*Nx+Ny*Ny+Nz*Nz)*(shadowray.direction[0]*shadowray.direction[0]+shadowray.direction[1]*shadowray.direction[1]+shadowray.direction[2]*shadowray.direction[2])));
 if(ans>0)
 return ans;
 //else 
 return 0;
}


pixel trace(Ray curr,int depth)//recursive function
{
 pixel pix={0.,0.,0.};//pix will store the color of the pixel.Here initialised to black
 
 double nearest=INF;//The distance of the nearest point from the origin of the curr ray
 //printf("nearest=%lf\n",nearest);
 
 int flag=0;//this is flag for primary ray intersection
 int flag2=0;//this is flag for shadow ray intersection
 double xi,yi,zi;//nearest intersection point
 
 double Nx,Ny,Nz;//outward normal at the intersection point.For spheres outward normal is by default
 //double Tx,Ty,Tz;//outward normal for triangle
 Sphere S=SphereIntersection(curr,xi,yi,zi,Nx,Ny,Nz,flag,nearest);//WORKING CORRECTLY
 
 //after leaving for loop either flag=0 which means no intersection,else yes.
 //Check for intersection with triangles
 Triangle T=TriangleIntersection(curr,xi,yi,zi,Nx,Ny,Nz,flag,nearest);//WORKING CORRECTLY
 if(flag==0)
 	return pix;//Black.Since no object was intersected
 	
pix.r=ambient_light[0];
pix.g=ambient_light[1];
pix.b=ambient_light[2];  	
	//return pix;//Just made a change  

 //I = (ka * Ia) + Ii * (kd * (L . N) + ks * (R . V)n)
 Ray shadowray;//this ray would go to different light sources from the intersection point
 shadowray.origin[0]=xi;
 shadowray.origin[1]=yi;
 shadowray.origin[2]=zi;
 //printf("Nz=%lf\n",Nz);
 //printf("Reflected ray:(%lf,%lf,%lf)\n",Reflect.direction[0],Reflect.direction[1],Reflect.direction[2]);
 vector<Light>::iterator it;
 for(it=lights.begin();it!=lights.end();it++)
 	{
 	 shadowray.direction[0]=it->position[0]-xi;
 	 shadowray.direction[1]=it->position[1]-yi;
 	 shadowray.direction[2]=it->position[2]-zi;
 	 double distance=sqrt(shadowray.direction[0]*shadowray.direction[0]+shadowray.direction[1]*shadowray.direction[1]+shadowray.direction[2]*shadowray.direction[2]);
 	 shadowray.normalize();
 	 Ray Reflect=ComputeReflectedRay(shadowray,Nx,Ny,Nz,xi,yi,zi);
         flag2=0;
 	 //nearest=INF;
 	 double xsi=INF,ysi=INF,zsi=INF;//intersection point of shadow ray
 	 double Nxs,Nys,Nzs;//normal for shadow ray
 	 double nearestshadow=distance;//initializing the distance to be distance between point of intersection and light source origin
 	 //double nearestshadow=INF;
 	 double color_diffuse[3],color_specular[3],shininess;//interpolated color_diffuse,color_specular and shininess for triangle intersection
 	 SphereIntersection(shadowray,xsi,ysi,zsi,Nxs,Nys,Nzs,flag2,nearestshadow);
 	 Triangle u=TriangleIntersection(shadowray,xsi,ysi,zsi,Nxs,Nys,Nzs,flag2,nearestshadow);
 	 if(flag2==0)//reaches light source without any intersection	    
 	 	{
 	 	 if(flag==1)//the primary light ray intersected the sphere.So colour properties of that sphere will be used
 	 	 	{shininess=S.shininess;
 	 	 	 double spec=CalcSpecularComponent(Reflect,curr,shininess);
 	 	 	 double diff=CalcDiffuseComponent(shadowray,Nx,Ny,Nz);
 	 	 	 pix.r+=(it->color[0])*(S.color_specular[0]*spec+S.color_diffuse[0]*diff);
 	 	 	 pix.g+=(it->color[1])*(S.color_specular[1]*spec+S.color_diffuse[1]*diff);
 	 	 	 pix.b+=(it->color[2])*(S.color_specular[2]*spec+S.color_diffuse[2]*diff);
 	 	 	}
 	 	 else//it intersected a triangle
 	 	 	{Interpolate(xi,yi,zi,T,color_diffuse,color_specular,shininess);
 	 	 	 double spec=CalcSpecularComponent(Reflect,curr,shininess);
 	 	 	 double diff=CalcDiffuseComponent(shadowray,Nx,Ny,Nz);
 	 	 	 pix.r+=(it->color[0])*(color_diffuse[0]*diff+color_specular[0]*spec);
 	 	 	 pix.g+=(it->color[1])*(color_diffuse[1]*diff+color_specular[1]*spec);
 	 	 	 pix.b+=(it->color[2])*(color_diffuse[2]*diff+color_specular[2]*spec);
 	 		}
 	 	  //Code for recursive reflection
 	 	  /*pixel pixref;
 	 	  if(shininess>10&&depth<MAX_RAY_DEPTH)//extra calculation for recursive reflective.Code for illumination by reflected ray
 			pixref=trace(Reflect,depth+1);
 		  pix.r+=.5*pixref.r;
 	          pix.g+=.5*pixref.g;
 		  pix.b+=.5*pixref.b;*/
	
 	 	}
 	  }
return pix; 	 
 

}
void plot_pixel_display(double x,double y,double z,double r,double g,double b)
{
  glColor3d(r,g,b);
  glVertex3d(x,y,z);
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  double x,y;
  double tanhalffov=tan((PI*fov)/360.0000);
  double camz=HEIGHT/(2*tanhalffov);//Distance of projection plane from camera.Camera is behind the xy plane
  for(x=-WIDTH/2;x<=WIDTH/2;x++)//Sending primary rays
  {
  glBegin(GL_POINTS);
    
    for(y=-HEIGHT/2;y<=HEIGHT/2;y++)
    {//printf("%lf\n",y);
    //compute primary ray direction for each pixel.Camera is at 0,0,0
    Ray primRay;
    primRay.origin[0]=primRay.origin[1]=0;
    primRay.origin[2]=0;
    primRay.direction[0]=x;
    primRay.direction[1]=y;
    primRay.direction[2]=-camz;
    primRay.normalize();//Normalized primary ray
    int depth=0;
    pixel pix=trace(primRay,depth);//trace returns the color of the pixel
    plot_pixel_display(x,y,0,pix.r,pix.g,pix.b);
    
    }
    glEnd();
    glFlush();
  
    }
  printf("Done!\n"); fflush(stdout);
}


void parse_check(const char *expected,const char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file,const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)//strcasecmp ignores the case
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  triangles.push_back(t);
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);
	  spheres.push_back(s);
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);
	  lights.push_back(l);
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{static int once=0;
  if(!once)
  {
      draw_scene();
       }
  once=1;
//draw_scene();
}

void init()
{ glViewport(-WIDTH/2,-HEIGHT/2,WIDTH,HEIGHT);
  glMatrixMode(GL_PROJECTION);
  //gluPerspective(180,double(WIDTH/HEIGHT),1,5);
  glOrtho(-WIDTH/2,WIDTH/2,-HEIGHT/2,HEIGHT/2,0,1000);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}


int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  init();
  glutMainLoop();
}
