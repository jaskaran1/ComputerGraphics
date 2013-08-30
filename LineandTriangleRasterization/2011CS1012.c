//Bresenham Algorithm to draw lines
//Triangle Rasterization
//For implementation details refer to README
#include <stdio.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include<math.h>
#include<stdlib.h>

int numLines,numTriangles;

typedef struct point{
float x,y;//for coordinates
float r,g,b,a;//for colors
}point;

typedef struct line{
struct point start,finish;
}line;
line *l;

typedef struct triangle{
struct point v[3];//3 vertices	
}tri;
tri *t;

int cmp(const void *e1,const void *e2)//compare function for the sort by x coordinate
{point u=*((point*)e1);
 point v=*((point*)e2);
 if(u.x<v.x)
 	return -1;
 if(u.x>v.x)
 	return 1;
 return 0;
}

int create_line(line l,point *list,int ind){//list is the list of all points interpolated on the line l1,ind is the starting index of list
 //if list is NULL.This means that the list of the points on the interpolated line isn't needed.
 //The function returns 1+the index of the final point on the line created.
 int x0=l.start.x;
 int y0=l.start.y;
 float r0=l.start.r;
 float g0=l.start.g;
 float b0=l.start.b;
 float a0=l.start.a;
 
 int x1=l.finish.x;
 int y1=l.finish.y;
 float r1=l.finish.r;
 float g1=l.finish.g;
 float b1=l.finish.b;
 float a1=l.finish.a;
 
 int flag=0;//this is set 1 when m>1.set -1 when m<-1.This is used to writepixel
 
 if(x0>x1){    //since we're dealing with only 4 octants where x>0.we'll draw the line from left to right.this is done to get dx positive
      int t=x0;//swap x
      x0=x1;
      x1=t;
  
      t=y0;//swap y
      y0=y1;
      y1=t;
  
      float tc=r0;//swap red
      r0=r1;
      r1=tc;
  
      tc=g0;//swap green
      g0=g1;
      g1=tc;
  
      tc=b0;//swap blue
      b0=b1;
      b1=tc;
  
      tc=a0;//swap a
      a0=a1;
      a1=tc;}
      
 int dx=x1-x0;//dx is positive
 int dy=y1-y0;
 
 if(dy-dx>0){//m>1.(x,y)->(y,x).This is done so as to project the line in the 1st octant where 0<=m<=1 by taking reflection in y=x
     int t=x0;
     x0=y0;
     y0=t;
     t=x1;
     x1=y1;
     y1=t;
     flag=1;}
     
 if(dy+dx<0){//m<-1(x,y)->(-y,-x).This is done so as to project the line in the octant where 0>m>=-1  by taking reflection in y=-x
     int t=-x0;
     x0=-y0;
     y0=t;
     t=-x1;
     x1=-y1;
     y1=t;
     flag=-1;}
     
 dy=y1-y0;//projected dy
 dx=x1-x0;//projected dx.positive 
          //such that -1<dy/dx<1
 int d=dy*2-dx;
 int incrE=dy*2;
 int incrNE=(dy-dx)*2;
 
 if(dy<0){//-1<m<0
    d=dy*2+dx;
    incrE=(dy+dx)*2;
    incrNE=2*dy;}
 
 int x=x0;
 int y=y0;
 float r=r0;
 float g=g0;
 float b=b0;
 float a=a0;
 glColor4f(r,g,b,a);
   if(flag==1){//m>1
      glVertex2d(y,x);//bringing it back in the m>1 octant from m<1 octant
      if(list!=NULL)
      {list[ind].x=y;
       list[ind].y=x;}
      }
      
   else if(flag==-1){//m<-1
      glVertex2d(-y,-x);//bringing it back from the m>-1 to m<-1 octant
      if(list!=NULL)
      {list[ind].x=-y;
       list[ind].y=-x;}
      }
 
   else{
       glVertex2d(x,y);
       if(list!=NULL)
       {list[ind].x=x;
        list[ind].y=y;}
       }
       if(list!=NULL)
       {list[ind].r=r;
        list[ind].g=g;
        list[ind].b=b;
        list[ind].a=a;
        ind++;}
 x++;//since x0 is already plotted we start from next vertex
 int c=1;//used for color interpolation of the line
 
 for(;x<=x1;x++,c++){//Color interpolation
   float u=(float)c/(x1-x0);
   r=r0*(1.0-u)+r1*u;
   g=g0*(1.0-u)+g1*u;
   b=b0*(1.0-u)+b1*u;
   a=a0*(1.0-u)+a1*u;
   
   if(d<=0){//when midpoint is above the line or on it
       d+=incrE;
       if(dy<0)
           y--;}
   
   else{//when midpoint is below the line
       d+=incrNE;
       if(dy>0)
           y++;}
  glColor4f(r,g,b,a);
  int xp=x,yp=y;//(xp,yp) point which is actually plotted 
  
  if(flag==1){
       xp=y;
       yp=x;}
  
  if(flag==-1){
   xp=-y;
   yp=-x;}
  
  glColor4f(r,g,b,a); 
  glVertex2d(xp,yp);
  if(list!=NULL){
  list[ind].x=xp;
  list[ind].y=yp;
  list[ind].r=r;
  list[ind].g=g;
  list[ind].b=b;
  list[ind].a=a;
  ind++;}
      
  
 }
 return ind;
 }

void create_triangle(tri t){//creates triangle with rasterization
 
  qsort(t.v,sizeof(t.v)/sizeof(*t.v),sizeof(*t.v),cmp);//sorts the points with x coordinate as key
  
  line l1={t.v[0],t.v[2]};//line joining from min to max
  
  int size=(t.v[2].x-t.v[0].x)+1;//size=largest-smallest+1 points
  
  point arr1[size];//will include all the points on line l1
  point arr2[size];//will include points on line l2 and l3
  
  create_line(l1,arr1,0);
  int i;
  
  line l2={t.v[0],t.v[1]};//line joining from min to medium
  int curr=create_line(l2,arr2,0);//curr is the index which will be starting for points on l3
  
  line l3={t.v[1],t.v[2]};//line joining from medium to max
  create_line(l3,arr2,curr);//creates the boundary of the triangle
  
  for(i=0;i<size;i++){
     line l={arr1[i],arr2[i]};
     create_line(l,NULL,-1);}
} 

void display(void)
{   glClear (GL_COLOR_BUFFER_BIT);
    
    glBegin(GL_POINTS);
    
    int i;
    
    for(i=0;i<numLines;i++)
     create_line(l[i],NULL,-1); 
    
    for(i=0;i<numTriangles;i++)
     create_triangle(t[i]);  
    
    glEnd();
    
    glFlush();//execute the commands in the buffer.don't wait for more drawing commands
}

void init (void) 
{   glClearColor (0.0, 0.0, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1000.0, 1000.0, -1000.0, 1000.0,-1.0,1.0);
}

int main(int argc, char** argv)
{   FILE *f1=fopen(argv[1],"r");

    if(!f1){
     printf("Error opening file\n");
     return 1;}
    
    glutInit(&argc, argv);
    
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    
    glutInitWindowSize (500, 500); 
    
    glutInitWindowPosition (100,100);
    
    glutCreateWindow (argv[0]);
    
    init ();
    
    glutDisplayFunc(display); 
    
    fscanf(f1,"%d",&numLines);
    
    l=(line*)malloc(numLines*sizeof(line));
    
    int i;
    
    for(i=0;i<numLines;i++)
    {fscanf(f1,"%f%f%f%f%f%f",&l[i].start.x,&l[i].start.y,&l[i].start.r,&l[i].start.g,&l[i].start.b,&l[i].start.a);
     fscanf(f1,"%f%f%f%f%f%f",&l[i].finish.x,&l[i].finish.y,&l[i].finish.r,&l[i].finish.g,&l[i].finish.b,&l[i].finish.a);}
     
    fscanf(f1,"%d",&numTriangles);
    
    t=(tri*)malloc(numTriangles*sizeof(tri));
    
    for (i = 0; i <numTriangles; ++i){
     fscanf(f1,"%f%f%f%f%f%f",&t[i].v[0].x,&t[i].v[0].y,&t[i].v[0].r,&t[i].v[0].g,&t[i].v[0].b,&t[i].v[0].a);
     fscanf(f1,"%f%f%f%f%f%f",&t[i].v[1].x,&t[i].v[1].y,&t[i].v[1].r,&t[i].v[1].g,&t[i].v[1].b,&t[i].v[1].a);
     fscanf(f1,"%f%f%f%f%f%f",&t[i].v[2].x,&t[i].v[2].y,&t[i].v[2].r,&t[i].v[2].g,&t[i].v[2].b,&t[i].v[2].a);}
    
    glutMainLoop();
    
    return 0;   /* ISO C requires main to return int. */
}
