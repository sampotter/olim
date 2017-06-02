/* Ordered Line Integral Method */
/* solves Eq. |grad g| = f(x,y); Line Integral: int_0^L|g(\phi)|ds */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
#define BETA 0.0
#define NX 21
#define NY 21
#define XMIN -1.0
#define XMAX 1.0
#define YMIN -1.0
#define YMAX 1.0


struct mycurve {
	double x;
	double y;
	struct mycurve *next;
	struct mycurve *prev;
};


struct myvector {
  double x;  
  double y;
};


struct mymatrix {
  double a11;
  double a12;
  double a21;
  double a22;
};

struct mysol {
  double g;
  char c;
};  



int main(void);
double slowness(double x,double y); 
double exact_solution(double x,double y);
double angle(double x,double y);
double length(double x,double y);
void param(void);
void fast_marching_method(void);
struct mysol triangle_update(long ind,long ind0,long ind1);
double one_pt_update(long ind,long ind0);			   
void addtree(long ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(long ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
struct mymatrix matrix_inverse(struct mymatrix matr);
struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b);
struct mymatrix transpose(struct mymatrix matr);
struct myvector matr_vec(struct mymatrix matr,struct myvector vec);
double dot_product(struct myvector a,struct myvector b); 
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
struct myvector getpoint(long ind);
double length_vec(struct myvector x);

struct mysol nonlin_solver(double u0,double u1,
			struct myvector x0,struct myvector x1,double f,struct myvector x);
double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,double b,
						struct myvector x1mx0,struct myvector x);
double geometric_action_line(double x0,double y0,double x1,double y1);
/***************************************/

const long nx1 = NX - 1, ny1 = NY - 1, nxy = NX*NY, nx2 = 2*NX, ny2 = 2*NY;
long count; /* # of considered points */
double h,hx,hy;
int ms[NX*NY]; /* 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double Uexact[NX*NY], g[NX*NY]; /* the exact solution and the computed solution */
double slo[NX*NY]; /* slowness computed on the twice refined mesh */
long pos[NX*NY]; /* pos(index of mesh pt) = position in binary tree */
long tree[NX*NY]; /* tree(position in the tree) = index of mesh pt */
const long neii[4]={1, NX, -1, -NX }; /* neighbor's indices */
double xa, ya; /* (xa,ya) is the initial point */

const char *f_output_name = "gfun.txt";

long INDUPDATE;

/**************************************/
double slowness(double x,double y) {
  double aux,aux0,aux1,vx,vy,v;
  
  aux = 0.5*(x + y);
  aux0 = cos(aux);
  aux1 = 0.5*sin(aux);
  vx = (x + aux0)*(1.0 - aux1);
  vy = y - aux0*aux1;
  v = sqrt(vx*vx + vy*vy);
  
  return 1.0;
}

/*************************************/
double exact_solution(double x,double y) {
  double aux,aux0;
  
  aux = 0.5*(x + y);
  aux0 = x + cos(aux);
  
  return length(x,y); //0.5*(aux0*aux0 + y*y);
}



/*************************************/

void param() {
  long i,j,ind;  
  struct myvector xy;
  
  xa = 0.0; //-0.900367222589747;
  ya = 0.0;

  count = 0; /* set the number of Considered points to zero */	
  printf("in param()\n");
  hx = (XMAX - XMIN)/(NX-1);
  hy = (YMAX - YMIN)/(NY-1);
  h = sqrt(hx*hx + hy*hy);
  for( i=0; i<NX; i++ ) {
    for( j=0; j<NY; j++ ) {
      ind = i + NX*j;
	  ms[ind] = 0;
	  g[ind] = INFTY;
	  xy = getpoint(ind);
	  Uexact[ind] = exact_solution(xy.x,xy.y);
    }
  }
  for( i=0; i<NX; i++ ) {
    for( j=0; j<NY; j++ ) {
      ind = i + NX*j;
	  slo[ind] = slowness(XMIN + hx*i,YMIN + hy*j);
	}
  }
}

/************************************/

void ipoint() {
  long i,j,ind,ind0,m,n;
  double gtemp;
  const long isur[4]={0, 1, NX+1, NX};
  struct myvector x,l,b;

  i = floor((xa - XMIN)/hx);
  j = floor((ya - YMIN)/hy);
  ind0 = i+j*NX;
  for( m=0; m<4; m++ ) {
    ind = ind0+isur[m];
	x = getpoint(ind);
    gtemp = geometric_action_line(xa,ya,x.x,x.y);
	g[ind] = gtemp;
	ms[ind] = 1;
	addtree(ind);
  }
}



/**********************************************/
/*** ordered line integral method ***/

void fast_marching_method(void) {
  long i0,j0,k,m,n,ind,ind0,ind1;
  int mycount = 0; /* counter for accepted points */
  double gold,gtemp;
  long Naf, AFneib[4], Nc, Cneib[4]; /* accepted front neighbors and new considered of the newly accepted point */
  struct mysol sol; 
  struct myvector x,x0,x1; 
  double b;
  
  int NAC = 0;
  
  printf("in OLIM()\n");

  while( count > 0 ) {
  	NAC++;
    ind0 = tree[1];
    x0 = getpoint(ind0);
    j0 = ind0/NX;
    i0 = ind0%NX;
    /* x and y of the newly accepted point */
    ms[ind0] = 2;
    deltree();
	mycount++;
    if( i0==1 || i0==nx1-1 || j0==1 || j0== ny1-1 || g[ind0] >= INFTY-1) {
      printf("The boundary is reached:\n");
	  printf("%d\t(%ld\t%ld) is accepted, g=%.4e\n",mycount,i0,j0,g[ind0]);
	  break; /* quit if we reach the boundary of the computational domain */
	}

//   	printf("%d\t(%ld\t%ld) is accepted, g=%.4e\n",mycount,i0,j0,g[ind0]); 

    /* Inspect the neighbors of the new Accepted point */
    Naf = 0;
    Nc = 0;
    for( k=0; k<4; k++ ) {
	  ind1 = ind0 + neii[k];
	  if( ms[ind1] <= 1 ) { /*  Considered neighbors */
		Cneib[Nc] = ind1;
		Nc++;
	  }  
	}
//  	printf("New Accepted: (%li,%li)\n",i0,j0);
// 	printf("Nc = %li\n",Nc);
	/* update all neighbors of the new Accepted point that are not Accepted */
	for( k = 0; k < Nc; k++ )  {
		ind = Cneib[k];
		gold = g[ind];
		g[ind]  = min(one_pt_update(ind,ind0),g[ind]);
		x = getpoint(ind);
		b = slo[ind];
		for( m = 0; m < 4; m++ ) {
			ind1 = ind + neii[m];
			if( ms[ind1] == 2 && ind1 != ind0 && labs(i0 - ind1%NX) <= 1 && labs(j0 - ind1/NX) <= 1 ) {			
				x1 = getpoint(ind1);
				sol = nonlin_solver(g[ind0],g[ind1],x0,x1,b,x);
// 				printf("back from nonlin_solver: %.4e\t%c\n",sol.g,sol.c);
				if( sol.c == 'y' ) {
					gtemp = sol.g;
					g[ind] = min(g[ind],gtemp);
// 					printf("(%li,%li),(%li,%li): gtemp = %.4e\t",i0,j0,ind1%NX,ind1/NX,gtemp);	
				}
			}	
		}
		if( ms[ind] == 1 && gold > g[ind] ) {
			updatetree(ind);
		} 
		else if( ms[ind] == 0 ) {
			ms[ind] = 1;
			addtree(ind);
		}
// 		printf("k = %li: indupdate= %li: (%li,%li): gold = %.4e, gnew = %.4e\n",
// 			k,ind,ind%NX,ind/NX,gold,g[ind]);		  
	}		
  } /* end while ( count > 0 ) */
}



  
/*********************************************/

double one_pt_update(long ind,long ind0) {
  struct myvector x,x0;
  double gtemp,l;
//   long idiff = labs(ind - ind0),i,j;
  
  x = getpoint(ind);
  x0 = getpoint(ind0);
  l = length_vec(vec_difference(x,x0));
  gtemp = g[ind0] + l*slo[ind];
//   printf("one_pt_update((%li,%li),(%li,%li):  gtemp = %.4e\n",
//   		ind%NX,ind/NX,ind0%NX,ind0/NX,gtemp);
  
//  fprintf(fup,"%li\t%.6e\t%i\t%li\t%i\t%i\t%.4e\n",ind,gtemp,1,ind0,-1,-1,length(l.x,l.y));

  return gtemp;
}
/*-------------*/

double geometric_action_line(double x0,double y0,double x1,double y1) {
  struct myvector l;
  double f;
  
  l.x = x1 - x0;
  l.y = y1 - y0;
  f = slowness(x1,y1);

  return f*length_vec(l);
}


/*-------------*/  
  
struct myvector getpoint(long ind) {
  struct myvector l;
  
  l.x = hx*(ind%NX) + XMIN;
  l.y = hy*(ind/NX) + YMIN;
  return l;
}

/***** N o n l i n e a r   1D    s o l v e r *****/


struct mysol nonlin_solver(double u0,double u1,
			struct myvector x0,struct myvector x1,double f,struct myvector x) {

	double A,B,C,du = u0 - u1,du2,a,b,c,f2 = f*f,s,s1,discr;
	struct myvector xmx1,x1mx0,xs;
	struct mysol sol;
	
	sol.c = 'y';
	sol.g = INFTY;
	xmx1 = vec_difference(x,x1);
	x1mx0 = vec_difference(x1,x0);
	A = x1mx0.x*x1mx0.x + x1mx0.y*x1mx0.y;
	B = dot_product(xmx1,x1mx0);
	C = xmx1.x*xmx1.x + xmx1.y*xmx1.y;
	du2 = du*du;
	a = A*(A*f2 - du2);
	b = B*(f2*A - du2);
	c = f2*B*B - du2*C;
	discr = b*b - a*c;
	if( fabs(a) > 1e-12 && discr >= 0.0 ) {
// 		printf("nonlin_solver: s1 = %.4e\ts2 = %.4e\n",s = (-b + sqrt(discr))/a,s = (-b - sqrt(discr))/a);
		s = (-b + sqrt(discr))/a;
		if( sgn(A*s + B) == sgn(du) ) s = (-b - sqrt(discr))/a;
		if( s >= 0.0 && s <= 1.0 ) {
			s1 = 1 - s;
			xs = vec_difference(x,vec_lin_comb(x0,x1,s,s1));
			sol.g = u0*s + u1*s1 + f*length_vec(xs);
			sol.c = 'y';
		}
	}
	else if( fabs(a) <= 1e-12 && fabs(b) > 1e-12 ) {
		s = -0.5*(c/b);
		if( s >= 0.0 && s <= 1.0 ) {
			s1 = 1 - s;
			xs = vec_lin_comb(x0,x1,s,s1);
			sol.g = u0*s + u1*s1 + f*length_vec(xs);
			sol.c = 'y';
		}
	}			
	
// 	printf("s = %.6e\t ch = %c\n",sol.g,sol.c);
	
	return sol;
}


/*--------------*/

double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,double b,
						struct myvector x1mx0,struct myvector x) {
	double s1 = 1.0 - s,ls,bs;
	struct myvector xs,xmxs;
	
	xs = vec_lin_comb(x0,x1,s,s1);
	xmxs = vec_difference(x,xs);
	ls = length_vec(xmxs);
	
	return u0 - u1 + b*dot_product(xmxs,x1mx0)/ls;
	
}
			



/**********************************************/
/*** linear algebra ***/

double angle(double x,double y) {
  double ang;
  if( y >= 0.0 ) ang = acos(x/sqrt(x*x + y*y));
  else ang = 2.0*PI - acos(x/sqrt(x*x + y*y));
  return ang;
} 

double length(double x,double y) {
  return sqrt(x*x+y*y);
}

double length_vec(struct myvector x) {
  return sqrt(x.x*x.x+x.y*x.y);
}

struct mymatrix matrix_inverse(struct mymatrix matr) {
  struct mymatrix mi;
  double rdet;

  rdet=1.0/(matr.a11*matr.a22-matr.a12*matr.a21);
  mi.a11=matr.a22*rdet;
  mi.a12=-matr.a12*rdet;
  mi.a21=-matr.a21*rdet;
  mi.a22=matr.a11*rdet;
  return mi;
}

struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b) {
  struct mymatrix c;

  c.a11=a.a11*b.a11+a.a12*b.a21;
  c.a12=a.a11*b.a12+a.a12*b.a22;
  c.a21=a.a21*b.a11+a.a22*b.a21;
  c.a22=a.a21*b.a12+a.a22*b.a22;
  return c;
}

struct mymatrix transpose(struct mymatrix matr) {
  struct mymatrix mt;

  mt.a11=matr.a11;
  mt.a22=matr.a22;
  mt.a12=matr.a21;
  mt.a21=matr.a12;
  return mt;
}

struct myvector matr_vec(struct mymatrix matr,struct myvector vec) {
  struct myvector c;

  c.x=matr.a11*vec.x+matr.a12*vec.y;
  c.y=matr.a21*vec.x+matr.a22*vec.y;
  return c;
}

double dot_product(struct myvector a,struct myvector b) {
   return a.x*b.x+a.y*b.y;
}

struct myvector ga_plus_b(double lam,struct myvector a,struct myvector b) {
  struct myvector c;

  c.x=lam*a.x+b.x;
  c.y=lam*a.y+b.y;
  return c;
}

double solve_quadratic(double a,double b,double c) {
  double discr,sol;

  discr=b*b-4.0*a*c;
  if( discr < 0.0 ) sol=INFTY;
  else sol=0.5*(-b+sqrt(discr))/a;
  return sol;
}

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b) {
	struct myvector v;
	
	v.x = a*v1.x + b*v2.x;
	v.y = a*v1.y + b*v2.y;
	
	return v;
}
	
struct myvector vec_difference(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	
	return v;
}
	
struct myvector vec_sum(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	
	return v;
}

/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(long ind) {
  long loc, ptemp;
  long indp, indc;
  char ch;

//  printf("addtree(%li.%li)\n",ind%NX,ind/NX);
  count++;
  tree[count]=ind;
  pos[ind]=count;
  if( count > 1 ) {
    loc=count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( g[indc] < g[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( g[indc] < g[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }  
}

/*------------------------------------------------------------------*/

void updatetree(long ind) {
  long loc, lcc;
  double g0,g1,g2;

//  printf("updatetree(%li.%li)\n",ind%NX,ind/NX);

  g0=g[ind];
  loc=pos[ind];
  while( loc > 1 && g0 < g[tree[loc/2]] ) {
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }  
  g1=g[tree[loc*2]];
  g2=g[tree[loc*2+1]];
  lcc=count;
  while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind; 
    pos[tree[loc]]=loc;
  }
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
void deltree() {
  long loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
  char chd, ch='n';;

//  printf("deltree(%li.%li)\n",tree[1]%NX,tree[1]/NX);

  mind=tree[1];
  pos[tree[1]]=0;
  tree[1]=tree[count];
  pos[tree[1]]=1;
  count--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < count )  {
    ic1=tree[lcc];
    ic2=tree[lcc+1];
    if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
      if( (g[ic1]) <= (g[ic2]) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=tree[lcc];
    if( (g[ind]) > (g[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {    
    ptemp=pos[ind];
    pos[ind]=pos[ic];
    tree[loc]=ic;
    pos[ic]=ptemp;
    tree[lcc]=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
        if( (g[ic1]) <= (g[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%li)=%li, t1=%.12e\n",ic1,lcc,g[ic1]);
      if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
}


/********************************************************/		    
/*** main ***/

 int main() {
  long i,j,ind,k; 
  double dd,errmax = 0.0,erms = 0.0;
  clock_t CPUbegin;
  double cpu;
  FILE *fg, *fu, *ft;
  
  param();
  CPUbegin=clock();
//  initial_curve();
  ipoint();
  fast_marching_method();
  cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
  printf("cputime = %g\n",cpu);
  /* compute the intencity of the reactive current */
  fg=fopen(f_output_name,"w");
  ind=0;
  k = 0;
  for( j=0; j<NY; j++ ) {
    for( i=0; i<NX; i++ ) {
    	if( ms[ind] < 2 ) g[ind] = INFTY;
    	else {
    		dd = fabs(g[ind] - Uexact[ind]);
    		errmax = max(errmax,dd);
    		erms += dd*dd;
    		k++;
    	}
      fprintf(fg,"%.12e\t",g[ind]);
      ind++;
    }
    fprintf(fg,"\n");
  }
  fclose(fg);
  printf("NX = %i, NY = %i, errmax = %.4e, erms = %.4e\n",NX,NY,errmax,sqrt(erms/k));
  return 0;
}  
