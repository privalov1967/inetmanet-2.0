// Calculates Input M/G/Infinity Source parameters from discrete traffic 
// statistics, which this source will try to simulate. 


#include <omnetpp.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

class MGI_LOG_Par : public cSimpleModule
{
  protected:
    simtime_t timeslot;      // time between packets
    int cellsize;            // cell size in bytes
    
    cXMLElement *doc;        // XML file with traffic statistics    

    // Fields and function to calculate ACF shape parameters
    double H;
    std::vector<double>  nacf;
    double differ(double x, double y);
    double amoeba2D(double F[3],double X[3], double Y[3],
                    double ftol,int *nfunk );
    double amotry2d(double F[3], double X[3], double Y[3],
                    double c);
    
    // Function to calculate MGI individual source parameters
    std::vector<double>  dst;
  
    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);
};

// The module class needs to be registered with OMNeT++
Define_Module(MGI_LOG_Par);


#define gDEBP(A,B) (A->getDocumentElementByPath(A,B))->getNodeValue()


void MGI_LOG_Par::initialize()
{
    // Initialize is called at the beginning of the simulation.

    EV<<"Source parameters calculation started"<<endl;

    // traffic data input
    doc=par("data");
    timeslot = atof(gDEBP(doc,"/TIMESLOT"));
    cellsize = atoi(gDEBP(doc,"/CELLSIZE"));     
    
    H        = atof(gDEBP(doc,"//HURST"));

    cStringTokenizer tok1(gDEBP(doc,"//ACF-VALUES"));
    nacf = tok1.asDoubleVector();
  
    cStringTokenizer tok2(gDEBP(doc,"/TRAFFIC-DISTRIBUTION"));
    dst  = tok2.asDoubleVector();

    double av;
    av = atof(gDEBP(doc,"//CELLSPERSLOT"));
   
    // MGI parameters calculations and output
    std::ofstream out( par("profile").stringValue() );
    out<<"<?xml version=\'1.0' ?>"<<endl<<endl;
    out<<"<MGI_PARAMETERS>"<<endl;
    out<<" <TIMESLOT> " << timeslot << " </TIMESLOT>"<<endl;
    out<<" <CELLSIZE> " << cellsize << " </CELLSIZE>"<<endl;

    out<<" <TRAFFIC-MEAN> " << av << " </TRAFFIC-MEAN>"<<endl;
    
    // Normalized ACF shape parameters calculation
    double X[3]={0.5, 0.5, 0.1}, Y[3]={0.1, 0.5, 0.5};
    double f[3], tol, alpha=4-2*H;
    int count=0;
    f[0]= differ(X[0], Y[0]);
    f[1]= differ(X[1], Y[1]);
    f[2]= differ(X[2], Y[2]);
    tol = amoeba2D(f, X, Y, 1.e-8, &count );
    double px=X[0], py=Y[0];
    
    int i;
    double s, xx;    
    for(s=py,i=2; (xx=pow(px+i,-alpha))/s > 1.0e-9; s+=xx,i++);
    int imax=i; 
    double A=(1-py)/(s-py);
    
    out<<" <ACF-SHAPE-PARAMETERS>  ";
    // alpha - tail power, px - tail shift, py - prob of 1
    out<< alpha <<"  "<< px <<"  "<< py <<"  "<< A <<"  "<<imax;  
    out<<"  </ACF-SHAPE-PARAMETERS>";

    // MGI individual sources parameters output
    out<<" <INDIVIDUAL-SOURCES>"<<endl;
    for(s=py,i=2; i<=imax; i++) s+=i*A*pow(px+i,-alpha);
    out<<"  <LIFETIME-MEAN> "<< s << " </LIFETIME-MEAN>"<<endl;
    out<<"  <SOURCE-BORN-RATE> "<< av/s;
    out<<" </SOURCE-BORN-RATE>"<<endl;

    out<<"  <MAXIMAL-SPEED> "<< 1 <<" </MAXIMAL-SPEED>"<<endl;
    out<<"  <SPEED-DISTRIBUTION> "<<endl<<"0 1"<<endl;
    out<<endl<<"  </SPEED-DISTRIBUTION>"<<endl;

    out<<" </INDIVIDUAL-SOURCES>"<<endl;
    
    out<<"</MGI_PARAMETERS>"<<endl;
    out.close();

    EV<<"Source parameters calculation complete"<<endl;      
}

void MGI_LOG_Par::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
    delete msg;
}


/* --------------------------------------------------------------- */
#define TINY 1.0e-10  /* A small number. */
#define NMAX 5000
#define SWAP(i,j) {(swap)=F[i];F[i]=F[j];F[j]=(swap);\
(swap)=X[i];X[i]=X[j];X[j]=(swap);(swap)=Y[i];Y[i]=Y[j];Y[j]=(swap);}

double MGI_LOG_Par::amotry2d(double F[3], double X[3], double Y[3],
                double c)
{
   double Ftry, Xtry, Ytry;

   Xtry=(1-c)*0.5*(X[0]+X[1]) + c*X[2];
   Ytry=(1-c)*0.5*(Y[0]+Y[1]) + c*Y[2];
   Ftry=differ(Xtry,Ytry);
   if( Ftry<F[2] ) {F[2]=Ftry; X[2]=Xtry; Y[2]=Ytry;}
   return Ftry;
}


double MGI_LOG_Par::amoeba2D(double F[3],double X[3], double Y[3],
                               double ftol,int *nfunk )
{
double swap, rtol, Ftry, Fsave;

*nfunk=0;
while(1) {

   /* Sorting to find best (F[0]), next to worst (F[1]) and worst (F[2])
      simplex point */
   if( F[0]>F[1] ) SWAP(0,1)
   if( F[1]>F[2] ) SWAP(1,2)
   if( F[0]>F[1] ) SWAP(0,1)
   /* Compute the fractional range from highest to lowest and return
      if satisfactory (return also if too many iterations made). */
   rtol=2.0*fabs(F[2]-F[0])/(fabs(F[2])+fabs(F[0])+TINY);
   if (rtol < ftol || *nfunk >= NMAX) return rtol;

   /* Begin a new iteration. First extrapolate by a factor -1 through
      the face of the simplex across from the high point, i.e., reflect
      the simplex from the high point. */
   Ftry = amotry2d(F,X,Y, -1);  (*nfunk)++;
   /* Gives a result better than the best point, so try an additional
      extrapolation by a factor 2. */
   if (Ftry <= F[0]) { Ftry = amotry2d(F,X,Y, 2);  (*nfunk)++; }
   else if (Ftry >= F[1]) {
      /* The reflected point is worse than the second-highest, so look for
         an intermediate lower point, i.e., do a one-dimensional contraction. */
      Fsave= F[2];
      Ftry = amotry2d(F,X,Y, 0.5);  (*nfunk)++;
      if (Ftry >= Fsave) {
         /* Can’t seem to get rid of that high point.
            Better contract around the lowest (best) point. */
         F[1]=differ( X[1]=0.5*(X[0]+X[1]), Y[1]=0.5*(Y[0]+Y[1]) );
         F[2]=differ( X[2]=0.5*(X[0]+X[2]), Y[2]=0.5*(Y[0]+Y[2]) );
         *nfunk+=2;
         }
      }
   }  /* end of while(1) */
}


#define FIT_POINTS 4
#define IMAX 100

double MGI_LOG_Par::differ(double x, double y)
{
int i, k, n, imax;
double h, s, s2, ss[FIT_POINTS+1], nr[FIT_POINTS+1], alpha=4-2*H;

for(k=0; k<=FIT_POINTS; k++) ss[k]=nr[k]=0;

for(s=h=pow(2+x,-alpha),n=3; h/s>1.e-6; n++) s+=h=pow(n+x,-alpha);
imax=n;
s2=s;

for(ss[1]=s, i=2; i<=IMAX && i<imax; i++) {
   s-=pow(i+x,-alpha);
   for(k=1; k<=FIT_POINTS; k++) if(i>k) ss[k]+=s;
   }  
for(h=0, k=1; k<=FIT_POINTS; k++) {
   nr[k]=(1-y)*ss[k]/(s2+(1-y)*ss[1]);
   h+=(nacf[k]-nr[k])*(nacf[k]-nr[k]);
   }  

return h; 
}
