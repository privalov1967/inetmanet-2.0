// Packet traffic discretization and statistical characteristics of 
// discrete traffic calculaton, such as probability distribution,
// normalized autocorrelation function and Hurst parameter.
// Hurst parameter is calculated by partial variations method.

#include <omnetpp.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

class Traffic_Stats : public cSimpleModule
{

  protected:
    // Basic parameters
    double timeslot;
    int cellsize;
    char* outFileName;  

    // Fields for traffic discretization
    simtime_t currentSlotBegin, currentSlotEnd;
    long int cellsArrived;

    // Fields to control discrete traffic output as vector
    bool SaveDiscreteTraffic;
    cOutVector traffic;     

    // Field for PDF calculation
    cDoubleHistogram dst;
    int numHistCells;

    // Fields for ACF calculation
    unsigned int num_of_points;
    std::vector<cStdDev> covar;
    std::deque<long int> traff;
    std::vector<double>  nacf;

    // Fields for Self-Similar parameters calculation
    int num_of_sizes;
    int minSize, maxSize;
    std::vector<cStdDev>  meansqr;
    long int *blk_size;  
    long int *counter;
    double   *sum;
    double H;             // Hurst parameter  
    
    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);
    virtual void finish();
};

// The module class needs to be registered with OMNeT++
Define_Module(Traffic_Stats);


void Traffic_Stats::initialize()
{
    // Initialize is called at the beginning of the simulation.
    timeslot = par("timeslot");
    cellsize = par("cellsize");

    SaveDiscreteTraffic = par("SaveDiscreteTraffic");
    if(SaveDiscreteTraffic) traffic.setName("Traffic");
    
    // For PDF calculation
    numHistCells=par("numHistCells");
    dst.setNumCells(numHistCells);
    dst.setRange(-0.5, -0.5+numHistCells);

    // For ACF calculation
    num_of_points=(unsigned int)par("num_of_points");
    covar.resize(num_of_points);
    nacf.resize(num_of_points);
    
    // For Self-similar parameters calculation
    num_of_sizes = par("num_of_sizes");
    minSize = par("minSize");
    maxSize = par("maxSize");
    meansqr.resize(num_of_sizes);
    counter = new long[num_of_sizes];
    sum     = new double[num_of_sizes];
    blk_size= new long[num_of_sizes];
    for(int i=0; i<num_of_sizes; i++) sum[i]=counter[i]=0;    
    blk_size[0] = minSize; blk_size[num_of_sizes-1] = maxSize;
    for(int i=1; i<num_of_sizes-1; i++) 
       blk_size[i] = (long int)ceil(exp(((num_of_sizes-1-i)*log(minSize)+i*log(maxSize))/(num_of_sizes-1)));

    currentSlotBegin  =0; 
    currentSlotEnd    =(simtime_t)timeslot;
    cellsArrived      =0;
}


void Traffic_Stats::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
    cPacket *pkt=check_and_cast<cPacket *>(msg);
    simtime_t curr_time = simTime();

    // traffic discretization and processing
    if( curr_time >= currentSlotEnd ) 
       do{
          if(SaveDiscreteTraffic) traffic.record(cellsArrived);
          
          // for PDF calculation 
          dst.collect(cellsArrived);

          // for ACF calculation
          traff.push_back(cellsArrived);
          if( traff.size() == num_of_points ) {
             for(unsigned int i=0; i<num_of_points; i++) 
                covar[i].collect(traff[0]*traff[i]);
             traff.pop_front();
             }

          // for Self-similar parameters calculation
          for(int i=0; i<num_of_sizes; i++) {
             counter[i]++; sum[i]+=cellsArrived;
             if( counter[i] == blk_size[i] ) {
                 meansqr[i].collect((sum[i]/blk_size[i])*(sum[i]/blk_size[i]));
                 sum[i]=counter[i]=0;
                 }
             }
           
          currentSlotBegin = currentSlotEnd;
          currentSlotEnd  += (simtime_t)timeslot;
          cellsArrived     = 0; 
       } while(curr_time >= currentSlotEnd);

    cellsArrived += (long int)ceil(1.0*pkt->getByteLength()/cellsize);     

    send(msg, "out");
}


void Traffic_Stats::finish()
{
    double a, b;
    int i, imax;
    
    // results output begin
    std::ofstream out( par("outFileName").stringValue() );
    out<<"<?xml version=\'1.0' ?>"<<endl<<endl;
    out<<"<STATISTICS>"<<endl;
    out<<" <TIMESLOT> " << timeslot << " </TIMESLOT>"<<endl;
    out<<" <CELLSIZE> " << cellsize << " </CELLSIZE>"<<endl;
    out<<" <TRAFFIC-MEAN>"<<endl;
    out<<"  <CELLSPERSLOT> "<< dst.getMean() << " </CELLSPERSLOT>"<<endl;
    out<<"  <BYTESPERSEC> "<< dst.getMean()*cellsize/timeslot <<" </BYTESPERSEC>"<<endl; 
    out<<" </TRAFFIC-MEAN>"<<endl;
    out<<" <TRAFFIC-VARIANCE> "<< dst.getVariance() <<" </TRAFFIC-VARIANCE>"<<endl;

    for(imax=dst.getNumCells()-1; imax>=0; imax--) 
       if( dst.getCellPDF(imax)!=0 ) break;
    out<<" <TRAFFIC-MAXVALUE> "<< imax << " </TRAFFIC-MAXVALUE>"<<endl;

    out<<" <TRAFFIC-DISTRIBUTION>"<<endl;
    for(i=0; i<=imax; i++)  
       out << dst.getCellPDF(i)<<"  ";
    out<<endl<<" </TRAFFIC-DISTRIBUTION>"<<endl;
    out<<" <OUT-OF-RANGE> "<<dst.getOverflowCell()+dst.getUnderflowCell()<<" </OUT-OF-RANGE>"<<endl;
    
    out<<" <NORM-ACF>"<<endl;
    out<<"  <ACF-RANGE> "<< num_of_points <<" </ACF-RANGE>"<<endl;
    out<<"  <ACF-VALUES> "<<endl;
    double ACF0 = covar[0].getMean() - dst.getMean()*dst.getMean();
    for(i=0; i<(int)num_of_points; i++) 
       out<< (nacf[i]=(covar[i].getMean() - dst.getMean()*dst.getMean())/ACF0) <<"  ";  
    out<<endl<<"  </ACF-VALUES> "<<endl;
    out<<" </NORM-ACF>"<<endl;
       
    // Regression for Self-similar parameters calculation
    double variance[num_of_sizes];
    double sigma[num_of_sizes];
    int good_point=num_of_sizes;
    cStdDev x;
    cStdDev y;
    cStdDev xy;
    cStdDev s2;
    for(int i=0; i<num_of_sizes; i++) 
       if(meansqr[i].getCount()<32) { good_point=i; break; } 
    for(i=0; i<good_point; i++) {
       variance[i]=meansqr[i].getMean() - dst.getMean()*dst.getMean();
       x.collect(log(blk_size[i]));
       y.collect(log(variance[i]));
       xy.collect(log(blk_size[i])*log(variance[i]));
       }
    a = (xy.getMean() - x.getMean()*y.getMean())/x.getVariance();
    b = y.getMean() - a*x.getMean();
    H = 1-fabs(a)/2;

    out<<" <SELFSIM>"<<endl;
    out<<"  <HURST> "<< H <<" </HURST> "<<endl;
    out<<"  <A> "<< a <<" </A> "<<endl;
    out<<"  <B> "<< b <<" </B> "<<endl;
    for(i=0; i<good_point; i++) 
       s2.collect((sigma[i]=log(variance[i]) - a * log(blk_size[i]) - b)*sigma[i]);
    out<<"  <SIGMA-SQUARE> "<< s2.getMean()<< " </SIGMA-SQUARE>"<<endl;
    out<<"   <REGRESSION-POINTS>"<<endl;
    out<<"    <POINTS-NUMBER> "<< good_point << " </POINTS-NUMBER>" << endl;
    out<<"    <POINTS-LEGEND> " <<"Blksize  Variance  Sigma"<<" </POINTS-LEGEND>"<<endl;
    for(i=0; i<good_point; i++) {
       out<<"    <POINT> ";
       out<<blk_size[i]<<"  "<<variance[i]<<"  "<<sigma[i];
       out<<"    </POINT>"<<endl;
       }
    delete[] counter; delete[] sum; delete[] blk_size;
    out<<"  </REGRESSION-POINTS>"<<endl;
    out<<" </SELFSIM>"<<endl;
     
    out<<"</STATISTICS>"<<endl;

    out.close();
    EV<<"All statistical calculations complete";    
}
