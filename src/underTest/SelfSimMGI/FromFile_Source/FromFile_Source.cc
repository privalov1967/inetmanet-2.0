// Source of packet traffic restored from packet trace file 
// of TL or TCP format

#include <omnetpp.h>
#include <string>
#include <fstream>
#include <math.h>

class FromFile_Source : public cSimpleModule
{
  protected:
    // file with traffic trace
    std::ifstream  trace;
    std::string filetype;    

    // fields for packet record processing
    double ArrTime;
    int pktsize;
    int i1, i2, i3, i4;

    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);

    // File type dependent function
    void readRecord();
};

// The module class needs to be registered with OMNeT++
Define_Module(FromFile_Source);


void FromFile_Source::readRecord()
{
    if (filetype == "tcp") 
       trace >> ArrTime >> i1 >> i2 >> i3 >> i4 >> pktsize;
    else   // here it means "TL"
       trace >> ArrTime >> pktsize;
}


void FromFile_Source::initialize()
{
    // Initialize is called at the beginning of the simulation.
    
    filetype=par("filetype").stringValue();    
    trace.open( par("tracefile").stringValue() );
    if (trace != NULL )
       if ( !trace.eof() ) {
           readRecord();
           ArrTime -= floor(ArrTime);
           cMessage *msg = new cMessage("next");
           scheduleAt( (simtime_t)ArrTime, msg);
           }        
       else { 
           trace.close();
           EV<<"Trace file is empty, source is shut down"<<endl;
           }   
    else EV<<"No trace file, source is shut down"<<endl;
}


void FromFile_Source::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
 
    if ( !trace.eof() ) {
        cPacket *pkt = new cPacket("packet");
        pkt->setByteLength(pktsize);
        send(pkt, "out");
        readRecord();
        scheduleAt((simtime_t)ArrTime, msg);
        }        
    else trace.close();    
}
