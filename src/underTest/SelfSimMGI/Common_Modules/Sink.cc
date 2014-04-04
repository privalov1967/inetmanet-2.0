// Simplest Sink for arriving packets

#include <omnetpp.h>

class Sink : public cSimpleModule
{
  protected:
   
    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);
};


// The module class needs to be registered with OMNeT++
Define_Module(Sink);


void Sink::initialize()
{
    // Initialize is called at the beginning of the simulation.
    
}


void Sink::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
    
    // cPacket *pkt = check_and_cast<cPacket *>(msg);
    delete msg;    
}
