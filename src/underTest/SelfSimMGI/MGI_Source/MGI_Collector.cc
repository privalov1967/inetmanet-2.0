// Collects packets from all individual sources into one output stream

#include <omnetpp.h>

class MGI_Collector : public cSimpleModule
{
  protected:
    
    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);

};

// The module class needs to be registered with OMNeT++
Define_Module(MGI_Collector);

void MGI_Collector::initialize()
{
    // Initialize is called at the beginning of the simulation.
    
}

void MGI_Collector::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
    send(msg,"out");
}
