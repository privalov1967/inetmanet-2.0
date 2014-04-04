// Generates new individual sources during Input M/G/Infinity Source 
// working.

#include <omnetpp.h>
#include <fstream>
#include <string>

class MGI_Generator : public cSimpleModule
{

  protected:
    cModuleType *Ind_Source_Factory;
    double lambda;
    simtime_t timeslot;

    cXMLElement *doc;
 
    // The following redefined virtual function holds the algorithm.
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);
};

// The module class needs to be registered with OMNeT++
Define_Module(MGI_Generator);


#define gDEBP(A,B) (A->getDocumentElementByPath(A,B))->getNodeValue()


void MGI_Generator::initialize()
{
    // Initialize is called at the beginning of the simulation.

    // Initializing factory of individual sources
    Ind_Source_Factory = cModuleType::get("MGI_Source.Individual_MGI_Source");
    
    // basic parameters initialization
    doc=par("profile");
    timeslot = atof(gDEBP(doc,"/TIMESLOT"));
    lambda   = atof(gDEBP(doc,"//SOURCE-BORN-RATE"));
    
    // Start message scheduling
    cMessage *msg = new cMessage("newsource");
    scheduleAt( simTime(), msg);
}

void MGI_Generator::handleMessage(cMessage *msg)
{
    // The handleMessage() method is called whenever a message arrives
    // at the module.
    
    // Creating new CBR source and planning next
    Ind_Source_Factory->createScheduleInit("Ind_Src", this->getParentModule());
    scheduleAt(simTime()+timeslot*exponential(1/lambda), msg);
}
