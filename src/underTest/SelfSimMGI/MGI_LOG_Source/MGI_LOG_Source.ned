// Source of Input M/G/Infinity type of selfsimilar traffic. 
// Parameters of generated traffic are given in xml file. 
// Also there is possibility to calculate these parameters from 
// discrete traffic statistics, which this Input M/G/Infinity Source 
// will try to simulate.
// ATTENTION! : Calculation of these parameters in initializaton phase 
// is very time consuming! Please, be patient!

package MGI_LOG_Source;

simple MGI_LOG_Par
{
    parameters:
        xml data = default(xmldoc("traffic_statistics.xml"));
        string profile = default("mgi_log_profile.xml");
}

simple Ind_MGI_LOG_Source
{
    parameters:
        xml profile = default(xmldoc("mgi_log_profile.xml"));
}


simple MGI_LOG_Collector 
{
    gates:
        input  in @directIn;
        output out;
}


simple MGI_LOG_Generator
{
    parameters:
        xml profile = default(xmldoc("mgi_log_profile.xml"));
}
      

module MGI_LOG_Source
{ 
    parameters:
        bool calcRequired = default(true);
    gates:
        output out;
    submodules:
        calc[calcRequired ? 1:0] : MGI_LOG_Par;
        newsrcs[calcRequired ? 0:1] : MGI_LOG_Generator;
        collect : MGI_LOG_Collector;
    connections:
        collect.out --> out;       
}