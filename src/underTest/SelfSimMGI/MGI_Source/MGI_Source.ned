// Source of Input M/G/Infinity type of selfsimilar traffic. 
// Parameters of generated traffic are given in xml file. 
// Also there is possibility to calculate these parameters from 
// discrete traffic statistics, which this Input M/G/Infinity Source 
// will try to simulate.
// ATTENTION! : Calculation of these parameters in initializaton phase 
// is very time consuming! Please, be patient!

package MGI_Source;

simple MGI_Parameters_Calculator
{
    parameters:
        xml data = default(xmldoc("traffic_statistics.xml"));
        string profile = default("mgi_profile.xml");
}

simple Individual_MGI_Source
{
    parameters:
        xml profile = default(xmldoc("mgi_profile.xml"));
}


simple MGI_Collector 
{
    gates:
        input  in @directIn;
        output out;
}


simple MGI_Generator
{
    parameters:
        xml profile = default(xmldoc("mgi_profile.xml"));
}
      

module MGI_Source
{ 
    parameters:
        bool calcRequired = default(true);
    gates:
        output out;
    submodules:
        calc[calcRequired ? 1:0] : MGI_Parameters_Calculator;
        newsrcs [calcRequired ? 0:1]: MGI_Generator;
        collect : MGI_Collector;
    connections:
        collect.out --> out;       
}