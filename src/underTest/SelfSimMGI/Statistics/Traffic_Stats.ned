// Packet traffic discretization and statistical characteristics of 
// discrete traffic calculaton, such as probability distribution,
// normalized autocorrelation function and Hurst parameter.
// Hurst parameter is calculated by partial variations method.

package Statistics;

simple Traffic_Stats
{
    parameters:
        double timeslot @unit(s) = default(0.015 ms);   
        int cellsize @unit(byte) = default(100 B);
        
        int numHistCells = default(300);     // for traffic distribution histogramm collection
          
        int num_of_points = default(50);     // for traffic autocorrelation function calculation

        int num_of_sizes = default(8);       // for Hurst parameter calculation
        int minSize = default(12);
        int maxSize = default(3072);

        string outFileName = default("traffic_statistics.xml");

        bool SaveDiscreteTraffic = default(false);
         
    gates:
        input  in;
        output out;
}