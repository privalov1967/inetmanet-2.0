// Source of packet traffic restored from packet trace file 
// of TL or TCP format

package FromFile_Source;

simple FromFile_Source
{
    parameters:
        string filetype  = default("tcp");
        string tracefile; 
    gates:
        output out;
}
