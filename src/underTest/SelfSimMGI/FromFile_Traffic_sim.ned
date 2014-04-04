import Common_Modules.*;
import FromFile_Source.*;
import Statistics.*;

network FromFile_Traffic_sim
{
    submodules:
        src:   FromFile_Source;
        stat:  Traffic_Stats;
        dst:   Sink;
    connections:
        src.out   --> stat.in;
        stat.out  --> dst.in; 
}

