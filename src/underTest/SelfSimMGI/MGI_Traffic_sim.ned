import Common_Modules.*;
import MGI_Source.*;
import Statistics.*;

network MGI_Traffic_sim
{
    submodules:
        src:   MGI_Source;
        stat:  Traffic_Stats;
        dst:   Sink;
    connections:
        src.out  --> stat.in;
        stat.out --> dst.in;
}

