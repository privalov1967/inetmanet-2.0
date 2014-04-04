import Common_Modules.*;
import MGI_LOG_Source.*;
import Statistics.*;

network MGI_LOG_Traffic_sim
{
    submodules:
        src:   MGI_LOG_Source;
        stat:  Traffic_Stats;
        dst:   Sink;
    connections:
        src.out  --> stat.in;
        stat.out --> dst.in;
}

