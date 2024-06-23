#include "GP_Param.h"

std::string get_model_name(int id)
{
    std::string path = STRINGIFY(INPUT_DIR) ;
    if (id == 0)
    {
        //return path + "IGA/konglong.xml";
        //return path + "IGA/H.xml";
        //return path + "IGA/cube324.xml";
        return path + "IGA/konglong.xml";
        //return path + "IGA/clef_1300.xml";
        //for kitty, requires 28119 iterations to get converged result
        //return path + "IGA/kitty_1800.xml";
        //for patch1536, requires 4406 iterations to get converged result
        //return path + "IGA/patches1536.xml";
    }
    else if (id == 1)
    {
        return  path + "IGA/cube324.xml";
    }
}

bool using_pcg()
{
    return true;
}