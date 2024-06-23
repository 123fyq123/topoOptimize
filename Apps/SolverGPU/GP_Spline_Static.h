#pragma once
#include <assert.h>
#include <fstream>
#include <memory>
#include <thread>

#include <cuda_profiler_api.h>

#include "GP_Attribute.h"
#include "GP_Handle.h"
#include "Kernels/GP_For_Each.cuh"
#include "GP_Launch_Box.h"
#include "GP_Types.h"
#include "Util/GP_Report.h"
#include "GP_Debug.h"
#include "Util/GP_Import_Obj.h"
#include "Util/GP_Log.h"
#include "Util/GP_Timer.h"

#include <filesystem>

#include <thread>

#pragma once
#include <assert.h>
#include <fstream>
#include <memory>

#include "GP_Spline.h"

#include "GP_IGA_kernel.cuh"

namespace gpmesh {

class GPSplineSolidStatic : public GPSplineSolid
{
   public:
    GPSplineSolidStatic(const GPSplineSolidStatic&) = delete;

    // debug part
    explicit GPSplineSolidStatic(const bool quite = false) : GPSplineSolid() {};

    void debug_IGA_test()
    {

    }
};

}