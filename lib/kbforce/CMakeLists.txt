include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( kbforce_SRC
        openmmapi/src/KBForce.cpp
        openmmapi/include/KBForce.h
        openmmapi/src/KBForceImpl.cpp
        openmmapi/include/KBKernels.h
        openmmapi/include/internal/KBForceImpl.h
        platforms/reference/src/ReferenceKBKernelFactory.cpp 
        platforms/reference/include/ReferenceKBKernelFactory.h
        platforms/reference/src/ReferenceKBKernels.h
        platforms/reference/src/ReferenceKBKernels.cpp
)

include_directories( platforms/reference/include )
include_directories( openmmapi/include )

add_library(kbforce OBJECT ${kbforce_SRC})
