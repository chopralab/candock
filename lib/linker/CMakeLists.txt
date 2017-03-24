include_directories($CMAKE_CURRENT_SOURCE_DIR)


if(CUDA_FOUND)

set( linker_SRC
        dockedconformation.cpp
        dockedconformation.hpp
        genericlinker.cpp
        cuda_iterativelinker.cpp
        iterativelinker.cpp
        linker.cpp
        linker.hpp
        partial.cpp
        partial.hpp
        poses.cpp
        poses.hpp
        seed.cpp
        seed.hpp
        segment.cpp
        segment.hpp
        state.cpp
        state.hpp
        staticlinker.cpp
        cuda_linker.cu
)

include_directories( ../modeler )

cuda_add_library(linker OBJECT ${linker_SRC})


else(CUDA_FOUND)

set( linker_SRC
        dockedconformation.cpp
        dockedconformation.hpp
        genericlinker.cpp
        cuda_iterativelinker.cpp
        iterativelinker.cpp
        linker.cpp
        linker.hpp
        partial.cpp
        partial.hpp
        poses.cpp
        poses.hpp
        seed.cpp
        seed.hpp
        segment.cpp
        segment.hpp
        state.cpp
        state.hpp
        staticlinker.cpp
)

include_directories( ../modeler )

add_library(linker OBJECT ${linker_SRC})

endif(CUDA_FOUND)
