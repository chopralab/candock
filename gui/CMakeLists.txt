set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUI  ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt5Widgets)

if(WIN32)
    set(GUI_TYPE WIN32)
elseif (APPLE)
	set(GUI_TYPE MACOSX_BUNDLE)
endif ()

qt5_wrap_ui(CANDOCK_GUI_UIS settings.ui mainwindow.ui)
set(CANDOCK_GUI_CPP main.cpp settings.cpp settings.hpp settings.ui
                             mainwindow.cpp mainwindow.hpp mainwindow.ui)

# Tell CMake to create the helloworld executable
add_executable(candock_gui ${GUI_TYPE} ${CANDOCK_GUI_CPP} ${CANDOCK_GUI_UIS})

# Use the Widgets module from Qt 5.
target_link_libraries(candock_gui Qt5::Widgets        
        candock
        ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES}
)

