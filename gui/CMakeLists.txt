set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUI  ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt5Widgets)

qt5_wrap_ui(CANDOCK_GUI_UIS settings.ui)
set(CANDOCK_GUI_CPP main.cpp settings.cpp settings.hpp settings.ui)

# Tell CMake to create the helloworld executable
add_executable(candock_gui ${CANDOCK_GUI_CPP} ${CANDOCK_GUI_UIS})

# Use the Widgets module from Qt 5.
target_link_libraries(candock_gui Qt5::Widgets)

