cmake_minimum_required(VERSION 2.6)

PROJECT(HighOrder)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

add_custom_command(OUTPUT shapeFunctions.h
    COMMAND bunzip2 -k shapeFunctions.h.bz2
    DEPENDS shapeFunctions.h.bz2
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Unziping shapeFunctions.h"
    VERBATIM
    )

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

ADD_PARAVIEW_PLUGIN(
    HighOrder
    "0.1"
    SERVER_MANAGER_XML HighOrder.xml
    SERVER_MANAGER_SOURCES vtkHighOrder.cxx shapeFunctions.h
    GUI_RESOURCES HighOrder.qrc
    GUI_RESOURCE_FILES HighOrderGUI.xml
)

INSTALL(TARGETS HighOrder DESTINATION lib)
