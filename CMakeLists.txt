cmake_minimum_required(VERSION 2.6)

PROJECT(HighOrder)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

ADD_PARAVIEW_PLUGIN(
    HighOrder
    "0.1"
    SERVER_MANAGER_XML HighOrder.xml
    SERVER_MANAGER_SOURCES vtkHighOrder.cxx
    GUI_RESOURCES HighOrder.qrc
    GUI_RESOURCE_FILES HighOrderGUI.xml
)
