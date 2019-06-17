QT += core
QT -= gui

CONFIG += c++11

TARGET = HRBFQI
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

DEFINES = "THREADS_NUM=4"

HEADERS += DataStructure/*.h \
           HRBFQI/FileManager.h \
           HRBFQI/HRBF.h \
           HRBFQI/MeshCleaner.h \
           numericalC/*.h \
           numericalC/*.H \
           polygonizer/*.h \
           PQP/include/*.h \
           PQP/src/*.h

SOURCES += HRBFQI/HRBFQI.cpp \
           HRBFQI/FileManager.cpp \
           HRBFQI/HRBF.cpp \
           HRBFQI/MeshCleaner.cpp \
           DataStructure/AxisKdTree.cpp \
           DataStructure/OctTree.cpp \
           DataStructure/PointSet.cpp \
           DataStructure/PolygonalMesh.cpp \
           polygonizer/ImplicitFunction.cpp \
           polygonizer/Polygonizer.cpp \
           PQP/src/Build.cpp \
           PQP/src/BV.cpp \
           PQP/src/PQP.cpp \
           PQP/src/TriDist.cpp \

QMAKE_CXXFLAGS += -std=c++11 -fpermissive -fopenmp
QMAKE_LFLAGS +=  -fopenmp
