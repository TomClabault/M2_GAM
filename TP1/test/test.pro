#-------------------------------------------------
#
# Project created by QtCreator 2018-08-28T10:55:17
#
#-------------------------------------------------

#QT       += core gui opengl
QT += core gui opengl openglwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Mesh_Computational_Geometry
TEMPLATE = app

SOURCES += main.cpp\
    TP1/offreader.cpp \
    TP1/vector.cpp \
        mainwindow.cpp \
    gldisplaywidget.cpp \
    mesh.cpp

HEADERS  += mainwindow.h \
    TP1/color.h \
    TP1/offreader.h \
    TP1/vector.h \
    gldisplaywidget.h \
    mesh.h \
    point.h

FORMS    += mainwindow.ui

#---- Comment the following line on MacOS
#---- Uncomment it on Windows and Linux
# LIBS = -lGLU

#---- Uncomment the following line on Windows
#---- Comment it on Linux and MacOS
LIBS += -lglu32
LIBS += -lOpengl32

