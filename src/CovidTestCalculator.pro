#-------------------------------------------------
#
# Project created by QtCreator 2020-08-26T18:22:28
#
#-------------------------------------------------

QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = CovidTestCalculator
TEMPLATE = app

INCLUDEPATH += ../../eigen

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    simulation.cpp

HEADERS += \
        mainwindow.h \
    simulation.h

FORMS +=

RESOURCES +=
