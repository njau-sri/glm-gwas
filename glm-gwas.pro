TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -lgfortran -lquadmath
QMAKE_CXXFLAGS += -std=c++11 -fopenmp
QMAKE_LFLAGS += -static -fopenmp

SOURCES += \
        main.cpp \
    glm_gwas.cpp \
    cmdline.cpp \
    lapack.cpp \
    lsfit.cpp \
    pheno.cpp \
    statsutil.cpp \
    vcf.cpp

HEADERS += \
    cmdline.h \
    lapack.h \
    lsfit.h \
    pheno.h \
    split.h \
    statsutil.h \
    vcf.h
