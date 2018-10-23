TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += D:/boost_1_67_0

LIBS += -Ld:/mingw-w64-i686-lapack-3.8.0-3-any.pkg/mingw32/lib -llapack -lblas -lgfortran -lquadmath

QMAKE_LFLAGS += -static

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
