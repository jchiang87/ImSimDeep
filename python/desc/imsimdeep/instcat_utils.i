// -*- c++ -*-
%define instcat_utils_DOCSTRING
"Swig-exposed functions for ImSimDeep"
%enddef

%feature("autodoc", "1");
%module(package="instcat_utils", docstring=instcat_utils_DOCSTRING) instcat_utils

%include "lsst/p_lsstSwig.i"
%lsst_exceptions()

%{
#include "lsst/afw.h"
#include "desc/imsimdeep/instcat_utils.h"
%}

%include "desc/imsimdeep/instcat_utils.h"
