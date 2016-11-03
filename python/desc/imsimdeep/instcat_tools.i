// -*- c++ -*-
%define instcat_tools_DOCSTRING
"Swig-exposed functions for ImSimDeep"
%enddef

%feature("autodoc", "1");
%module(package="instcat_tools", docstring=instcat_tools_DOCSTRING) instcat_tools

%include "lsst/p_lsstSwig.i"
%lsst_exceptions()

%{
#include "lsst/afw.h"
#include "desc/imsimdeep/InstcatTools.h"
%}

%include "desc/imsimdeep/InstcatTools.h"
