#!/bin/bash

#
# This file contains function(s) to build the doxygen documentation
#

export cmd_DOXYGEN="doxygen"
export cmd_DOT="dot"


#
# This function checks that we are in the root directory and that
# doxygen and dot are available. 
#
checks() {
  if test ! -d source -o ! -d include -o ! -d doc ; then
    echo "*** This script must be run from the top-level directory."
    exit 1
  fi

  if ! [ -x "$(command -v ${cmd_DOXYGEN})" ]; then
    echo "***"
    echo "***   No doxygen program found. Install form your package manager."
    echo "***"
    exit 1
  fi
  
  if ! [ -x "$(command -v ${cmd_DOT})" ]; then
    echo "***"
    echo "***   No dot program found (part of graphviz package). Install form your package manager."
    echo "***"
    exit 1
  fi
}


#
# Make the docs.
#
build_documentation()
{
  cd doc/
  file="Doxyfile"
  pwd
  echo ${file}
  "${cmd_DOXYGEN}" "${file}"
}
