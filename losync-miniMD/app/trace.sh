#!/bin/bash

export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so
export EXTRAE_CONFIG_FILE=extrae.xml

$@

