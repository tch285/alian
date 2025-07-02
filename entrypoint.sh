#!/usr/bin/bash

. /etc/profile.d/z00_lmod.sh
/opt/yasp/yaspenv.sh "module load bundle/hepbase heppyy alian && $*"