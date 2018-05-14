#!/bin/bash
# This file is for pushing data into the cluster
# Please, to send a file use:
#   $ ./transfer.sh push local files.rar
#
# For sending use:
#   $ ./transfer.sh pull local files.rar
#
# If connected from outside of Yachay Tech please use:
#   $ ./transfer.sh pull pub files.rar
#

if [ $1 == "push" ]
    then
        scp -P 22022 $3 190.15.128.35:/home/joshua/test/
    else
        if [ $1 == "pull" ]
            then
                scp -P 22022 190.15.128.35:/home/joshua/test/$3 .
        fi
fi
