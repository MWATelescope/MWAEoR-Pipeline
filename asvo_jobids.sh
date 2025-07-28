#!/bin/bash
find /data/curtin_mwaeor/asvo -name '*.metafits' | while read -r path; do
    export dirname=$(dirname $path) metafits=$(basename $path);
    export jobid=$(basename $dirname) obsid=${metafits%.metafits};
    [ ! -d /data/curtin_mwaeor/nfdata/$obsid ] && mkdir -p /data/curtin_mwaeor/nfdata/$obsid
    [ ! -d /data/curtin_mwaeor/nfdata/$obsid/raw ] && ln -s /data/curtin_mwaeor/asvo/$jobid /data/curtin_mwaeor/nfdata/$obsid/raw
done