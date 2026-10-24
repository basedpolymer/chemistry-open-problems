#!/bin/bash
# helper for tar --to-command: decompress one census member (alltrans<n>_k=<k>.gz) and test it with kscheck -t
# writes one summary line per member to $CENSUS_LOG and survivors to $CENSUS_OUT
name="$TAR_FILENAME"
res=$(gzip -dc | "$KSCHECK" -t 2>&1 >> "$CENSUS_OUT" | head -1)
echo "$name $res" >> "$CENSUS_LOG"
