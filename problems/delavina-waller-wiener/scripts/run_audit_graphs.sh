#!/bin/bash
# Graph-level audit: all connected graphs for d=3,4; for d=5 all connected graphs with at most 15 edges.
# Also aggregates the exhaustive n=11 enumeration into certificates/enumeration-summary.txt.
source ~/envs/sci/bin/activate
P=/mnt/c/Users/demor/chemistry-open-problems/problems/delavina-waller-wiener
L=$P/certificates/audit-graphs.log
: > $L
geng -cq 7 | python $P/scripts/audit_graphs.py 3 >> $L 2>&1
geng -cq 9 | python $P/scripts/audit_graphs.py 4 >> $L 2>&1
geng -cq 11 10:15 | python $P/scripts/audit_graphs.py 5 >> $L 2>&1
cat $L
python $P/scripts/enum_summary.py > $P/certificates/enumeration-summary.txt 2>&1
cat $P/certificates/enumeration-summary.txt
