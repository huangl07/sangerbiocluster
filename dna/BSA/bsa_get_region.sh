#!/bin/bash -e

for i in $1/**/pop.*.region
do
    bn=$(basename $i)
    bn=${bn/.region/}
    bn=${bn/pop./}
    REGION=$(readlink -f $i)
    if [ -s ${REGION} ]; then
        echo -e ${bn}'\t'${REGION}
    fi
done
