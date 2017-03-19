#!/bin/bash

sh testScript.sh | tee testScript.log
diff testScript.log baseline.log

RES=$?

if [ ${RES} != 0 ]; then
    echo "Test failed"
    exit ${RES}
else
    echo "Test passed"
fi

exit 0

