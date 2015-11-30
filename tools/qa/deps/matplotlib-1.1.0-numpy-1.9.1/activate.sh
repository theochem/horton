NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
source tools/qa/deps/common.sh
if [ -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages" ]; then
    echo -e "${COLOR}Activating ${NAMEVER}${RESET}"
    export PYTHONPATH=${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages:${PYTHONPATH}
    export MATPLOTLIBRC=${QAWORKDIR}
    echo "backend: agg" > $MATPLOTLIBRC/matplotlibrc
fi
