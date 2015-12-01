NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
source tools/qa/deps/common.sh
if [ -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages" ]; then
    echo -e "${GREEN}Activating ${NAMEVER}${RESET}"
    export PYTHONPATH=${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages:${PYTHONPATH}
fi
