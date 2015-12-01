NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
source tools/qa/deps/common.sh
if [ -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib" ]; then
    echo -e "${GREEN}Activating ${NAMEVER}${RESET}"
    export LD_LIBRARY_PATH=${QAWORKDIR}/depinstall/${NAMEVER}/lib:${LD_LIBRARY_PATH}
    export CPATH=${QAWORKDIR}/depinstall/${NAMEVER}/include:${CPATH}
    # Copy file for documentation generator
    mkdir -p depends/${NAMEVER}/src/
    cp ${QAWORKDIR}/depinstall/${NAMEVER}/funcs_key.c depends/${NAMEVER}/src/
fi
