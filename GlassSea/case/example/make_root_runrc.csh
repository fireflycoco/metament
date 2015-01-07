#!/bin/csh -f

set RUNRC = "${G3S}/runpool/${RUNNAME}/G3S_runrc.csh"

cat > ${RUNRC} << EOF
setenv G3S "$G3S"
setenv CASENAME "$CASENAME"
setenv RUNTIME "$RUNTIME"
setenv RUNNAME "$RUNNAME"

#TODO: add other env vars

EOF

chmod +x ${RUNRC}
