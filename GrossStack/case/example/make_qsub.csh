#!/bin/csh -f

set RUN_ROOT="${G3S}/runpool/${RUNNAME}"
set QSUB_FILE="${RUN_ROOT}/for_qsub.csh"

cat > ${QSUB_FILE} << EOF
#!/bin/csh -f
#PBS -l nodes=1:ppn=8

cd ${RUN_ROOT}
source ./G3S_runrc.csh

\${G3S}/sys/util/RUN >& ${RUN_ROOT}/log.log

EOF

chmod +x ${QSUB_FILE}
