voms-proxy-destroy;
kx509;
voms-proxy-init -rfc --voms=fermilab:/fermilab/minerva/Role=Analysis --noregen -valid 24:0;
#voms-proxy-init -rfc -noregen;
export EXPERIMENT=minerva
export IFDH_DEBUG=0
