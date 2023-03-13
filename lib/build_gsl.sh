curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd gsl-2.3
./configure --prefix=$PWD/../gsl_local/
make -j4
make install
# Also need to set LD_LIBRARY_PATH env variable.  One way to do this is to put
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/tjhladish/work/AbcSmc/gsl_local/lib
# at the end of .bashrc, or something else that gets sourced reliably.
