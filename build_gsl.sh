curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd gsl-2.3
./configure --prefix=$PWD/../gsl_local/
make -j4
make install
