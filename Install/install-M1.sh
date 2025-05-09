#!/bin/bash

export CC=gcc
export CXX=g++

if [ -z $INSTALL_PATH ]; then
    INSTALL_PATH=/usr/local/heplib
fi
export INSTALL_PATH

if [ -z $jn ]; then
    export jn=8
fi
export jn

export CWD=$PWD
export LOGS=$CWD/logs
mkdir -p logs

#================================================================
# brew
#================================================================
#brew install gmp
#brew install qhull
#brew install texinfo

#================================================================
# Install MPFR
#================================================================
export pkg="mpfr-4.2.2"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.gz
cd $pkg
./configure --prefix=$INSTALL_PATH --with-gmp=/opt/homebrew --enable-float128 --enable-thread-safe >>$LOG 2>>$LOG
make -j $jn >>$LOG 2>>$LOG
make install >>$LOG 2>>$LOG
cd $CWD
rm -rf $pkg
echo ""

#================================================================
# Install CLN
#================================================================
export pkg="cln-1.3.7"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.bz2
cd $pkg
./configure --prefix=$INSTALL_PATH --with-gmp=/opt/homebrew >>$LOG 2>>$LOG
make -j $jn >>$LOG 2>>$LOG
make install >>$LOG 2>>$LOG
cd $CWD
rm -rf $pkg
echo ""

#================================================================
# Install GiNaC - Modified version
#================================================================
export pkg="GiNaC"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.gz
cd $pkg
./configure --prefix=$INSTALL_PATH CLN_CFLAGS="-I$INSTALL_PATH/include" CLN_LIBS="-L$INSTALL_PATH/lib -lcln" >>$LOG 2>>$LOG
make -j $jn >>$LOG 2>>$LOG
make install >>$LOG 2>>$LOG
cd $CWD
rm -rf $pkg
echo ""


#================================================================
# Install FLINT
#================================================================
export pkg="flint-3.2.1"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.gz
cd $pkg
./configure CFLAGS="-Wno-incompatible-pointer-types" --enable-static=no --prefix=$INSTALL_PATH --with-gmp=/opt/homebrew --with-mpfr=$INSTALL_PATH >>$LOG 2>>$LOG
make -j $jn >>$LOG 2>>$LOG
make install >>$LOG 2>>$LOG
cd $CWD
rm -rf $pkg
echo ""

#================================================================
# Install JeMalloc
#================================================================
export pkg="jemalloc-5.3.0"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.bz2
cd $pkg
./configure --prefix=$INSTALL_PATH >>$LOG 2>>$LOG
make -j $jn >>$LOG 2>>$LOG
make install >>$LOG 2>>$LOG
cd $CWD
rm -rf $pkg
echo ""

#================================================================
# Install Fermat
#================================================================
echo "Installing Fermat ..."
uo="$(uname -s)"
case "${uo}" in
    Linux*)     pkg="Ferl7";;
    Darwin*)    pkg="Ferm7a";;
esac
export pkg
rm -rf $pkg
tar xf $pkg.tar.gz
rm -rf $INSTALL_PATH/$pkg
mv $pkg $INSTALL_PATH/
cd "$INSTALL_PATH/bin"
ln -s -f ../$pkg/fer64 .
cd $CWD
echo ""

#================================================================
# Install Form
#================================================================
uo="$(uname -s)"
case "${uo}" in
    Linux*)     pkg="form-4.3.1-x86_64-linux";;
    Darwin*)    pkg="form-4.3.1-x86_64-osx";;
esac
export pkg
echo "Installing $pkg ..."
tar xf $pkg.tar.gz
cp -rf $pkg/form "$INSTALL_PATH/bin/"
cp -rf $pkg/tform "$INSTALL_PATH/bin/"
rm -rf $pkg
cd $CWD
echo ""

#================================================================
# Install FIRE - Modified version
#================================================================
export pkg="FIRE"
export LOG=$LOGS/$pkg.log
echo "Installing $pkg ..."
cd "$INSTALL_PATH"
tar xf "$CWD/$pkg.tar.gz"
cd $pkg
make -j $jn INSTALL_PATH="$INSTALL_PATH" >>$LOG 2>>$LOG
make clean >>$LOG 2>>$LOG
cd $CWD
echo ""

#================================================================
# Install HepLib
#================================================================
export pkg="HepLib"
echo "Installing $pkg ..."
rm -rf $pkg
tar xf $pkg.tar.gz
cd $pkg
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -DINC_PATH=/opt/homebrew/include -DLIB_PATH=/opt/homebrew/lib ..
make -j $jn
make install
cd $CWD
rm -rf $pkg
echo ""

echo ""
echo "Installation Completed!"
echo ""


