#!/bin/bash
#
# See the documentation in README.md. When you're ready, just execute ./run.sh
# and the output will be placed in a new directory named 'build'.
#

# these tools will be used
gcc=`which gcc`
make=`which make`
patch=`which patch`
head=`which head`
tar=`which tar`
gpp=`which g++`
makeParallel=-j8

# get the location of this script
scriptDir=$(cd $(dirname $0) && pwd)

# create the output directory
outDir=$scriptDir/build
mkdir -p $outDir

# define some useful variables
filesDir=$scriptDir/files
pyversion=2.7
pythonPackageName=Python-2.7.3
numpyPackageName=numpy-1.6.2

pythonSourceDir=$outDir/$pythonPackageName
numpySourceDir=$outDir/$numpyPackageName
installDir=$outDir/install
python=$installDir/bin/python
freezeOutputDir=$outDir/freezeOutput

# remove python related environment variables
unset PYTHONHOME
unset PYTHONPATH
unset PYTHONSTARTUP


############################
# define the script routines

extractAndBuildPython()
{
  cd $outDir
  $tar -zxf $filesDir/$pythonPackageName.tgz
  cd $pythonSourceDir

  # patch python import
  $patch -p0 < $filesDir/patch-python-import.diff

  # setup builtin modules
  cp $filesDir/Setup.local ./Modules/

  ./configure install --disable-shared --enable-unicode=ucs4 --prefix $installDir || exit
  $make $makeParallel || exit
  $make install || exit
}


extractAndBuildNumpy()
{
  cd $outDir
  $tar -zxf $filesDir/$numpyPackageName.tar.gz
  cd $numpyPackageName
  $python ./setup.py install --prefix $installDir || exit
}


runFreezeTool()
{
  targetScript=$filesDir/freezeInputScript.py
  freezeScript=$pythonSourceDir/Tools/freeze/freeze.py
  ignores="-X codecs -X copy -X distutils -X encodings -X locale -X macpath -X ntpath -X os2emxpath -X popen2 -X pydoc  -X readline -X _warnings"

  export PYTHONHOME=$installDir
  export PYTHONPATH=$installDir/lib/python$pyversion/site-packages

  $python -S $freezeScript -o $freezeOutputDir -p $pythonSourceDir $ignores $targetScript || exit

  unset PYTHONHOME
  unset PYTHONPATH

  # convert frozen.c to frozen.h by removing the "int main()" function at the end
  cd $freezeOutputDir
  numberOfLines=$(cat frozen.c | wc -l)
  numberOfLinesToKeep=$(expr $numberOfLines - 9)
  $head -n $numberOfLinesToKeep frozen.c > frozen.h
}

buildloader()
{
  cd $outDir
  cp -r $filesDir/loader ./
  cd loader

  pythonInclude=$installDir/include/python$pyversion
  pythonLibrary=$installDir/lib/libpython$pyversion.a
  numpyBuildDir=$numpySourceDir/build

  loaderFiles="linker_generator.c linux_functions.c cpp-wrappers.o"
  boincIncludes="-I$BOINC_DIR/api -I$BOINC_DIR/lib -I$BOINC_DIR/zip"

  # glob numpy .o files and frozen .c source files
  numpyObjects=$(find $numpyBuildDir -name \*.o)
  frozenSources=$(find $freezeOutputDir -name M_\*.c)

  ln -s `$gpp -print-file-name=libstdc++.a`

  $gpp -c cpp-wrappers.cpp $boincIncludes -o cpp-wrappers.o

  # osx and linux take slightly different linking flags
  if [ "$(uname)" == "Darwin" ]; then
    linkFlags=""
    requiredBlasLibrary="-framework Accelerate"
  else
    linkFlags="-Xlinker -export-dynamic"
    requiredBlasLibrary="/usr/lib/atlas-base/atlas/liblapack.a /usr/lib/atlas-base/atlas/libblas.a libstdc++.a"
    ARCH=$(uname -m)
    if [ ${ARCH} == 'x86_64' ]; then
        requiredBlasLibrary="$requiredBlasLibrary /usr/lib/gcc/x86_64-linux-gnu/4.4/libgfortran.a"
        loaderFiles="$loaderFiles  res/x86_64/libboinc_api.a res/x86_64/libboinc_zip.a res/x86_64/libboinc.a"
    else
        requiredBlasLibrary="$requiredBlasLibrary /usr/lib/gcc/i486-linux-gnu/4.4/libgfortran.a"
        loaderFiles="$loaderFiles  res/i486/libboinc_api.a res/i486/libboinc_zip.a res/i486/libboinc.a"
    fi
  fi

  requiredLibs="-lpthread -lm -ldl -lutil -static-libgcc"

  # this build includes the numpy object files and frozen source files
  $gcc $linkFlags -o linker_gen $loaderFiles $boincIncludes -I$pythonInclude -I$freezeOutputDir \
        $frozenSources \
        $numpyObjects \
        $pythonLibrary \
        $requiredBlasLibrary \
        $requiredLibs

}

buildPythonModule()
{
  cd $outDir
  cp -r $filesDir/python_module ./
  cd python_module
  cp $outDir/install/lib/libpython2.7.a .
#  nuitka --module --standalone --python-version=2.7 --nofreeze-stdlib --recurse-not-to=numpy --recurse-not-to=loader --recurse-all --show-progress --show-scons linker_generator.py | tee nuitka.log
  ln -s `$gpp -print-file-name=libstdc++.a`
#  linkingStep=$(cat nuitka.log | tail -n 3 | head -n 1)
#  linkingStep=$(echo $linkingStep | sed 's/-lpython2.7/libpython2.7.a/')
#  linkingStep= "$linkingStep -static-libgcc -L."
#  echo $linkingStep
#  $linkingStep
}

#makeDeploy()
#{
#  mkdir deploy
#  cd deploy
#  cp $outDir/loader/linker_gen .
#  cp $outDir/python_module/linker_generator.dist/linker_generator.so .
#  zip -r linkergen_res.zip linker_generator.so
#}

##########################
# call the script routines

set -x
extractAndBuildPython
extractAndBuildNumpy
runFreezeTool
buildloader
#buildPythonModule
