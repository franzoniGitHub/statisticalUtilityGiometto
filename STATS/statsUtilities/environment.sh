#!/bin/sh
# Variables for the name of the input/output files
    inputFileb='b'
    inputFileU='U'
    inputFileNUT='nut'
    timeDir='*/'
# Prepare the environment
# -----------------------
echo " "
echo '-------------------'
echo "Creating the input and output directories"
cd ../
mkdir stats
mkdir -p output/autocorrelation/alongX
mkdir output/autocorrelation/alongY
mkdir output/firstOrder
mkdir output/secondOrder
cd ../

# Enter the time folders
n=0
for h in $timeDir
  do
  if [[ ( $h != "0/" ) && ( -e $h/$inputFileb ) && ( -e $h/$inputFileU ) && ( -e $h/$inputFileNUT ) ]]; then
    echo 'Copying data from '$h
    cd $h
     cp $inputFileb ../STATS/stats/"$inputFileb$n"
     cp $inputFileU ../STATS/stats/"$inputFileU$n"
     cp $inputFileNUT ../STATS/stats/"$inputFileNUT$n"
    cd ..
    ((n++))
  fi
done

cd STATS/stats
cat > "nSteps.dat" << EOF
$n
EOF
cd ../statsUtilities
echo '-------------------'
echo " "
