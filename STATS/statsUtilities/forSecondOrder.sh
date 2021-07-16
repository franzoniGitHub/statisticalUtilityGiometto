#!/bin/sh
# Variables
    inputFileDIS="DIS"
    inputFileU="U"
    timeDir='*/'
# Prepare the environment
# -----------------------

echo " "
echo '-------------------'
echo "Looking for dissipation files..."
if [[ ! ( -e  ../stats/"$inputFileDIS"0 ) ]]; then
  cd ../../
for h in $timeDir
  do
  if [[ ( -e $h/$inputFileU ) ]]; then
     cd $h
     cp ../STATS/Umean .
     cd ..
  fi
done

  echo "Computing the dissipation field"
  cd STATS/computeDissipation
  rm -rf Make/linux*
  wmake > log.wmake
  cd ../../
  ./STATS/computeDissipation.exe > STATS/log.dissipation &
  echo "Now waiting for the results.. It may take a while (30 minutes or more).."
  wait
  cd STATS

echo "Deleting old dissipation files in stats"
rm stats/"$inputFileDIS"*

# Enter the time folders
cd ..
echo "Copying the dissipation fields into STATS/stats"
n=0
for h in $timeDir
  do
  if [[ ( $h != "0/" ) && ( -e $h/$inputFileU ) && ( -e $h/$inputFileDIS )]]; then
    echo ' ... copying data from '$h
    cd $h
     cp $inputFileDIS ../STATS/stats/"$inputFileDIS$n"
    cd ..
    ((n++))
  fi
done

cd STATS/stats
cat > "nSteps2.dat" << EOF
$n
EOF
cd ../statsUtilities
fi
echo '----------------------------------'
echo " "
