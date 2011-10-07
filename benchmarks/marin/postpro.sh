#!/bin/bash

proc=`ls marin_p*.h5 | wc -l`

# Create capstone xmf file 
echo "Gather xmf data ...."
rm -f marin_p_all*
$PROTEUS/proteusModule/scripts/gatherArchives.py -s $proc -f marin_p  > log

# Create directory
mkdir -p solution
cd solution

# Extract pressure and heights
echo "Extract Pressures ... "
$PROTEUS_MPRANS/scripts/marinExtractPressure.py -f ../marin_p_all* >> ../log
echo "Extract Heights ..."
$PROTEUS_MPRANS/scripts/marinExtractHeight.py -f ../marin_p_all*   >> ../log

# Repackage full xmf in lean xmf and vtk  
proc=`ls ../marin_p*.h5 | wc -l`
end=`cat pressure.txt | wc -l`
$PROTEUS_MPRANS/scripts/extractSolution.py -f ../marin_p -s 0 -e $end -i 10 -n $proc  

rm -rvf sol.*.h5
