#!/bin/bash

proc=`ls dambreak_p*.h5 | wc -l`

# Create capstone xmf file 
echo "Gather xmf data ...."
rm -f dambreak_p_all*
$PROTEUS/proteusModule/scripts/gatherArchives.py -s $proc -f dambreak_p  > log

# Create directory
mkdir -p solution
cd solution

# Extract pressure and heights
echo "Extract Heights ..."
$PROTEUS_MPRANS/scripts/dambreakExtractHeightAndWidth.py -f ../dambreak_p_all*   >> ../log

# Repackage full xmf in lean xmf and vtk  
end=`cat height.txt | wc -l`
$PROTEUS_MPRANS/scripts/extractSolution.py -f ../dambreak_p -s 0 -e $end -i 10 -n $proc  

rm -rvf sol.*.h5
