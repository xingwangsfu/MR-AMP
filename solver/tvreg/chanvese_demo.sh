#! /bin/sh
# Chan-Vese two-phase segmentation demo, Pascal Getreuer 2010
#
# This is a script that should run on most UNIX systems.  To 
# run the script, open a terminal and enter the command
#     chmod +x chanvese_demo.sh
# to make this script executable.  Then run the script with
#     ./chanvese_demo.sh


echo
echo '   -----------------------------'
echo '    Chan-Vese segmentation demo '
echo '   -----------------------------'
echo
echo 'The input image wrench.bmp is segmented using the Chan-Vese'
echo 'active contours without edges method with parameter mu = 0.2.'
echo

echo './chanvese mu:0.2 wrench.bmp chanvese.bmp'
./chanvese mu:0.2 wrench.bmp chanvese.bmp

echo 
echo 'Please compare wrench.bmp with chanvese.bmp.'
echo
