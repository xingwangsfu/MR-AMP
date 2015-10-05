#! /bin/sh
# TV deconvolution demo using tvrestore, Pascal Getreuer 2010
#
# This is a script that should run on most UNIX systems.  To 
# run the script, open a terminal and enter the command
#     chmod +x tvdeconv_demo.sh
# to make this script executable.  Then run the script with
#     ./tvdeconv_demo.sh


echo
echo '   -----------------------------------'
echo '    TV-regularized deconvolution demo '
echo '   -----------------------------------'
echo
echo 'The input image blurry.bmp is deconvolved using a'
echo 'disk-shaped point spread function.'
echo

echo './tvrestore lambda:3e3 K:disk:3 blurry.bmp deconv.bmp'
./tvrestore lambda:3e3 K:disk:3 blurry.bmp deconv.bmp

echo 
echo 'Please compare blurry.bmp with deconv.bmp.'
echo
