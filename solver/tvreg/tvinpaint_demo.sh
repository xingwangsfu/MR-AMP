#! /bin/sh
# TV inpainting demo using tvrestore, Pascal Getreuer 2010
#
# This is a script that should run on most UNIX systems.  To 
# run the script, open a terminal and enter the command
#     chmod +x tvinpaint_demo.sh
# to make this script executable.  Then run the script with
#     ./tvinpaint_demo.sh


echo
echo '   --------------------------------'
echo '    TV-regularized inpainting demo '
echo '   --------------------------------'
echo
echo 'The input image signal-f.bmp is inpainted on the domain'
echo 'defined by signal-D.bmp.'
echo

echo './tvrestore D:signal-D.bmp lambda:1e3 signal-f.bmp inpainted.bmp'
./tvrestore D:signal-D.bmp lambda:1e3 signal-f.bmp inpainted.bmp

echo 
echo 'Please compare signal-f.bmp with inpainted.bmp.'
echo
