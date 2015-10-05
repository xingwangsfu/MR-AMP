#! /bin/sh
# TV denoising demo using tvrestore, Pascal Getreuer 2010
#
# This is a script that should run on most UNIX systems.  To 
# run the script, open a terminal and enter the command
#     chmod +x tvdenoise_demo.sh
# to make this script executable.  Then run the script with
#     ./tvdenoise_demo.sh


echo
echo '   -----------------------------------'
echo '    TV-regularized noise removal demo '
echo '   -----------------------------------'
echo
echo 'The input image lighthouse.bmp is denoised using three'
echo 'different denoising strengths: lambda = 5, 20, and 80.'


echo
echo 'Calling tvrestore with lambda = 5 (strong denoising)'
read -p 'Press any key to continue...'

echo './tvrestore lambda:5 lighthouse.bmp denoise1.bmp'
./tvrestore lambda:5 lighthouse.bmp result1.bmp


echo
echo 'Calling tvrestore with lambda = 20 (moderate denoising)'
read -p 'Press any key to continue...'

echo './tvrestore lambda:20 lighthouse.bmp denoise2.bmp'
./tvrestore lambda:20 lighthouse.bmp result2.bmp


echo
echo 'Calling tvrestore with lambda = 80 (light denoising)'
read -p 'Press any key to continue...'

echo './tvrestore lambda:80 lighthouse.bmp denoise3.bmp'
./tvrestore lambda:80 lighthouse.bmp result3.bmp


echo 
echo 'Please compare lighthouse.bmp with the denoising results'
echo 'denoise1.bmp, denoise2.bmp, denoise3.bmp.'
echo
