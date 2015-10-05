@echo off
rem TV deconvolution demo using tvrestore, Pascal Getreuer 2010
rem
rem This is a script that should run on Windows systems.  To run
rem the script, open a command (MS-DOS) prompt, use cd to reach
rem the directory containing this file, and enter
rem     tvdeconv_demo

echo.
echo    -----------------------------------
echo     TV-regularized deconvolution demo
echo    -----------------------------------
echo.
echo The input image blurry.bmp is deconvolved using a
echo disk-shaped point spread function.
echo.

echo tvrestore lambda:3e3 K:disk:3 blurry.bmp deconv.bmp
tvrestore lambda:3e3 K:disk:3 blurry.bmp deconv.bmp

echo.
echo Please compare blurry.bmp with deconv.bmp.
echo.
