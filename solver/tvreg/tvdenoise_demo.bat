@echo off
rem TV denoising demo using tvrestore, Pascal Getreuer 2010
rem
rem This is a script that should run on Windows systems.  To run
rem the script, open a command (MS-DOS) prompt, use cd to reach
rem the directory containing this file, and enter
rem     tvdenoise_demo


echo.
echo    -----------------------------------
echo     TV-regularized noise removal demo
echo    -----------------------------------
echo.
echo The input image lighthouse.bmp is denoised using three
echo different denoising strengths: lambda = 5, 20, and 80.


echo.
echo Calling tvrestore with lambda = 5 (strong denoising)
echo Press any key to continue...
pause>nul

echo tvrestore lambda:5 lighthouse.bmp denoise1.bmp
tvrestore lambda:5 lighthouse.bmp denoise1.bmp


echo.
echo Calling tvrestore with lambda = 20 (moderate denoising)
echo Press any key to continue...
pause>nul

echo tvrestore lambda:20 lighthouse.bmp denoise2.bmp
tvrestore lambda:20 lighthouse.bmp denoise2.bmp


echo.
echo Calling tvrestore with lambda = 80 (light denoising)
echo Press any key to continue...
pause>nul

echo tvrestore lambda:80 lighthouse.bmp denoise3.bmp
tvrestore lambda:80 lighthouse.bmp denoise3.bmp


echo.
echo Please compare lighthouse.bmp with the denoising results
echo denoise1.bmp, denoise2.bmp, denoise3.bmp.
echo.
