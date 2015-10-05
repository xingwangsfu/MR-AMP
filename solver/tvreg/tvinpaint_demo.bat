@echo off
rem TV inpainting demo using tvrestore, Pascal Getreuer 2010
rem
rem This is a script that should run on Windows systems.  To run
rem the script, open a command (MS-DOS) prompt, use cd to reach
rem the directory containing this file, and enter
rem     tvinpaint_demo


echo.
echo    --------------------------------
echo     TV-regularized inpainting demo 
echo    --------------------------------
echo.
echo The input image signal-f.bmp is inpainted on the domain
echo defined by signal-D.bmp.
echo.

echo tvrestore D:signal-D.bmp lambda:1e3 signal-f.bmp inpainted.bmp
tvrestore D:signal-D.bmp lambda:1e3 signal-f.bmp inpainted.bmp

echo.
echo Please compare signal-f.bmp with inpainted.bmp.
echo.
