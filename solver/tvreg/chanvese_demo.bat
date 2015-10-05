@echo off
rem Chan-Vese two-phase segmentation demo, Pascal Getreuer 2010
rem
rem This is a script that should run on Windows systems.  To run
rem the script, open a command (MS-DOS) prompt, use cd to reach
rem the directory containing this file, and enter
rem     chanvese_demo


echo.
echo    -----------------------------
echo     Chan-Vese segmentation demo
echo    -----------------------------
echo.
echo The input image wrench.bmp is segmented using the Chan-Vese
echo active contours without edges method with parameter mu = 0.2.
echo.

echo chanvese mu:0.2 wrench.bmp chanvese.bmp
chanvese mu:0.2 wrench.bmp chanvese.bmp

echo.
echo Please compare wrench.bmp with chanvese.bmp.
echo.
