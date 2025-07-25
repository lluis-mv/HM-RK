#!/bin/bash


alias matlab='export MESA_LOADER_DRIVER_OVERRIDE=i965; /usr/local/MATLAB/R2024a/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd NewExample1/

matlab "MakeSnapShotsOfTheSolution; quit" > output 2>&1 &


matlab "PlotTheMeshes; quit" > output 2>&1 &

matlab "Example0; quit" > output 2>&1 &

matlab "ExampleOneUgly; quit" > output2 2>&1 &

cd ..


cd NewExample3/

matlab "CompareElements; quit" > output 2>&1 &

matlab "CompareElementsElastic; quit" > output1 2>&1 &

matlab "CompareElementsElasticDrained; quit" > output1 2>&1 

matlab "PreExampleThree; quit" > output2 2>&1 &

matlab "ExampleThree; quit" > output3 2>&1 &

matlab "ExampleThreeElastic; quit" > output4 2>&1 &

matlab "ExampleThreeElasticDrained; quit" > output4 2>&1 &


