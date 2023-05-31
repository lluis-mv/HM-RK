#!/bin/bash


alias matlab='export MESA_LOADER_DRIVER_OVERRIDE=i965; /home/lmonforte/Matlab/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd ModifiedCamClay/VP

matlab "MakeFigureCases; quit" > output1 2>&1 
matlab "FigureSourceTerm; quit" > output1 2>&1 

matlab "NewExample1; quit" > output1 2>&1 & 
matlab "NewExample1a; quit" > output2 2>&1 & 
matlab "NewExample1b; quit" > output3 2>&1 & 
matlab "NewExample1c; quit" > output4 2>&1 & 
matlab "NewExample2; quit" > output5 2>&1 & 

cd ../../
cd VP-Footing

matlab "ExampleThree; quit" > output1 2>&1   &
matlab "ExampleSolutionCost; quit" > output1 2>&1   &
#matlab "ExampleParallel; quit" > output2 2>&1 &


