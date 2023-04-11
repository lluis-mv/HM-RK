#!/bin/bash


alias matlab='export MESA_LOADER_DRIVER_OVERRIDE=i965; /home/lmonforte/Matlab/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd ModifiedCamClay/VP

matlab "FirstExample; quit" > output2 2>&1 & 
matlab "FirstExampleSecondPart; quit" > output2 2>&1 & 
matlab "ThirdExample; quit" > output3 2>&1 & 

cd ../../
cd VP-Footing

#matlab "ExampleThree; quit" > output1 2>&1   &
#matlab "ExampleParallel; quit" > output2 2>&1 &


