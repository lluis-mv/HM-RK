#!/bin/bash


alias matlab='export MESA_LOADER_DRIVER_OVERRIDE=i965; /home/lmonforte/Matlab/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd SFEM/

matlab "SFEM_Figure; quit" > output 2>&1  
matlab "SFEM_Figure_Quad; quit" > output 2>&1  

cd ../

cd SFEM-Example1

matlab "ExampleOne; quit" > output 2>&1 &


matlab "ExampleTwo; quit" > output2 2>&1 & 

cd ..


cd SFEM-Footing/

matlab "ExampleThree; quit" > output1 2>&1 & 

matlab "ExampleFour; ExampleFour2; quit" > output2 2>&1 &

matlab "ExampleFive; quit" > output3 2>&1  &

cd ..

cd SFEM-Biaxial-CASM/

matlab "ExampleThree; quit" > output1 2>&1 &

cd ..


cd QFEM-Example1

matlab "ExampleOne; quit" > output 2>&1 &


matlab "ExampleTwo; quit" > output2 2>&1 & 

cd ..


cd QFEM-Footing/

matlab "ExampleThree; quit" > output1 2>&1 & 

matlab "ExampleFour; ExampleFour2; quit" > output2 2>&1 &


cd ..



