#!/bin/bash


alias matlab='/usr/local/MATLAB/R2020a/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd SFEM/

matlab "SFEM_Figure; quit" > output 2>&1 &

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

