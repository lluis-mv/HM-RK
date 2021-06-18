#!/bin/bash


alias matlab='/home/lluis/matlab2021a/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd Example1/
matlab "ExampleOne; quit" > output 2>&1 &

matlab "InfluenceStab; quit" > output2 2>&1 &

cd ..

cd Example1Bis/


matlab "ExampleOneBis; quit" > output 2>&1 &

cd ..

cd Example3/

matlab "CompareElements; quit" > output 2>&1

matlab "CompareElementsElastic; quit" > output1 2>&1 &

matlab "PreExampleThree; quit" > output2 2>&1

matlab "ExampleThree; quit" > output3 2>&1 &

matlab "ExampleThreeElastic; quit" > output4 2>&1 &



