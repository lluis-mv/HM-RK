#!/bin/bash


alias matlab='/home/lluis/matlab2021a/bin/matlab -nodesktop -nosplash -noFigureWindows -r'

cd Example3/


matlab "ExampleThree; quit" > output3 2>&1 &




