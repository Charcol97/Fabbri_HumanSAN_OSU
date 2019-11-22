%%%%%%%%%%%%%   run this under a linux computer %%%%%%%%%%%%%%%

%% unzip the test.zip file into the same directory as the c code.

%% Compile c code to executable and outputs with a.out 
gcc ExampleCcode.c -mcmodel=large -lm

%% run the executable to produce *.vtk files for outputs with each time step
./a.out

%% using free opensoftware paraview to open the *.vtk outputs (https://www.paraview.org/download/) and tutorial for paraview is attached here


