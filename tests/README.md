## Description
Here is developed the testing unit for all the lina functions that need numerical testing.

## Usage
For each function in the lina library named in the form of _lina_something()_, here is defined a folder named _something_. In each folder there are many tests, each one identified by a ti.txt file, for i=1,...,n.
Each test file is defined as follows: The first matrix/matrices are the inputs of the function (depending on the function, for example: lina_add() has two inputs A,B and one output C=A+B. A and B have to be the first two matrices in the test file), the last matrix/matrices are the output of the function and after there are input scalar values (ordered in the same order of the function under test) of the function represented as a 1x1 matrix.

For example, a scale test file, that is a test for the lina_scale() function is defined as follows:

[1 1 1,1 1 1,1 1 1]    
[2 2 2,2 2 2,2 2 2]   
[2]

Where the first matrix is the input, the second the output and the last is the scalar value.

By default, executing the test file will generate all the testing and provide the results on the stdout.
