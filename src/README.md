# README
## INTRODUCTION
This repository stores the code for sIM-CoMeT C++ implementation developed by Jayghosh Rao for his Master's thesis work. The program calculates the steady state velocity and curve radius of the trajectory of a melting probe in a given Phase Change Material (PCM). With inputs of thermo-physical properties of the PCM, and the applied thermal flux/temperature, the program uses the CCMSOLVE algorithm to calculate the required variables. 

## Dependencies [To be manually installed]
1. suitesparse/umfpack 
	Ubuntu: `sudo apt-get install libsuitesparse-dev`
	Arch : `sudo pacman -S suitesparse`
2. gcc version 5 and above. 
	* Due to variations in the UMFPACK package headers in different OSes, there was a need to use the '\_\_has\_include' preprocessor command to first check the system include directory before including them. This command is only available in gcc 5 and above. 
3. gnuplot

## Other Packages Used
1. inih - https://github.com/benhoyt/inih
2. muParser

Note: these packages are integrated into the source code. There is no need to compile them from the upstream sources. 

## Usage

```
make clean
make 
./simcomet
```

### Notes: 
- Input values required are to be filled in the inputs.ini file. If any entry is missing,  default values will be used by the program.
- The outputs are stored in plain-text format in the directory 'outputs'. The format followed in the data is: <co-ordinates> <value> ; and as such may directly be used to construct plots using various tools.

contact email : jayghosh.rao@rwth-aachen.de
