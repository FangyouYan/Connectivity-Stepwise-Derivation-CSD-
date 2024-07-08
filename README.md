## Project introduction
This project aims to provide two methods for calculating the Full Step Matrix (distance matrix) of a molecule. It uses the Floyd-Warshall algorithm and the CSD algorithm to calculate the molecular structure information and compares it with the distance matrix provided by the RDKit library.

## File structure（python_code）
- 'cal_full_step_matrix.py' : The main program file, which contains the main logic and algorithm implementation for computing MSF.
- 'data/ben.mol' : The "ben.mol" file of the sample molecule, used to test and demonstrate the algorithm's functionality.

## dependence
- Python (3.8.10)
- numpy (1.23.1)
- rdkit (2023.9.5)

## Installation
1. Make sure you have Python 3.8 or later installed.
2. Use pip to install the required libraries:
```sh
pip install numpy rdkit
```

## Usage
1. Place the sample .mol file in the data directory.
2. Run the cal_full_step_matrix.py file:
```sh
python cal_full_step_matrix.py
```

if you are using cpp, referring the "cpp_code" directory.

## File structure（cpp_code）
- 'main/getMSFbyCSD.py' : The program file, which contains the main logic and algorithm implementation for computing MSF using CSD.
- 'main/getMSFbyFloydWarshall.py' : The program file, which contains the main logic and algorithm implementation for computing MSF using Floyd-Warshall.
- 'mol_data/benzene.mol' : The "benzene.mol" file of the sample molecule, used to test and demonstrate the algorithm's functionality.

## dependence
- C++ (standard 17)
- Make sure you have the necessary C++ compilers and build toolchains installed on your system, such as CMake and make.

## Usage(windows)
1.Complete the pre-work
```shell
# enter the cpp code directory
cd cpp_code
# creat build directory
mkdir build
# enter the build directory
cd build
```

2.config the project
```shell
cmake -G "MinGW Makefiles" ..
```

3.build the project
```shell
# make sure you have installed this environment.
mingw32-make
# if you installed make,you can also use blow command.
# make
```
4.run the executable program
```shell
getMSFbyCSD.exe
getMSFbyFloydWarshall.exe
```