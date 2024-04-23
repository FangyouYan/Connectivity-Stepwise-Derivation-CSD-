//
// Created by king on 2024/3/2.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <sstream>
#include <iterator>

#include <algorithm>
#include <tuple>
#include <ostream>

namespace fs = std::filesystem;
using namespace std;

// define the matrix
using Matrix = std::vector<std::vector<float>>;

std::tuple<std::vector<std::string>, int, int, int> betterFastStrFromMol(const std::vector<std::string> &ISI);

std::tuple<Matrix>
FastStepBondFromMol(int n_atom, int n_atom_start, int n_adj, const std::vector<std::string> &ISI);

// generate the MSA faster
std::tuple<std::vector<std::vector<int>>>
FastAdjacentLists(int n_atom, std::vector<std::vector<float>> Sa);

//  FullStepMatrix implements
std::vector<std::vector<float>>
FastFullStepMatrix(int n_atom, std::vector<std::vector<float>> &Sa, std::vector<std::vector<int>> &Ladj);

vector<vector<float>> floyd(vector<vector<float>> MSa);

std::vector<std::string> readMolFile(const std::string& filename);

std::vector<std::vector<float>> getMSFbyCSD(std::vector<std::string> mol_content);

std::vector<std::vector<float>> getMSFbyFloydWarshall(const vector<std::string> &mol_content);


void printMSF(vector<std::vector<float>> &MSF);
