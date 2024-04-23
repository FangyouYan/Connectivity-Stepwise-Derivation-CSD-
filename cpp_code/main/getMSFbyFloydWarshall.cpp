//
// Created by king on 2024/4/15.
//
#include <iostream>
#include <vector>
#include "../common_utils/MsUtils.h"


int main() {
    std::string filename = std::filesystem::current_path().parent_path().string() +"/mol_data/benzene.mol";
    std::vector<std::string> mol_content=readMolFile(filename);
    std::vector<std::vector<float>> MSF=getMSFbyFloydWarshall(mol_content);
    // print the MSF into console
    printMSF(MSF);
    return 0;
}

