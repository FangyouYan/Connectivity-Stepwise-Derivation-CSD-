//
// Created by code_king on 2024/4/23.
//
#include <vector>
#include <string>
#include "../common_utils/MsUtils.h"


int main() {
    std::string filename = std::filesystem::current_path().parent_path().string() +"/mol_data/benzene.mol";
    std::vector<std::string> mol_content=readMolFile(filename);
    std::vector<std::vector<float>> MSF= getMSFbyCSD(mol_content);
    // print the MSF into console
    printMSF(MSF);
    return 0;
}

