//
// Created by code_king on 2024/3/2.
//

#include "MsUtils.h"
#include <cfloat>


std::tuple<int, int, int> betterFastStrFromMol(const std::vector<std::string> &ISI) {
    int n_atom_start = 4;
    int n_atom = std::stoi(ISI[3].substr(0, 3));
    int n_adj = std::stoi(ISI[3].substr(3, 3));
    return {n_atom, n_adj, n_atom_start};
}


// get adjacency matrix
std::tuple<Matrix>
FastStepBondFromMol(int n_atom, int n_atom_start, int n_adj, const std::vector<std::string> &ISI) {
    // init Sa matrix,it is adjacency matrix
    Matrix Sa(n_atom, std::vector<float>(n_atom, 0.0));
    for (int j = n_atom_start + n_atom; j < n_atom_start + n_atom + n_adj; ++j) {
        const std::string &row = ISI[j];
        int num1 = std::stoi(row.substr(0, 3));
        int num2 = std::stoi(row.substr(3, 3));
        int i_s = std::min(num1 - 1, num2 - 1);
        int i_g = std::max(num1 - 1, num2 - 1);
        Sa[i_s][i_g] = 1.0;
        Sa[i_g][i_s]= 1.0;
    }
    return {Sa};
}

// get atom connectivity information
std::tuple<std::vector<std::vector<int>>,std::vector<int>,std::vector<int>>
FastAdjacentLists(int n_atom,  std::vector<std::vector<float>> Sa) {
    std::vector<std::vector<int>> Ladj;
    // ms row and column
    std::vector<int> w_ms_r;
    std::vector<int> w_ms_c;
    for (int j = 0; j < n_atom; ++j) {
        std::vector<int> Ladj_;
        for (size_t i = 0; i < Sa[j].size(); ++i) {
            if (Sa[j][i] > 0) {
                Ladj_.push_back(i + 1);
                w_ms_r.push_back(j);
                w_ms_c.push_back(i);
            }
        }
        Ladj.push_back(Ladj_);
    }
    // w_ms_r and w_ms_c are used for traversing the first step
    return std::make_tuple(Ladj,w_ms_r,w_ms_c);
}

//  full step matrix implements
std::vector<std::vector<float>>
FastFullStepMatrix(int n_atom, std::vector<std::vector<float>> &Sa, std::vector<std::vector<int>> &Ladj,std::vector<int> w_ms_r,std::vector<int> w_ms_c) {
    std::vector<std::vector<float>> SF = Sa;

    // Iterate m from 1 to n_atom-2
    for (int m = 1; m < n_atom - 1; ++m) {
        if (w_ms_r.empty()) {
            break;
        }
        std::vector<int> w_ms_r_;
        std::vector<int> w_ms_c_;
        std::vector<int> w_msi;

        for (int k = 0; k < w_ms_r.size(); ++k) {
            w_msi = Ladj[w_ms_c[k]];
            for (int j: w_msi) {
                if (SF[w_ms_r[k]][j - 1] == 0 && w_ms_r[k] != j - 1) {
                    SF[w_ms_r[k]][j - 1] = m + 1;
                    w_ms_r_.push_back(w_ms_r[k]);
                    w_ms_c_.push_back(j - 1);
                }
            }
        }
        w_ms_r = w_ms_r_;
        w_ms_c = w_ms_c_;
    }
    return SF;
}

std::vector<std::string> readMolFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::string> lines;
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return lines;
    }
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();
    return lines;
}

vector <vector<float>> floyd(vector <vector<float>> MSa) {
    int n = MSa.size();
    vector <vector<float>> dis = MSa;
    // Use transform and lambda expressions to replace all zeros with FLT_MAX
    std::transform(dis.begin(), dis.end(), dis.begin(),
                   [](std::vector<float> &row) {
                       std::transform(row.begin(), row.end(), row.begin(),
                                      [](float val) { return val == 0.0f ? FLT_MAX : val; });
                       return row;
                   });
    for (int i = 0; i < n; i++) {
        dis[i][i] = 0;
    }
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dis[i][k] + dis[k][j] < dis[i][j]) {
                    dis[i][j] = dis[i][k] + dis[k][j];
                    // path[i][j] = k;
                }
            }
        }
    }
    return dis;
}

std::vector<std::vector<float>> getMSFbyCSD(std::vector<std::string> mol_content){
    std::vector<std::vector<float>> SF;
    auto [n_atom, n_adj, n_atom_start] = betterFastStrFromMol(mol_content);
    auto [Sa] = FastStepBondFromMol(n_atom, n_atom_start, n_adj, mol_content);
    auto [Ladj,w_ms_r,w_ms_c] = FastAdjacentLists(n_atom,Sa);
    SF = FastFullStepMatrix(n_atom, Sa, Ladj,w_ms_r,w_ms_c);
    return SF;
}

std::vector<std::vector<float>> getMSFbyFloydWarshall(const vector<std::string> &mol_content) {
    vector<vector<float>> SF;
    auto [n_atom, n_adj, n_atom_start] = betterFastStrFromMol(mol_content);
    auto [Sa] = FastStepBondFromMol(n_atom, n_atom_start, n_adj, mol_content);
    SF = floyd(Sa);
    return SF;
}
void printMSF(vector<std::vector<float>> &MSF) {
    for (const auto &item: MSF) {
        for (const auto &e: item) {
            cout << e << " ";
        }
        cout<<"\n";
    }
}
