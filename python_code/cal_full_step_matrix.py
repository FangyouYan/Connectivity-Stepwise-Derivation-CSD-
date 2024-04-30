# @author: code_king
# @date: 2024/4/9 pm：3:19
from typing import Dict, List, Set, Tuple

import numpy as np
from numpy import mat
from rdkit import Chem


class Algorithm_MS:
    @staticmethod
    def FullStepMatrixbyFloydWarshall(M_S_A):
        """
        Floyd-Warshall algorithm
        :param M_S_A: adjacency matrix
        :return: dis: MSF    path: Routing matrix
        """
        # dis = copy.deepcopy(M_S_A)  # adjacency matrix
        dis = M_S_A  # adjacency matrix
        dis = np.where(dis == 0, np.inf, dis)  # makes 0 infinite, otherwise 0 is the smallest by default
        n = np.shape(M_S_A)[0]
        for i in range(n):
            dis[i, i] = 0
        # path = np.zeros((n, n), dtype=int)  # Initialize the routing matrix so that it is all zeros
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if dis[i][k] + dis[k][j] < dis[i][j]:
                        # Find a shorter path through point k, and accept this shorter path length
                        dis[i][j] = dis[i][k] + dis[j][k]
                        # path[i][j] = k  # The routing matrix records paths
        # return dis, path  # dis:MSF      path: Routing matrix
        return dis

    @staticmethod
    def FullStepMatrixbyCSD(n_atom: int, Sa: [[int]], Ladj: Dict[int, List[int]], step_xy: Set[Tuple[int, int]]):
        """
        CSD algorithm
        :param n_atom: number of atoms
        :param Sa: adjacency matrix
        :param Ladj: atom connectivity
        :return: MSF
        """
        SF = Sa
        for m in range(1, n_atom - 1):
            if len(step_xy) == 0: break
            temp_step_xy = []
            for item in step_xy:
                r = item[0]
                for c in Ladj[item[1]]:
                    if SF[r - 1, c - 1] == 0 and r != c:
                        SF[r - 1, c - 1] = m + 1
                        temp_step_xy.append((r, c))
            step_xy = temp_step_xy
        return SF


class CommonUtils:
    @staticmethod
    def fastStrFromMol(ISI: [[str]]):
        si = []
        n_atom_start = 4
        n_atom = int(ISI[3][0:3])
        n_adj = int(ISI[3][3:6])
        for j in range(n_atom_start, n_atom_start + n_atom):
            si.append(ISI[j].split())
        si = np.array(si)
        return si, n_atom, n_adj, n_atom_start

    @staticmethod
    def getMSa(n_atom: int, n_atom_start: int, n_adj: [], ISI: [[str]]):
        Sa = mat(np.zeros((n_atom, n_atom)))
        for j in range(n_atom_start + n_atom, n_atom_start + n_atom + n_adj):
            i_s = min([int(ISI[j][0:3]) - 1, int(ISI[j][3:6]) - 1])
            i_g = max([int(ISI[j][0:3]) - 1, int(ISI[j][3:6]) - 1])
            Sa[i_s, i_g] = 1
            Sa[i_g, i_s] = 1
        return Sa

    @staticmethod
    def getMSaAndLadj(n_atom: int, n_atom_start: int, n_adj: [], ISI: [[str]]):
        Sa = mat(np.zeros((n_atom, n_atom)))
        Ladj = {}
        # get adjacency matrix and atom connectivity
        # define wsc ,wsr
        step_xy = set()
        for j in range(n_atom_start + n_atom, n_atom_start + n_atom + n_adj):
            i_s = int(ISI[j][0:3]) - 1
            i_g = int(ISI[j][3:6]) - 1
            Sa[i_s, i_g] = 1
            Sa[i_g, i_s] = 1
            # get Ladj
            Ladj.setdefault(i_s + 1, set()).add(i_g + 1)
            Ladj.setdefault(i_g + 1, set()).add(i_s + 1)
            # wsc,wsr
            step_xy.add((i_s + 1, i_g + 1))
            step_xy.add((i_g + 1, i_s + 1))
        return Sa, Ladj, step_xy

    @staticmethod
    def FastAdjacentLists(n_atom: int, Sa: [[int]]):
        Ladj = []
        for j in range(0, n_atom):
            Ladj_ = np.where(Sa[j, :] > 0)
            Ladj.append(list(1 + Ladj_[1]))
        return Ladj


class StructureInformation:
    @staticmethod
    def fast_mol(ISI):
        si, n_atom, n_adj, n_atom_start = CommonUtils.fastStrFromMol(ISI)
        Sa, Ladj, step_xy = CommonUtils.getMSaAndLadj(n_atom, n_atom_start, n_adj, ISI)
        # Ladj = CommonUtils.FastAdjacentLists(n_atom, Sa)
        SF = Algorithm_MS.FullStepMatrixbyCSD(n_atom, Sa, Ladj, step_xy)
        SI = {'n_atom': n_atom, 'SF': SF}
        return SI

    @staticmethod
    def fast_floyd_mol(ISI):
        si, n_atom, n_adj, n_atom_start = CommonUtils.fastStrFromMol(ISI)
        Sa = CommonUtils.getMSa(n_atom, n_atom_start, n_adj, ISI)
        SF = Algorithm_MS.FullStepMatrixbyFloydWarshall(Sa)  # Plus routing matrix
        SI = {'n_atom': n_atom, 'SF': SF}
        return SI


if __name__ == '__main__':
    with open('./data/ben.mol', 'r') as f:
        ISI = f.readlines()
        # get full step matrix from CSD and Floyd-Warshall
        SF_CSD = StructureInformation.fast_mol(ISI)
        SF_floyd = StructureInformation.fast_floyd_mol(ISI)
        # print MSF into console
        print(SF_CSD['SF'])
        print(SF_floyd['SF'])

    with open('./data/ben.mol', 'r') as f:
        content = f.read()
        # get full step matrix from RDKit
        mol = Chem.MolFromMolBlock(content, removeHs=False)
        SF_RDKit = Chem.GetDistanceMatrix(mol)
        # print MSF into console
        print(SF_RDKit)
