# @author: code_king
# @date: 2024/4/9 pm：3:19
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
        dis = M_S_A  # 邻接矩阵
        dis = np.where(dis == 0, np.inf, dis) # makes 0 infinite, otherwise 0 is the smallest by default
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
    def FullStepMatrixbyCSD(n_atom:int, Sa:[[int]], Ladj:[[int]]):
        """
        CSD algorithm
        :param n_atom: number of atoms
        :param Sa: adjacency matrix
        :param Ladj: atom connectivity
        :return: MSF
        """
        SF = Sa
        w_ms = np.where(SF == 1)
        w_ms_r = w_ms[0]
        w_ms_c = w_ms[1]
        for m in range(1, n_atom - 1):
            if len(w_ms_r) == 0:
                break
            w_ms_r_ = []
            w_ms_c_ = []
            # w_ms_r is constant level, looking for atoms with step size n, up to n-1
            for k, i in enumerate(w_ms_r):
                w_msi = Ladj[w_ms_c[k]]
                # w_msi Adjacent atom. It's a constant, maximum n-1
                for j in w_msi:
                    if SF[i, j - 1] == 0 and i != j - 1:
                        SF[i, j - 1] = m + 1
                        # SF[j - 1, i] = m + 1
                        w_ms_r_.append(i)
                        w_ms_c_.append(j - 1)
                        # w_ms_r_.append(j - 1)
                        # w_ms_c_.append(i)
            w_ms_r = w_ms_r_
            w_ms_c = w_ms_c_
        return SF

class CommonUtils:
    @staticmethod
    def fastStrFromMol(ISI:[[str]]):
        si = []
        n_atom_start = 4
        n_atom = int(ISI[3][0:3])
        n_adj = int(ISI[3][3:6])
        for j in range(n_atom_start, n_atom_start + n_atom):
            si.append(ISI[j].split())
        si = np.array(si)
        return si, n_atom, n_adj, n_atom_start

    @staticmethod
    def FastStepBondFromMol(n_atom:int, n_atom_start:int, n_adj:[], ISI:[[str]]):
        Sa = mat(np.zeros((n_atom, n_atom)))
        for j in range(n_atom_start + n_atom, n_atom_start + n_atom + n_adj):
            i_s = min([int(ISI[j][0:3]) - 1, int(ISI[j][3:6]) - 1])
            i_g = max([int(ISI[j][0:3]) - 1, int(ISI[j][3:6]) - 1])
            Sa[i_s, i_g] = 1
            Sa[i_g,i_s] = 1
        return Sa

    @staticmethod
    def FastAdjacentLists(n_atom:int, Sa:[[int]]):
        Ladj = []
        for j in range(0, n_atom):
            Ladj_ = np.where(Sa[j, :] > 0)
            Ladj.append(list(1 + Ladj_[1]))
        return Ladj

class StructureInformation:
    @staticmethod
    def fast_mol(ISI):
        si, n_atom, n_adj, n_atom_start = CommonUtils.fastStrFromMol(ISI)
        Sa = CommonUtils.FastStepBondFromMol(n_atom, n_atom_start, n_adj, ISI)
        Ladj = CommonUtils.FastAdjacentLists(n_atom, Sa)
        SF = Algorithm_MS.FullStepMatrixbyCSD(n_atom, Sa, Ladj)
        SI = {'n_atom': n_atom, 'SF': SF}
        return SI

    @staticmethod
    def fast_floyd_mol(ISI):
        si, n_atom, n_adj, n_atom_start = CommonUtils.fastStrFromMol(ISI)
        Sa = CommonUtils.FastStepBondFromMol(n_atom, n_atom_start, n_adj, ISI)
        SF = Algorithm_MS.FullStepMatrixbyFloydWarshall(Sa)  # Plus routing matrix
        # SI = {'n_atom': n_atom, 'SF': SF[0]}
        SI = {'n_atom': n_atom, 'SF': SF}
        return SI

if __name__ == '__main__':
    with open('./data/ben.mol', 'r') as f:
        ISI =f.readlines()
        SF_CSD=StructureInformation.fast_mol(ISI)
        SF_floyd=StructureInformation.fast_floyd_mol(ISI)
        # print MSF into console
        print(SF_CSD['SF'])
        print(SF_floyd['SF'])

    with open('./data/ben.mol', 'r') as f:
        content=f.read()
        mol = Chem.MolFromMolBlock(content, removeHs=False)  # removeHs是控制去除H原子
        SF_RDKit = Chem.GetDistanceMatrix(mol)
        print(SF_RDKit)