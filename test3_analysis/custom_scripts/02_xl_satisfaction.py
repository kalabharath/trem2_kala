import sys, os, glob
import csv
import collections
from itertools import combinations
from scipy.spatial.distance import cdist
from multiprocessing import  Pool

def read_xls(xls_path):
    tdata = []
    with open(xls_path) as csvfile:
        data = csv.reader(csvfile)
        for row in data:
            try:
                row[1] = int(row[1])
                row[3] = int(row[3])
                tdata.append(row)
            except ValueError:
                pass
    return tdata


def read_pdb(pdb_file):
    coor_dict = collections.defaultdict(list)
    additional_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
                         'T', 'U', 'V',
                         'W', 'X', 'Y', 'Z', ]

    with open(pdb_file) as fin:
        lines = fin.readlines()
    for line in lines:
        if line[21:22] == ' ':
            coor_dict['trem2'].append(line)
        elif line[21:22] in additional_chains:
            coor_dict[line[21:22]].append(line)
        else:
            pass

    new_coor_dict = collections.defaultdict(list)
    chains = coor_dict.keys()
    for chain in chains:
        tcoors = coor_dict[chain]
        last_line = tcoors[-1]
        last_line = last_line.split()
        if chain == 'trem2':
            no_of_resi = int(last_line[4])
        else:
            no_of_resi = int(last_line[5])
        temp_coors = [[999.999, 999.999, 999.999]] * (no_of_resi + 1)
        for t in tcoors:
            t = t.split()
            if chain == 'trem2':
                temp_coors[int(t[4])] = [float(t[5]), float(t[6]), float(t[7])]
            else:
                temp_coors[int(t[5])] = [float(t[6]), float(t[7]), float(t[8])]
        new_coor_dict[chain] = temp_coors
    return new_coor_dict

def get_dist(coor1, coor2):
    pass
    return True


def check_xls(tpdb):
    coors = read_pdb(tpdb)
    xls_path = "../../data/xl/bs3_nshsf_mixed.csv"
    xls = read_xls(xls_path)
    percent = 0.0
    total_xls = float(len(xls))
    offset = 0
    xl_cutoff = 25.0  # Angstrom
    pdbids = []
    additional_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
                         'T', 'U', 'V',
                         'W', 'X', 'Y', 'Z', ]
    beta_fibril_combi = list(combinations(additional_chains, 3))
    percent_satisfaction = collections.defaultdict(list)
    for fibrl in beta_fibril_combi:
        xl_count = 0
        xl_found = []
        for chain in fibrl:
            for xl in xls:
                coor1 = coors['trem2'][xl[-1] - offset]
                coor2 = coors[chain][xl[1]-offset]
                dist = cdist([coor1], [coor2], 'euclidean')
                # print xl, dist
                if dist < xl_cutoff:
                    if xl in xl_found:
                        pass
                    else:
                        xl_count += 1
                        xl_found.append(xl)
                else:
                    pass
        percent_satisfaction[xl_count / total_xls].append(fibrl)

    tnum = percent_satisfaction.keys()
    tnum.sort()
    tnum.reverse()
    if tnum[0] > 0.6 and len(percent_satisfaction[tnum[0]]) > 20:
        print tpdb, percent_satisfaction[tnum[0]], tnum[0]
        pdbids.append(tpdb)

    return True


def run():

    tpdb = "./h1_run32_1118.pdb"
    pdbs = glob.glob("../sampled_pdbs/h1*.pdb")

    p = Pool(35)
    p.map(check_xls, pdbs)

    """

    for tpdb in pdbs:

        chains, percent = check_xls(tpdb, xls)

        if percent > 0.5:
            print tpdb, percent, chains
        else:
            print tpdb, percent
            pass
    return True
    """

if __name__ == '__main__':
    run()
