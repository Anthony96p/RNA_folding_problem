#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import time
import math
import argparse
import numpy as np

from TP_RNA_mainFonctions import calcul_dist
from TP_RNA_mainFonctions import pair_res_format


def linear_interpolation(directory, pair_res, dist, Gibbs):
    if dist <= 20:

        train_pair_res = np.genfromtxt(directory + '/' + pair_res + '.txt')

        x1, x2, y1, y2 = math.floor(dist), math.ceil(dist), train_pair_res[math.floor(dist), 1], train_pair_res[
            math.ceil(dist), 1]
        # print(x1,x2,y1,y2)

        E_scores = y1 + (dist - x1) * ((y2 - y1) / (x2 - x1))  # x2-x1=1
        # print (E_scores)

        Gibbs = Gibbs + E_scores
        # print (Gibbs)

    return Gibbs


def main_script3(train_dir, test_file):
    start = time.time()

    Gibbs_free_energy = 0

    # The script opens and browses the file 2 times in order to browse all possible peers
    with open(test_file, 'r') as file1:
        for line1 in file1:
            if (re.match('^ATOM', line1)) and line1[12:16] == " C3'":
                # --------------------------------------------------------------------------- #
                with open(test_file, 'r') as file2:
                    for line2 in file2:
                        if (re.match('^ATOM', line2)) and line2[12:16] == " C3'":
                            if (line1[17:20] and line1[17:20] in ['  A', '  U', '  G', '  C']) \
                                    and float(line2[22:26]) > float(line1[22:26]) + 3 \
                                    and line1[21:22] == line2[21:22]:
                                # Allows to select only RNA residues.
                                # and allows only residues separated by a distance of 3 or more to be taken and avoids testing the same pair of residues twice.
                                # and allows to consider only the intra-chain distance.
                                # --------------------------------------------------------------------------- #

                                pair_res = pair_res_format(line1, line2)  # Formalise the studied pair
                                # print(pair_res)

                                dist = calcul_dist(line1, line2)  # Calculates the distance between the two residues studied
                                # print(dist)

                                Gibbs_free_energy = linear_interpolation(train_dir, pair_res, dist, Gibbs_free_energy)

    # print(cpt)  # 12258
    end = time.time()
    elapsed = end - start

    print(f'Execution time script 3 : {elapsed:.1f}s')

    print(f"\t\nPredicted gibbs free energy of {test_file} :", Gibbs_free_energy)


def argparse_script3():
    parser = argparse.ArgumentParser(description='RNA folding problem; Output log ratio tabulate')
    parser.add_argument('-tr', '-train', type=str, required=True,
                        help='Input directory of RNA train file(s) (output directory script 1) [str]')
    parser.add_argument('-te', '-test', type=str, required=True,
                        help='Input RNA test file [pdb]')

    args = parser.parse_args()
    print('Input train directory:', args.tr)
    print('Input test file:', args.te)

    directory_train = args.tr
    file_test = args.te
    main_script3(directory_train, file_test)


if __name__ == '__main__':
    argparse_script3()
