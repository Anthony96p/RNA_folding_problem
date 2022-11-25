#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import math
import time
import shutil
import argparse
from os import listdir

from TP_RNA_mainFonctions import calcul_dist
from TP_RNA_mainFonctions import pair_res_format


def dict_obs_freq():
    # Creation of the dictionary that will group the counts of the number of residues per pair for observed probability.
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

    dict_obs = {}

    for i in base_pairs:
        dict_obs[i] = {}
        for j in range(21):
            dict_obs[i][j] = int()
        dict_obs[i]['Nij'] = int()
    # print(dict_obs)

    return dict_obs


def compute_obs_freq(pair_res, dist, obs_freq):
    # Count the number of residues per pair for each distance.
    if dist <= 21:  # Restriction of distances counted to 21 Å.
        obs_freq[pair_res][int(math.floor(dist))] = obs_freq[pair_res][int(math.floor(dist))] + 1

    return obs_freq


def end_obs_freq(fin_obs_freq):
    # Will add the Nij values for each peer and calculate the observed probability.
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

    for i in base_pairs:
        for j in range(21):
            fin_obs_freq[i]['Nij'] = fin_obs_freq[i]['Nij'] + fin_obs_freq[i][j]

    # Displays the final observed count table :
    # for i in base_pairs:
    #     print(i, '=>', fin_obs_freq[i])

    for i in base_pairs:
        for j in range(21):
            # print ('NijR/Nij',float(obs_freq[i][j]) / float(obs_freq[i]['Nij']))
            fin_obs_freq[i][j] = float(fin_obs_freq[i][j]) / float(fin_obs_freq[i]['Nij'])

    # Displays the final observer probability table :
    # for i in base_pairs:
    #     print(i, '=>', fin_obs_freq[i])

    return fin_obs_freq


def dict_ref_freq():
    # Creation of the dictionary that will group the counts of the number of residues per pair for reference probability.
    dict_ref = {}

    for j in range(21):
        dict_ref[j] = int()
    dict_ref['Nxx'] = int()
    # print(dict_ref)

    return dict_ref


def compute_ref_freq(dist, ref_freq):
    # Count the number of residues for each distance.
    if dist <= 21:  # Restriction of distances counted to 21 Å.
        ref_freq[int(math.floor(dist))] = ref_freq[int(math.floor(dist))] + 1

    return ref_freq


def end_ref_freq(fin_ref_freq):
    # Will add the Nxx value for each distance referenced probability.
    c = 0
    for j in range(21):
        fin_ref_freq['Nxx'] = fin_ref_freq['Nxx'] + fin_ref_freq[j]
        c = c + fin_ref_freq[j]

    # Displays the final referenced count table :
    # print("Ref=>",fin_ref_freq)

    for j in range(21):
        fin_ref_freq[j] = float(fin_ref_freq[j]) / float(fin_ref_freq['Nxx'])

    # Displays the final reference probability table :
    # print("Ref=>", fin_ref_freq)

    return fin_ref_freq


def compute_log_ratio(obs_freq, ref_freq, dir):
    # Calculates all log ratios
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

    for i in base_pairs:
        for j in range(21):
            try:
                # print('fijOBS(r)', obs_freq[i][j])
                # print('fijREF(r)', ref_freq[j])
                log_ratio = -math.log(obs_freq[i][j] / ref_freq[j])
                if log_ratio > 10:
                    log_ratio = 10
                create_file(log_ratio, i, j, dir)
            except ZeroDivisionError:
                # print('fijREF(r)', ref_freq[j])
                # print('div0')
                ratio = 10
                create_file(ratio, i, j, dir)
            except ValueError:
                # print('fijOBS(r)', obs_freq[i][j])
                # print('log(0)')
                ratio = 10
                create_file(ratio, i, j, dir)


def create_dir(dir):
    # Delete directory if exist
    try:
        shutil.rmtree(dir)
    except:
        pass

    # Create new directory
    os.makedirs(dir)


def create_file(log_ratio, pair_res, dist_intervals, path_dir):
    # Created a scoring values file for each pair of residues
    # print (log_ratio)
    with open(path_dir + '/' + pair_res + '.txt', 'a+') as file_scoring:
        file_scoring.seek(0)
        data = file_scoring.read(100)
        if len(data) > 0:
            file_scoring.write("\n")
        file_scoring.write(str(dist_intervals) + '\t' + str(log_ratio))


def main_script1(pdb_file, Dir, timepdb):
    start = time.time()

    obs_freq = dict_obs_freq()  # Creation of the Observer Frequency Dictionary
    ref_freq = dict_ref_freq()  # Creation of the Referenced Frequency Dictionary

    # The script opens and browses the file 2 times in order to browse all possible peers
    for pdb in listdir(pdb_file + "/"):
        print(f'Current train RNA : {pdb}')
        if timepdb is True:
            startpdb = time.time()
        with open(pdb_file + "/" + pdb, 'r') as file1:
            for line1 in file1:
                if (re.match('^ATOM', line1)) and line1[12:16] == " C3'":
                    # --------------------------------------------------------------------------- #
                    with open(pdb_file + "/" + pdb, 'r') as file2:
                        for line2 in file2:
                            if (re.match('^ATOM', line2)) and line2[12:16] == " C3'":
                                if (line1[17:20] and line1[17:20] in ['  A', '  U', '  G', '  C']) \
                                        and float(line2[22:26]) > float(line1[22:26]) + 3 \
                                        and line1[21:22] == line2[21:22]:
                                    # Allows to select only RNA residues.
                                    # and allows only residues separated by a distance of 3 or more to be taken and avoids testing the same pair of residues twice.
                                    # and allows to consider only the intra-chain distance.
                                    # --------------------------------------------------------------------------- #

                                    pair_res = pair_res_format(line1, line2)  # Formalise the studied pair.
                                    # print(pair_res)

                                    dist = calcul_dist(line1, line2)  # Calculates the distance between the two residues studied.
                                    # print(dist)

                                    # Compute the observed frequencies: 10 × 20 distances intervals (0 to 20 Å) :
                                    obs_freq = compute_obs_freq(pair_res, dist, obs_freq)

                                    # Compute the reference frequency (= the “XX” pair) :
                                    ref_freq = compute_ref_freq(dist, ref_freq)
        if timepdb is True:
            endpdb = time.time()
            elapsedpdb = endpdb - startpdb

            print(f'File execution time  : {elapsedpdb:.1f}s\n')

    fin_obs_freq = end_obs_freq(obs_freq)  # Calculates the frequencies observe total
    fin_ref_freq = end_ref_freq(ref_freq)  # Calculates the frequencies reference total

    # Create output directory
    create_dir(Dir)

    # Compute the log-ratio of the two frequencies :
    compute_log_ratio(fin_obs_freq, fin_ref_freq, Dir)

    end = time.time()
    elapsed = end - start

    print(f'\nExecution time script 1 : {elapsed:.1f}s')


def argparse_script1():
    parser = argparse.ArgumentParser(description='RNA folding problem; Output log ratio tabulate')
    parser.add_argument('-i', '-input', type=str, required=True,
                        help='Input directory with RNA file (pdb) [str]')
    parser.add_argument('-o', '-output', type=str, required=True,
                        help='Output directory of log ratio file (tabular) [str]')
    parser.add_argument('-d', '-details', type=lambda x: (str(x).lower() == 'true'), default=False,
                        required=False,
                        help='An optional boolean switch to display the duration of each RNA files training [boolean]')

    args = parser.parse_args()

    print('Input directory:', args.i)
    print('Output directory:', args.o)
    print('Print time details:', args.d, "\n")

    file_dir = args.i
    directory = args.o
    timepdb = args.d
    main_script1(file_dir, directory, timepdb)


if __name__ == '__main__':
    argparse_script1()
