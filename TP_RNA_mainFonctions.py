#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math


def calcul_dist(line1, line2):
    # Calculates the distance in Ångström (Å) between two residues.
    dist = math.sqrt((float(line1[30:38]) - float(line2[30:38])) ** 2 +
                     (float(line1[38:46]) - float(line2[38:46])) ** 2 +
                     (float(line1[46:54]) - float(line2[46:54])) ** 2)
    # print(dist)
    return dist


def pair_res_format(line1, line2):
    # Allows the formalisation of residue pairs.
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

    res1 = line1[17:20].replace(" ", "")
    res2 = line2[17:20].replace(" ", "")

    pair_res = res1 + res2

    if not any(res in pair_res for res in base_pairs):
        # If the pair is not formalised as desired (according to the "Base_pair" dictionary), the order of the residues is reversed.
        pair_res = res2 + res1

    return pair_res
