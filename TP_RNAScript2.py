#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from colorhash import ColorHash


def main_script2(directory, save_png, print_png, results_file=False):
    start = time.time()

    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

    for i in base_pairs:
        # Affiche chaque le profile interaction pour chaque pair
        c = ColorHash(i)  # Génère une couleur aléatoire pour chaque graph.

        # Va lire le fichier de scoring value et afficher ou sauvegarder sous forme de graph
        Res = np.genfromtxt(str(directory) + '/' + i + '.txt')
        plt.plot(Res[:, 0], Res[:, 1], 'o-', color=c.hex, label=i)

        plt.legend(loc='upper right')
        plt.title(i + ' interaction profile of  RNA')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Energy (Pairwise log-ratio)')
        plt.xlim([-1, 20])

        if save_png is True:
            plt.savefig(directory + '/' + i + '.png')
        else:
            try:
                os.remove(directory + '/' + i + '.png')
            except:
                pass
        if print_png is True:
            plt.show()

    end = time.time()
    elapsed = end - start

    print(f'Execution time script 2 : {elapsed:.1f}s')

    if results_file is True:
        with open('Results_RNA_folding_problem.txt', 'a+') as file_scoring:
            file_scoring.write(f'\nExecution time script 2 : {elapsed:.1f}s')



def argparse_script2():
    parser = argparse.ArgumentParser(description='RNA folding problem')
    parser.add_argument('-i', '-input', type=str, required=True,
                        help='Input directory (tabular) [str]')
    parser.add_argument('-s', '-savepng' , type=lambda x: (str(x).lower() == 'true'), default=False,
                        required=False, help='An optional boolean switch to output matplolib plot (png) [boolean]')
    parser.add_argument('-p', '-printpng', type=lambda x: (str(x).lower() == 'true'), default=True,
                        required=False, help='An optional boolean switch to print matplolib plot [boolean]')

    args = parser.parse_args()

    print('Output directory:', args.i)
    print('Print plot png:', args.p)
    print('Save plot png:', args.s)

    directory = args.i
    save_png = args.s
    print_png = args.p

    main_script2(directory, save_png, print_png)


if __name__ == '__main__':
    argparse_script2()
