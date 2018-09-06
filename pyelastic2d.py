#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
    Pyezo2D - python program to calculate elastic coefficients and piezoelectric constant of two-dimensional materials
        File:              pyezo2d.py
        Author:            Leandro Seixas Rocha
        Email:             leandro.seixas@mackenzie.br
        Institution:       MackGraphe - Graphene and Nanomaterials Research Center,
                           Mackenzie Presbyterian University
        Created:           29 August 2018
        Last modification: 06 September 2018

        To do list:        1. Calculate R^2
                           2. Calculate piezoelectric constant
'''

import sys, getopt
import numpy as np

def polynom(coeff,x,y):
    return 0.5*coeff[0]*x**2+0.5*coeff[1]*y**2+coeff[2]*x*y


# function fit_elastic_energy
def fit_elastic_energy(dft_data):
    '''
      Fitting the function U = 0.5*C_11*eps_x^2 + 0.5*C_22*eps_y^2 + C_12*eps_x*eps_y from dft data.

        dft_data should be in csv format with 3 columns: eps_x, eps_y, u:

          u = (E_tot(eps_x, eps_y) - E_tot(eps_x=0,eps_y=0))/area0,

        where E_tot are the total energy, and area0 is the unit cell area with no strains.
    ''' 
    X = dft_data[:,0]
    Y = dft_data[:,1]
    elastic_energy = dft_data[:,2]
    A = np.array([0.5*X**2, 0.5*Y**2, X*Y]).T
    elastic_coeff, r, rank, s = np.linalg.lstsq(A, elastic_energy, rcond=None)

    #rsqr = rsquare(elastic_energy,coeff)
    #result = [elastic_coeff, rsqr]
    #return result
    return elastic_coeff

# function report_result
def report_result(coeff):
    conversion_factor = 0.0160217662   # convert meV/angstrom^2 to N/m
    print("Elastic coefficients (N/m):")
    print("    C_11: {:.4f}".format(coeff[0]*conversion_factor))
    print("    C_22: {:.4f}".format(coeff[1]*conversion_factor))
    print("    C_12: {:.4f}".format(coeff[2]*conversion_factor))


# Reading file from argument

def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"he:p:",["elastic="])
   except getopt.GetoptError:
      print('usage: pyezo2d.py -[e,p] <inputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('usage: pyezo2d.py -[e,p] <inputfile>')
         sys.exit()
      elif opt in ("-e", "--elastic"):
         inputfile = arg
         # Reading data from file
         dft_data = np.genfromtxt(inputfile, delimiter=',')

         # Fitting quadratic polynomial from elastic energy data
         elastic_coeff=fit_elastic_energy(dft_data)

         # Print elastic coefficients
         report_result(elastic_coeff)
      elif opt in ("-p", "--polarization"):
         inputfile = arg
         # Reading data from file
         dft_data = np.genfromtxt(inputfile, delimiter=',')
  

if __name__ == "__main__":
   main(sys.argv[1:])

