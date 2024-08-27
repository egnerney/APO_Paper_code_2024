#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 17:15:37 2024

@author: edne8319
"""

import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.optimize import differential_evolution, least_squares
#from scipy.integrate import newton_cotes
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import time

#import os


def int_tabulated(x, f):
    """
    Integrates a tabulated set of data (x, f) using a 5-point Newton-Cotes formula
    on the closed interval [min(x), max(x)].

    Parameters:
    x (array_like): The tabulated x-value data. This data may be irregularly gridded.
    f (array_like): The tabulated f-value data corresponding to x.
    sort (bool): If True, the x and f data will be sorted into ascending order. Default is False.
    double (bool): If True, computations are done in double precision. Default is False.

    Returns:
    float: The integral of f over x.
    """

    # Determine the number of segments as an integer multiple of 4
    num_segments = len(x) - 1
    if num_segments % 4 != 0:
        num_segments += 4 - num_segments % 4  # Adjust to the nearest multiple of 4

    # Define the uniform step size based on extended segments
    x_min, x_max = np.min(x), np.max(x)
    h = (x_max - x_min) / num_segments

    # Interpolate using a cubic spline, default to linear if not enough points
    #from scipy.interpolate import interp1d
    interp_kind = 'cubic' if len(x) > 3 else 'linear'
    interp_func = interp1d(x, f, kind=interp_kind, fill_value="extrapolate", bounds_error=False)
    
    # Evaluate the interpolates at regular intervals
    x_gridz = np.linspace(x_min, x_max, num_segments + 1)
    f_interp = interp_func(x_gridz)

    # Compute the integral using the 5-point Newton-Cotes formula
    integral = 0.0
    for i in range(0, num_segments, 4):
        integral += 2 * h * (7 * (f_interp[i] + f_interp[i + 4]) +
                             32 * (f_interp[i + 1] + f_interp[i + 3]) +
                             12 * f_interp[i + 2]) / 45

    return integral


def given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(Tec_a, nec_a, nsp_a, ns2p_a, nop_a, rho_max, rho_min, drho, yptsi_, nel_, tec_):
    num_emiss = 8
    num_elements = len(nec_a)
    epsilon_a = np.zeros((num_elements, num_emiss))

    for j in range(num_emiss):
        tr = np.array([nel_, tec_])
        epsilon_a[:, j] = griddata(tr.T, yptsi_[:, j].flatten(), (nec_a, Tec_a), method='linear')

    epsilon_a[:, 0] *= 1e-6 * ns2p_a
    epsilon_a[:, 1] *= 1e-6 * nop_a
    epsilon_a[:, 2] *= 1e-6 * nop_a
    epsilon_a[:, 3] *= 1e-6 * nsp_a
    epsilon_a[:, 4] *= 1e-6 * nsp_a
    epsilon_a[:, 5] *= 1e-6 * ns2p_a
    epsilon_a[:, 6] *= 1e-6 * nsp_a
    epsilon_a[:, 7] *= 1e-6 * nsp_a

    rayleighs = np.zeros((num_elements, 8))
    conversion = 7.1492e9
    n_elem_s = int(2 * rho_max / drho) + 1
    s_LOS = drho * np.arange(n_elem_s)
    n_elem_x0_LOS = num_elements#int((rho_max - rho_min) / drho) + 1
    y0 = -rho_max
    x_LOS = drho * np.arange(n_elem_x0_LOS) + rho_min
    y_LOS = y0 + s_LOS
    
    x_LOS =    x_LOS[::-1]
    rho_LOS = np.sqrt(x_LOS[:, np.newaxis] ** 2 + y_LOS[np.newaxis, :] ** 2)

 

    for idx in range(num_elements):
        idx_want = np.where((rho_LOS[idx, :] <= rho_max) & (rho_LOS[idx, :] >= rho_min))[0]
        if len(idx_want) >= 2:
            for linenum in range(8):
                emiss_s = interp1d(x_LOS, epsilon_a[:, linenum], kind='linear')(rho_LOS[idx, idx_want])
                rayleighs[idx, linenum] = conversion * int_tabulated(s_LOS[idx_want], emiss_s) #simps(emiss_s, s_LOS[idx_want])
        elif len(idx_want) == 1:
            path_length = s_LOS[idx_want[0] + 1] - s_LOS[idx_want[0]] if idx_want[0] + 1 < len(s_LOS) else drho
            for linenum in range(8):
                emiss_value = interp1d(x_LOS, epsilon_a[:, linenum], kind='linear')(rho_LOS[idx, idx_want])[0]
                rayleighs[idx, linenum] = conversion * emiss_value * path_length
        else:
            rayleighs[idx, :] = 0

    return rayleighs



def residuals_function(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, totchmix_for_apo):
    """
    Calculate residuals for the least_squares optimization.

    Parameters:
    p : array-like
        Parameters for the functional forms.
    x, y : array-like
        Independent and dependent variables for fitting.
    err : array-like
        Errors associated with y-values.
    x_grid : array-like
        The x-values of the data points.
    x_step_size : float
        The step size of x_grid.
    yptsi_in_2dto1d, nel2dto1d, tec2dto1d : array-like
        Additional data arrays for interpolation and fitting.
    nspmix, ns2pmix, nopmix : array-like
        Mixing ratios for S+, S++, and O+.

    Returns:
    residuals : array-like
        The residuals for the optimization.
    """
    nx_grid = len(x_grid)
    nec_funcform_test = np.zeros(nx_grid)
    nsp_funcform_test = np.zeros(nx_grid)
    ns2p_funcform_test = np.zeros(nx_grid)
    nop_funcform_test = np.zeros(nx_grid)
    flags2p = np.zeros(nx_grid)
    flagsp = np.zeros(nx_grid)

    rhocnsp = p[0]
    rhoRnsp = p[18]#rhoRnsp_
    rhoORnsp = p[23]#rhoORne_
    G = p[1]
    Hin = p[2]
    Hout = p[3]
    I = p[4]
    Jin = p[5]
    Jout = p[6]
    L = p[24]#3.37 + nec_pl  # assuming nec_pl is p[14]

    K = G * np.exp(-((rhoORnsp - rhocnsp)/Hout)**2.) + I * np.exp(-((rhoORnsp - rhoRnsp)/Jout)**2.)

    rhoRns2p = p[19]#rhoRne_
    rhoORns2p = p[23]#rhoORne_
    C2 = p[7]
    Din2 = p[8]
    Dout2 = p[9]
    F2 = p[25]# nec_pl + 1.13  # assuming nec_pl is p[14]

    E2 = C2 * np.exp(-((rhoORns2p - rhoRns2p)/Dout2)**2.)

    rhocnop = rhocnsp
    rhoRnop = p[20]#rhoRne_
    rhoORnop = p[23]#rhoORne_
    G2 = p[10]
    Hin2 = Hin
    Hout2 = Hout
    I2 = p[11]
    Jin2 = p[12]
    Jout2 = p[13]
    L2 = p[26]# nec_pl + 0.96  # assuming nec_pl is p[14]

    K2 = G2 * np.exp(-((rhoORnop - rhocnop)/Hout2)**2.) + I2 * np.exp(-((rhoORnop - rhoRnop)/Jout2)**2.)

    tec = np.zeros(nx_grid)
    A_TE1 = p[14]
    B_te1 = p[15]
    B_te2 = p[16]
    A_TE2 = p[17]

    beta_ = np.log(A_TE2 / A_TE1) / np.log(p[22] /p[21])

    for idx in range(nx_grid):
        nsp_funcform_test[idx] = K * (x_grid[idx] / rhoORnsp)**(-L)
        ns2p_funcform_test[idx] = E2 * (x_grid[idx] / rhoORns2p)**(-F2)
        nop_funcform_test[idx] = K2 * (x_grid[idx] / rhoORnop)**(-L2)

        if x_grid[idx] < rhoORnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jout)**2.)
        if x_grid[idx] < rhoORns2p:
            ns2p_funcform_test[idx] = C2 * np.exp(-((x_grid[idx] - rhoRns2p)/Dout2)**2.)
        if x_grid[idx] < rhoORnop:
            nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hout2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jout2)**2.)
            
        if x_grid[idx] < rhoRnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)
        if x_grid[idx] < rhoRns2p:
            ns2p_funcform_test[idx] = C2 * np.exp(-((x_grid[idx] - rhoRns2p)/Din2)**2.)
        if x_grid[idx] < rhoRnop:
            nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hout2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jin2)**2.)

        if x_grid[idx] < rhocnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hin)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)

        if x_grid[idx] < rhocnop:
            nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hin2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jin2)**2.)

        nec_funcform_test[idx] = (nsp_funcform_test[idx] + 2. * ns2p_funcform_test[idx] + nop_funcform_test[idx]) / totchmix_for_apo[idx]

        if x_grid[idx] <= p[21]:
            tec[idx] = A_TE1 * (x_grid[idx] / p[21])**B_te1
        elif x_grid[idx] > p[22]:
            tec[idx] = A_TE2 * (x_grid[idx] / p[22])**B_te2
        else:
            tec[idx] = A_TE1 * (x_grid[idx] / p[21])**beta_

        if nec_funcform_test[idx] < 0.000001:
            nspmixtemp = nsp_funcform_test[idx] / nec_funcform_test[idx]
            ns2pmixtemp = ns2p_funcform_test[idx] / nec_funcform_test[idx]
            nopmixtemp = nop_funcform_test[idx] / nec_funcform_test[idx]
            nec_funcform_test[idx] = 0.0000011
            nsp_funcform_test[idx] =  nspmixtemp * nec_funcform_test[idx]
            ns2p_funcform_test[idx] = ns2pmixtemp *  nec_funcform_test[idx]
            nop_funcform_test[idx] = nopmixtemp *  nec_funcform_test[idx]
            
        if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.05) and (x_grid[idx] <= 5.3)):
            flags2p[idx] = 1.
        if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.20) and (x_grid[idx] <= 5.73)):
            flags2p[idx] = 1.
        if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.19) and (x_grid[idx] <= 5.68)):
            flags2p[idx] = 1.
        if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.18) and (x_grid[idx] <= 5.65)):
            flags2p[idx] = 1.
        if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.16) and (x_grid[idx] <= 5.59)):
            flags2p[idx] = 1.
        if (((nsp_funcform_test[idx] / nec_funcform_test[idx]) < 0.55) and (x_grid[idx] <= 5.0)):
            flagsp[idx] = 1.

    model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
        tec, nec_funcform_test, nsp_funcform_test, ns2p_funcform_test, nop_funcform_test, max(x_grid), min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
    )

    idx_want = [0, 1, 2, 3, 4,  5, 6, 7]
    model = model[:, idx_want].flatten()

    penalt = 1e6
    flagtots2p = np.sum(flags2p)
    flagtotsp = np.sum(flagsp)


    if flagtots2p > 0:
        model = model +  penalt * flagtots2p / len(model)
    if flagtotsp > 0:
        model = model +  penalt * flagtotsp / len(model)

    result = np.sum(((y - model) / err) ** 2.)# (y - model) / err
    return result


def residuals_function_simult(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, nspmix, ns2pmix, nopmix, totchmix_for_apo):
    
    nx_grid = len(x_grid)
   # nec_tempi = p[0:nx_grid]*nspmix
    nsp_tempi = p[0*nx_grid:1 * nx_grid]
    ns2p_tempi = p[1 * nx_grid:2 * nx_grid]
    nop_tempi = p[2 * nx_grid:3 * nx_grid]
    tec_tempi = p[3 * nx_grid:4 * nx_grid]
    nec_tempi = (nsp_tempi + 2. * ns2p_tempi + nop_tempi) / totchmix_for_apo


    model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
        tec_tempi, nec_tempi, nsp_tempi, ns2p_tempi, nop_tempi, max(x_grid), min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
    )

    idx_want = [0, 1, 2, 3, 4, 5, 6, 7]
    model = model[:, idx_want].flatten()


    result = (y - model) / err
    return result


def parallel_jacobian(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, totchmix_for_apo, n_jobs=-1):
    """
    Calculate the Jacobian matrix using parallel processing.

    Parameters:
    p : array-like
        Parameters for the functional forms.
    x, y : array-like
        Independent and dependent variables for fitting.
    err : array-like
        Errors associated with y-values.
    x_grid : array-like
        The x-values of the data points.
    x_step_size : float
        The step size of x_grid.
    yptsi_in_2dto1d, nel2dto1d, tec2dto1d : array-like
        Additional data arrays for interpolation and fitting.
    nspmix, ns2pmix, nopmix : array-like
        Mixing ratios for S+, S++, and O+.
    n_jobs : int
        Number of parallel jobs to run.

    Returns:
    jacobian : array-like
        The Jacobian matrix.
    """
    eps = np.sqrt(np.finfo(float).eps)

    def compute_jacobian_column(i):
        p1 = p.copy()
        p1[i] += eps
        p2 = p.copy()
        p2[i] -= eps
        f1 = residuals_function(p1, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, totchmix_for_apo)
        f2 = residuals_function(p2, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, totchmix_for_apo)
        return (f1 - f2) / (2. * eps)

    jacobian_columns = Parallel(n_jobs=n_jobs)(delayed(compute_jacobian_column)(i) for i in range(len(p)))
    jacobian = np.column_stack(jacobian_columns)
    return jacobian



def parallel_jacobian_simult(p,  x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, nspmix, ns2pmix, nopmix, totchmix_for_apo, n_jobs=-1):
    """
    Calculate the Jacobian matrix using parallel processing.

    Parameters:
    p : array-like
        Parameters for the functional forms.
    x, y : array-like
        Independent and dependent variables for fitting.
    err : array-like
        Errors associated with y-values.
    x_grid : array-like
        The x-values of the data points.
    x_step_size : float
        The step size of x_grid.
    yptsi_in_2dto1d, nel2dto1d, tec2dto1d : array-like
        Additional data arrays for interpolation and fitting.
    nspmix, ns2pmix, nopmix : array-like
        Mixing ratios for S+, S++, and O+.
    n_jobs : int
        Number of parallel jobs to run.

    Returns:
    jacobian : array-like
        The Jacobian matrix.
    """
    eps = np.sqrt(np.finfo(float).eps)

    def compute_jacobian_column_simult(i):
        p1 = p.copy()
        p1[i] += eps
        p2 = p.copy()
        p2[i] -= eps
        f1 = residuals_function_simult(p1, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, nspmix, ns2pmix, nopmix, totchmix_for_apo)
        f2 = residuals_function_simult(p2, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, nspmix, ns2pmix, nopmix, totchmix_for_apo)
        return (f1 - f2) / (2. * eps)

    jacobian_columns = Parallel(n_jobs=n_jobs)(delayed(compute_jacobian_column_simult)(i) for i in range(len(p)))
    jacobian = np.column_stack(jacobian_columns)
    return jacobian




if __name__ == '__main__':
    yptsi_in = np.loadtxt('yptsi_8_APO_LINES_optical_vary_nec_tec_lookuptable_66x58x8_CHIANTI_10.1.csv', delimiter=',')
    yptsi_in = yptsi_in.reshape((8, 58, 66)).transpose(2, 1, 0)

    x_grid_coadded = np.loadtxt('x_grid_to_plot_60_points.csv', delimiter=',')
    xgrid_strings = ['x_grid_night1.csv','x_grid_night2.csv','x_grid_night3.csv','x_grid_night4.csv','x_grid_dusk_night5.csv','x_grid_night6.csv']
    nightname_strings = ['night1','night2','night3','night4','night5','night6']
    spec_string = ['Dusk_night1_grid_to_plot_71x8.csv','Dusk_night2_grid_to_plot_71x8.csv','Dusk_night3_grid_to_plot_68x8.csv','Dusk_night4_grid_to_plot_67x8.csv','Dusk_night5_grid_to_plot_66x8.csv','Dusk_night6_grid_to_plot_64x8.csv']
    
    err_spec_string =  ['err_Dusk_night1_grid_to_plot_71x8.csv','err_Dusk_night2_grid_to_plot_71x8.csv','err_Dusk_night3_grid_to_plot_68x8.csv','err_Dusk_night4_grid_to_plot_67x8.csv','err_Dusk_night5_grid_to_plot_66x8.csv','err_Dusk_night6_grid_to_plot_64x8.csv']
    
    out_string = ['outs_and_errors_dusk_night1_funcform_fit_71x16.csv','outs_and_errors_dusk_night2_funcform_fit_71x16.csv','outs_and_errors_dusk_night3_funcform_fit_68x16.csv','outs_and_errors_dusk_night4_funcform_fit_67x16.csv','outs_and_errors_dusk_night5_funcform_fit_66x16.csv','outs_and_errors_dusk_night6_funcform_fit_64x16.csv']
    
    idxss = [1]
    for w in idxss:
        # Load and preprocess data
        
        #dawn_grid1 = np.loadtxt('dusk_grid_to_plot_60x8.csv', delimiter=',')
        #dawn_grid = np.transpose(dawn_grid1)
        

        #err_dawn_grid1 = np.loadtxt('err_dusk_grid_to_plot_60x8.csv', delimiter=',')
       # err_dawn_grid = np.transpose(err_dawn_grid1)
       
        x_grid = np.loadtxt(xgrid_strings[w], delimiter=',')
        dawn_grid1 = np.loadtxt(spec_string[w], delimiter=',')

        dawn_grid = np.transpose(dawn_grid1)
        
        err_dawn_grid1 = np.loadtxt(err_spec_string[w], delimiter=',')

        err_dawn_grid = np.transpose(err_dawn_grid1)

        nel = np.concatenate(([0.000001, 0.00001, 0.0001, 0.001, 0.1], np.arange(9.) + 1., 10. * np.arange(9.) + 10., 100. * np.arange(4.) + 100., 250. * np.arange(39.) + 500.))
        Tec = np.concatenate(([0.01], 0.1 * np.arange(9.) + 0.1, 0.5 * np.arange(17.) + 1., np.arange(10.) + 10., 5. * np.arange(8.) + 20., 20. * np.arange(7.) + 60., 80. * np.arange(6.) + 200.))

        n_ne = len(nel)
        n_te = len(Tec)
        nel2dto1d = np.repeat(nel, n_te)
        tec2dto1d = np.tile(Tec, n_ne)
        yptsi_in_2dto1d = yptsi_in.reshape((n_ne * n_te, 8))

        XTITTLEE = 'Dusk $\\rho_c$ ($R_J$)'
        idx_want = [0, 1, 2, 3, 4, 5, 6, 7]
        dawn_grid_temp = dawn_grid[:, idx_want]
        err_dawn_grid_temp = err_dawn_grid[:, idx_want]

        ndawn_grid = dawn_grid_temp.size
        y = dawn_grid_temp.flatten()
        err = err_dawn_grid_temp.flatten()
        x = np.ones(ndawn_grid)

        x_step_size = round(abs(x_grid[1] - x_grid[0]), 5)
        x_grid_new = round(abs(x_grid[1] - x_grid[0]), 5) * np.arange(len(x_grid)) + min(x_grid)
        x_grid = x_grid_new[::-1]

        out_fits1 = np.loadtxt(out_string[w], delimiter=',')
        out_fits = np.transpose(out_fits1)

        totchmix_for_apo = np.loadtxt('totchmix_for_apo_7pointavg_dusk_unshifted.csv', delimiter=',')
        
        xtoplot = x_grid_coadded[20:46]
        aptochmix_to_plot =  interp1d(np.array([x_grid_coadded[20],x_grid_coadded[45]]), np.array([totchmix_for_apo[20],totchmix_for_apo[45]]), kind='linear')(xtoplot)
        totchmix_for_apo = np.concatenate([totchmix_for_apo[0:20],aptochmix_to_plot,totchmix_for_apo[46:]])
        
        interp_function = interp1d(x_grid_coadded, totchmix_for_apo, kind='linear', bounds_error=False, fill_value="extrapolate")

        # Interpolate to all bins except the last one
        totchmix_for_apo = interp_function(x_grid[:-1])

        # Append 1.0 for the last bin
        totchmix_for_apo = np.append(totchmix_for_apo, 1.0)
        Tec_fixed_found = out_fits[:, 0]
        nec_fixed_found = out_fits[:, 1]
        nsp_fixed_found = out_fits[:, 2]
        ns2p_fixed_found = out_fits[:, 3]
        nop_fixed_found = out_fits[:, 4]
        nspmix = out_fits[:, 5]
        ns2pmix = out_fits[:, 6]
        nopmix = out_fits[:, 7]
        Tec_error = out_fits[:, 8]
        nec_error = out_fits[:, 9]
        nsp_error = out_fits[:, 10]
        ns2p_error = out_fits[:, 11]
        nop_error = out_fits[:, 12]
        nspmix_error = out_fits[:, 13]
        ns2pmix_error = out_fits[:, 14]
        nopmix_error = out_fits[:, 15]

        model_spectra = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
           Tec_fixed_found, nec_fixed_found, nec_fixed_found * nspmix, nec_fixed_found * ns2pmix, nec_fixed_found * nopmix,
            np.max(x_grid), np.min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
        )

        fig, axs = plt.subplots(8, 1, figsize=(10, 20), sharex=True)

        for i, ax in enumerate(axs):
            ax.errorbar(x_grid, dawn_grid[:, i], yerr=err_dawn_grid[:, i], fmt='o', label=f'Line {i} Data')
            ax.plot(x_grid, model_spectra[:, i], label=f'Line {i} Fit')
            ax.set_ylabel(f'Rayleighs (Line {i})')
            ax.legend()
            ax.grid(True)

        axs[-1].set_xlabel(XTITTLEE)
        plt.tight_layout()
        #plt.savefig('initial_emiss_funcform_for_python_dusk.png')
        plt.show()


        nec_pl = 3.
        rhocne = 5.06
        rhoRne = 5.7
            
        rhoRns2p_ =  5.7
        Din2max = 0.2
            
        rhoORne = 5.82
        A = 1960.
        Bin = 0.3
        Bout = 0.3
        C = 3430.
        Din = 0.2
        Dout = 0.2
        F = nec_pl
        E = A*np.exp(-((rhoORne - rhocne)/Bout) **2.) + C*np.exp(-((rhoORne - rhoRne)/Dout) **2.)


        rhocnsp =    rhocne
        rhoRnsp = 5.66
        rhoORnsp = rhoORne
        G = 600.
        Hin= 0.11
        Hout= 0.33
        I = 400.
        Jin = 0.142
        Jout = 0.192
        L = 3.37 + nec_pl

        K = G*np.exp(-((rhoORnsp - rhocnsp)/Hout) **2.) + I*np.exp(-((rhoORnsp - rhoRnsp)/Jout) **2.)
            
            
        rhoRne_ = rhoRne
        rhoRnsp_  = rhoRnsp
        rhoORne_ = rhoORne
        rhocne_ = rhocne

        Bin_ = Bin
        Bout_ = Bout


        rhoRns2p = rhoRns2p_
        rhoORns2p = rhoORne_
        A = 1960.
        Bin = 0.3
        Bout = 0.3
        C2 = 520.
        Din2 = 0.09
        Dout2 = 0.91
        F2 = nec_pl  + 1.13
        E2 =  C2*np.exp(-((rhoORns2p - rhoRns2p)/Dout2) **2.)

        rhocnop = rhocnsp
        rhoRnop = rhoRne_ 
        rhoORnop = rhoORne_ 
        G2 = 400.
        Hin2 = Hin 
        Hout2 = Hout 
        I2 = 600.
        Jin2 = 0.119
        Jout2 = 0.27
        L2 = nec_pl + 0.96

        K2 = G2*np.exp(-((rhoORnop - rhocnop)/Hout2) **2.) + I*np.exp(-((rhoORnop - rhoRnop)/Jout2) **2.)
            
            
        A_TE1 = 2.1
           
        B_te1 = 6.76
        B_te2 = 1.823
            
        A_TE2 = 5.
            
        rho4069peak = 5.58
            
        beta_ = np.log(A_TE2 / A_TE1) / np.log(  rhoORne_ /   rho4069peak)

        nec_lower_bounds = np.full(len(x_grid), 0.000001)
        nec_upper_bounds = np.full(len(x_grid), 10000.0)
        tec_lower_bounds = np.full(len(x_grid), 0.1)
        tec_upper_bounds = np.full(len(x_grid), 8.0)

        #lower_bounds = np.concatenate([nec_lower_bounds, tec_lower_bounds])
        #upper_bounds = np.concatenate([nec_upper_bounds, tec_upper_bounds])


        #lower_bounds = [4.9,300., 0.1, 0.1, 100., 0.05, 0.05, 100.,0.05,0.05,100.,100.,0.05,0.05, 1., 0.01, 0.01, 3. ]
        #upper_bounds = [5.1,1400., 0.5,0.6, 1000.,0.5, 1.5, 1000.,0.8,2.5,1000.,1000.,0.5,2.5,   4.,10., 10.,8.   ]
        # Initialize parameters (these should be aligned with the p[] indices used in the residuals function)
        #p0 = [5.27, 1298, 0.2, 0.3, 743, 0.08, 0.1, 520, 0.09, 0.91, 470, 1106, 0.1, 0.2, 3.37, 1.0, 2.0, 3.0]

        #p0 = [rhocnsp, G,Hin,Hout,I,Jin,Jout,A_TE1,B_te1,B_te2, A_TE2,A,Bout,C,Din,Dout]
        p0 = [5.04, 800.,0.198,Hout,600.,Jin,Jout,650.,Din2,Dout2, 350., I2,Jin2,Jout2, A_TE1,B_te1 ,B_te2, A_TE2,5.66,5.7,5.7,5.58,5.82,5.82,6.37,4.13,3.96]
        start_time = time.time()

        bounds = [(5.03, 5.06), (600., 900.), (0.15, 0.5), (0.1, 0.6), (300.,600.), (0.05, 0.5), (0.05, 1.5), (300., 700.), (0.05, 0.8), (0.05, 2.5), (300., 700.), (300., 700.), (0.05, 0.5), (0.05, 2.5), (1., 3.), (0.01, 10.), (0.01, 10.), (3., 5.),(5.55, 5.63), (5.73, 5.78), (5.68, 5.75),(5.25, 5.65), (5.66, 5.95),(5.75,6.05),(3.,10.),(2.,8.),(0.95,8.)]




        result = differential_evolution(
                residuals_function,
                bounds=bounds,updating='deferred',
                args=(x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, totchmix_for_apo),
                workers=-1
            )


        popt = result.x
        #perr = np.sqrt(np.diag(result.jac.T @ result.jac))

        end_time = time.time()

        print("--- %s seconds ---" % (end_time - start_time))

        nx_grid = len(x_grid)
        nec_funcform_test = np.zeros(nx_grid)
        nsp_funcform_test = np.zeros(nx_grid)
        ns2p_funcform_test = np.zeros(nx_grid)
        nop_funcform_test = np.zeros(nx_grid)
        flags2p = np.zeros(nx_grid)

        rhocnsp = popt[0]
        rhoRnsp = popt[18] #rhoRnsp_
        rhoORnsp = popt[23]#rhoORne_
        G = popt[1]
        Hin = popt[2]
        Hout = popt[3]
        I = popt[4]
        Jin = popt[5]
        Jout = popt[6]
        L =popt[24]# 3.37 + nec_pl  

        K = G * np.exp(-((rhoORnsp - rhocnsp)/Hout)**2.) + I * np.exp(-((rhoORnsp - rhoRnsp)/Jout)**2.)

        rhoRns2p = popt[19] #rhoRns2p_
        rhoORns2p = popt[23]#rhoORne_
        C2 = popt[7]
        Din2 = popt[8]
        Dout2 = popt[9]
        F2 = popt[25]# nec_pl + 1.13  

        E2 = C2 * np.exp(-((rhoORns2p - rhoRns2p)/Dout2)**2.)

        rhocnop = rhocnsp
        rhoRnop = popt[20]#rhoRne_
        rhoORnop =popt[23]# rhoORne_
        G2 = popt[10]
        Hin2 = Hin
        Hout2 = Hout
        I2 = popt[11]
        Jin2 = popt[12]
        Jout2 = popt[13]
        L2 = popt[26]# nec_pl + 0.96  

        K2 = G2 * np.exp(-((rhoORnop - rhocnop)/Hout2)**2.) + I2 * np.exp(-((rhoORnop - rhoRnop)/Jout2)**2.)

        tec = np.zeros(nx_grid)
        A_TE1 = popt[14]
        B_te1 = popt[15]
        B_te2 = popt[16]
        A_TE2 = popt[17]

        beta_ = np.log(A_TE2 / A_TE1) / np.log(popt[22] / popt[21])

        for idx in range(nx_grid):
            nsp_funcform_test[idx] = K * (x_grid[idx] / rhoORnsp)**(-L)
            ns2p_funcform_test[idx] = E2 * (x_grid[idx] / rhoORns2p)**(-F2)
            nop_funcform_test[idx] = K2 * (x_grid[idx] / rhoORnop)**(-L2)

            if x_grid[idx] < rhoORnsp:
                nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jout)**2.)
            if x_grid[idx] < rhoORns2p:
                ns2p_funcform_test[idx] = C2 * np.exp(-((x_grid[idx] - rhoRns2p)/Dout2)**2.)
            if x_grid[idx] < rhoORnop:
                nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hout2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jout2)**2.)
                
            if x_grid[idx] < rhoRnsp:
                nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)
            if x_grid[idx] < rhoRns2p:
                ns2p_funcform_test[idx] = C2 * np.exp(-((x_grid[idx] - rhoRns2p)/Din2)**2.)
            if x_grid[idx] < rhoRnop:
                nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hout2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jin2)**2.)

            if x_grid[idx] < rhocnsp:
                nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hin)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)

            if x_grid[idx] < rhocnop:
                nop_funcform_test[idx] = G2 * np.exp(-((x_grid[idx] - rhocnop)/Hin2)**2.) + I2 * np.exp(-((x_grid[idx] - rhoRnop)/Jin2)**2.)

            nec_funcform_test[idx] = (nsp_funcform_test[idx] + 2. * ns2p_funcform_test[idx] + nop_funcform_test[idx]) / totchmix_for_apo[idx]

            if x_grid[idx] <= popt[21]:
                tec[idx] = A_TE1 * (x_grid[idx] /popt[21])**B_te1
            elif x_grid[idx] > popt[22]:
                tec[idx] = A_TE2 * (x_grid[idx] / popt[22])**B_te2
            else:
                tec[idx] = A_TE1 * (x_grid[idx] / popt[21])**beta_

            if nec_funcform_test[idx] < 0.000001:
                nec_funcform_test[idx] = 0.000001
                nsp_funcform_test[idx] *= nec_funcform_test[idx]
                ns2p_funcform_test[idx] *= nec_funcform_test[idx]
                nop_funcform_test[idx] *= nec_funcform_test[idx]
                
            if (((ns2p_funcform_test[idx] / nec_funcform_test[idx]) > 0.05) and (x_grid[idx] <= 5.3)):
                flags2p[idx] = 1.

        model_opt = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
            tec, nec_funcform_test, nsp_funcform_test, ns2p_funcform_test, nop_funcform_test, max(x_grid), min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
        )

        nec_fixed_found = nec_funcform_test
        nsp_fixed_found = nsp_funcform_test
        ns2p_fixed_found = ns2p_funcform_test
        nop_fixed_found = nop_funcform_test
        tec_fixed_found = tec

        nspmix = nsp_fixed_found / nec_fixed_found

        ns2pmix = ns2p_fixed_found / nec_fixed_found

        nopmix = nop_fixed_found / nec_fixed_found


        #p0 = [nec_fixed_found, nsp_fixed_found,ns2p_fixed_found,nop_fixed_found,tec_fixed_found]
        p0 = np.concatenate([ nsp_fixed_found,ns2p_fixed_found,nop_fixed_found,tec_fixed_found])

        nec_lower_bounds = np.full(len(x_grid), 0.000001)
        nec_upper_bounds = np.full(len(x_grid), 10000.0)

        nsp_lower_bounds = np.full(len(x_grid), 0.)
        nsp_upper_bounds = np.full(len(x_grid), 2000.0)

        ns2p_lower_bounds = np.full(len(x_grid), 0.)
        ns2p_upper_bounds = np.full(len(x_grid), 2000.0)

        nop_lower_bounds = np.full(len(x_grid), 0.)
        nop_upper_bounds = np.full(len(x_grid), 2000.0)


        tec_lower_bounds = np.full(len(x_grid), 0.1)
        tec_upper_bounds = np.full(len(x_grid), 8.0)

        lower_bounds_simult = np.concatenate([ nsp_lower_bounds ,ns2p_lower_bounds ,nop_lower_bounds , tec_lower_bounds])
        upper_bounds_simult = np.concatenate([ nsp_upper_bounds, ns2p_upper_bounds, nop_upper_bounds, tec_upper_bounds])

        pold = popt
        start_time = time.time()

        result_simult = least_squares(
            residuals_function_simult, p0, jac=parallel_jacobian_simult, bounds=(lower_bounds_simult, upper_bounds_simult), max_nfev=1,
            args=(x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, nspmix, ns2pmix, nopmix, totchmix_for_apo)
        )

        popt = result_simult.x
        perr = np.sqrt(np.diag(result_simult.jac.T @ result_simult.jac))





        end_time = time.time()

        print("--- %s seconds ---" % (end_time - start_time))

        #optimized_nec = popt[:len(x_grid)]
        optimized_nsp = popt[0*len(x_grid):1*len(x_grid)]
        optimized_ns2p = popt[1 * len(x_grid):2 * len(x_grid)]
        optimized_nop = popt[2 * len(x_grid):3*len(x_grid)]
        optimized_tec = popt[3 * len(x_grid): 4 * len(x_grid)] #popt[len(x_grid):2*len(x_grid)]
        optimized_nec = (optimized_nsp + 2.*optimized_ns2p + optimized_nop)/ totchmix_for_apo


        #optimized_nec_err = perr[:len(x_grid)]

        nsp_error = perr[0*len(x_grid):1*len(x_grid)]#optimized_nsp * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (nsp_error / optimized_nsp) ** 2) 
        ns2p_error = perr[1*len(x_grid):2*len(x_grid)]#optimized_ns2p * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (ns2p_error / optimized_ns2p) ** 2) 
        nop_error = perr[2*len(x_grid):3*len(x_grid)]#optimized_nop * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (nop_error / optimized_nop) ** 2) 
        optimized_tec_err = perr[3*len(x_grid):4*len(x_grid)]#perr[len(x_grid):2*len(x_grid)]
        
        optimized_nec_err = np.sqrt((nsp_error) ** 2. + (2. * ns2p_error) ** 2.  + (nop_error) ** 2.  )/totchmix_for_apo
        
        nspmix = optimized_nsp / optimized_nec
        ns2pmix = optimized_ns2p / optimized_nec
        nopmix = optimized_nop / optimized_nec
        
        nspmix_error = nspmix*np.sqrt((nsp_error/optimized_nsp)**2. + (optimized_nec_err/optimized_nec)**2. )
        ns2pmix_error = ns2pmix*np.sqrt((ns2p_error/optimized_ns2p)**2. + (optimized_nec_err/optimized_nec)**2. )
        nopmix_error = nopmix*np.sqrt((nop_error/optimized_nop)**2. + (optimized_nec_err/optimized_nec)**2. )



        output_data = np.column_stack((optimized_tec, optimized_tec_err, optimized_nec, optimized_nec_err, optimized_nsp, nsp_error, optimized_ns2p, ns2p_error, optimized_nop, nop_error, nspmix, ns2pmix, nopmix, nspmix_error, ns2pmix_error, nopmix_error))
        np.savetxt( nightname_strings[w] + '_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dusk_funcform_tec3partpl_python_60x16.csv', output_data, delimiter=',', comments='')
        np.savetxt( nightname_strings[w] + '_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_dusk_funcform_params_pold_tec3partpl_python_23vals.csv', pold, delimiter=',', comments='')
        plt.clf()
        fig, axs = plt.subplots(4, 2, figsize=(10, 8),sharex=True)

        axs[0,0].errorbar(x_grid, optimized_tec, yerr=optimized_tec_err, fmt='o', label='Tec')
        axs[0,0].set_xlabel(XTITTLEE)
        axs[0,0].set_ylabel('T$_{ec}$ (eV)')
        axs[0,0].set_ylim([0,10])
        axs[0,0].legend()
        axs[0,0].grid(True)

        axs[0,1].errorbar(x_grid, optimized_nec, yerr=optimized_nec_err, fmt='o', label='nec')
        axs[0,1].set_xlabel(XTITTLEE)
        axs[0,1].set_ylabel('n$_{ec}$ ($cm^{-3}$)')

        axs[0,1].set_ylim([0,4000])
        axs[0,1].legend()
        axs[0,1].grid(True)

        axs[1,0].errorbar(x_grid, optimized_nsp, yerr=nsp_error, fmt='o', label='nsp')
        axs[1,0].set_xlabel(XTITTLEE)
        axs[1,0].set_ylabel('n$_{S^{+}}$ ($cm^{-3}$)')
        axs[1,0].set_ylim([0,800])
        axs[1,0].legend()
        axs[1,0].grid(True)

        axs[1,1].errorbar(x_grid, optimized_nsp/optimized_nec, yerr=nspmix_error, fmt='o', label='nspmix')
        axs[1,1].set_xlabel(XTITTLEE)
        axs[1,1].set_ylabel('n$_{S^{+}}$ Mixing Ratio')

        axs[1,1].set_ylim([0,0.8])
        axs[1,1].legend()
        axs[1,1].grid(True)

        axs[2,0].errorbar(x_grid, optimized_ns2p, yerr=ns2p_error, fmt='o', label='ns2p')
        axs[2,0].set_xlabel(XTITTLEE)
        axs[2,0].set_ylabel('n$_{S^{++}}$ ($cm^{-3}$)')
        axs[2,0].set_ylim([0,800])
        axs[2,0].legend()
        axs[2,0].grid(True)

        axs[2,1].errorbar(x_grid, optimized_ns2p/optimized_nec, yerr=ns2pmix_error, fmt='o', label='ns2pmix')
        axs[2,1].set_xlabel(XTITTLEE)
        axs[2,1].set_ylabel('n$_{S^{++}}$ Mixing Ratio')
        axs[2,1].set_ylim([0,0.4])
        axs[2,1].legend()
        axs[2,1].grid(True)

        axs[3,0].errorbar(x_grid, optimized_nop, yerr=nop_error, fmt='o', label='nop')
        axs[3,0].set_xlabel(XTITTLEE)
        axs[3,0].set_ylabel('n$_{O^{+}}$ ($cm^{-3}$)')
        axs[3,0].set_ylim([0,800])
        axs[3,0].legend()
        axs[3,0].grid(True)


        axs[3,1].errorbar(x_grid, optimized_nop/optimized_nec, yerr=nopmix_error, fmt='o', label='nopmix')
        axs[3,1].set_xlabel(XTITTLEE)
        axs[3,1].set_ylabel('n$_{O^{+}}$ Mixing Ratio')
        axs[3,1].set_ylim([0,0.6])
        axs[3,1].legend()
        axs[3,1].grid(True)

        plt.tight_layout()
        plt.savefig( nightname_strings[w] + '_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_bestfit_params_python_dusk_funcform_tec3partpl.png')
        #os.utime('diff_ev_bestfit_params_python_dusk_funcform_tec3partpl.png', (time.time(), time.time()))  # Update timest
        plt.show()
        model_spectra = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
            optimized_tec, optimized_nec, optimized_nec * nspmix, optimized_nec * ns2pmix, optimized_nec * nopmix,
            np.max(x_grid), np.min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
        )

        np.savetxt( nightname_strings[w] + '_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_model_dusk_funcform_tec3partpl_python_60x8.csv', model_spectra, delimiter=',', comments='')
        plt.clf()
        fig, axs = plt.subplots(8, 1, figsize=(10, 20), sharex=True)

        for i, ax in enumerate(axs):
            ax.errorbar(x_grid, dawn_grid[:, i], yerr=err_dawn_grid[:, i], fmt='o', label=f'Line {i} Data')
            ax.plot(x_grid, model_spectra[:, i], label=f'Line {i} Fit')
            ax.set_ylabel(f'Rayleighs (Line {i})')
            ax.legend()
            ax.grid(True)

        axs[-1].set_xlabel(XTITTLEE)
        plt.tight_layout()
        #plt.savefig('newdiff_ev_bestfit_emiss_for_python_dusk_funcform_tec3partpl.png')
        #os.utime('diff_ev_bestfit_emiss_for_python_dusk_funcform_tec3partpl.png', (time.time(), time.time()))  # Update timest
        plt.show()
        
        colors = ['pink', 'purple', 'blue', 'lightgreen', 'darkgreen', 'orange', 'red', 'brown']
        labels = [r'$S^{++}$ 3722 Å', r'$O^{+}$ 3726 Å', r'$O^{+}$ 3729 Å', r'$S^{+}$ 4069 Å', r'$S^{+}$ 4076 Å', r'$S^{++}$ 6312 Å', r'$S^{+}$ 6716 Å', r'$S^{+}$ 6731 Å']
        plt.figure(figsize=( 6.67,10))
        idxs = [7,6,5,4,3,2,1,0]
        # Plot each simulated profile
        for i in idxs:
            plt.errorbar(x_grid, dawn_grid[:, i], yerr=err_dawn_grid[:, i], color=colors[i], label=labels[i])
            plt.plot(x_grid, model_spectra[:, i], color=colors[i])
        plt.xlabel(r'Dusk $\rho_c$ ($R_J$)')
        plt.ylabel('Rayleighs')
        plt.ylim([0,400])
        plt.legend(loc='best')
        plt.savefig( nightname_strings[w] + '_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_bestfit_emiss_for_python_dusk_funcform_tec3partpl.png')
        plt.show()
        
   

