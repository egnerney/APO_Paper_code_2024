#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 17:15:37 2024

@author: edne8319
"""

import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.optimize import differential_evolution, least_squares
from scipy.integrate import newton_cotes
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import time

def newton_cotes_integration(x, y, num_points=5):
    """
    Perform Newton-Cotes integration on tabulated data using a specified number of points.

    Parameters:
    x : array-like
        The x-values of the data points.
    y : array-like
        The y-values of the data points corresponding to x.
    num_points : int
        The number of points to use for the Newton-Cotes formula (default is 5).

    Returns:
    integral : float
        The estimated integral of the function.
    """
    n = len(x)
    if n < num_points:
        raise ValueError(f"At least {num_points} points are required for Newton-Cotes integration.")

    weights, _ = newton_cotes(num_points - 1)
    integral = 0.0
    h = (x[1] - x[0])  # Assuming uniform spacing

    for i in range(0, n - num_points + 1, num_points - 1):
        integral += h * np.dot(weights, y[i:i+num_points])

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
                rayleighs[idx, linenum] = conversion * newton_cotes_integration(s_LOS[idx_want], emiss_s) #simps(emiss_s, s_LOS[idx_want])
        elif len(idx_want) == 1:
            path_length = s_LOS[idx_want[0] + 1] - s_LOS[idx_want[0]] if idx_want[0] + 1 < len(s_LOS) else drho
            for linenum in range(8):
                emiss_value = interp1d(x_LOS, epsilon_a[:, linenum], kind='linear')(rho_LOS[idx, idx_want])[0]
                rayleighs[idx, linenum] = conversion * emiss_value * path_length
        else:
            rayleighs[idx, :] = 0

    return rayleighs



def residuals_function(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl):
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
    
    flagsp = np.zeros(nx_grid)
    
    rhocne = p[0]#;5.27
    rhoRne =  rhoRne_# ;p[1];5.66d;p[1];5.69
    rhoORne = rhoORne_#;5.82d;p[1];5.8 rhoRne_,rhoORne_,rhocne_
    A = p[11]#;1960
    Bin = p[2]#;p[12];0.2
    Bout = p[12]#;0.2
    C = p[13]#;3430
    Din = p[14]#;0.2
    Dout = p[15]#;0.2
    F = nec_pl#;5.4d;p[14];5.4d

    E = A*np.exp(-((rhoORne - rhocne)/Bout) ** 2.) + C*np.exp(-((rhoORne - rhoRne)/Dout) **2.)

    rhocnsp = p[0]
    rhoRnsp = rhoRnsp_
    rhoORnsp = rhoORne_
    G = p[1]
    Hin = p[2]
    Hout = p[3]
    I = p[4]
    Jin = p[5]
    Jout = p[6]
    L = 3.37 + nec_pl  # assuming nec_pl is p[14]

    K = G * np.exp(-((rhoORnsp - rhocnsp)/Hout)**2.) + I * np.exp(-((rhoORnsp - rhoRnsp)/Jout)**2.)

   
    tec = np.zeros(nx_grid)
    A_TE1 = p[7]
    B_te1 = p[8]
    B_te2 = p[9]
    A_TE2 = p[10]

    beta_ = np.log(A_TE2 / A_TE1) / np.log(rhoORne_ / rho4069peak)

    for idx in range(nx_grid):
        nec_funcform_test[idx] = E*(x_grid[idx]/rhoORne)**(-F)
        nsp_funcform_test[idx] = K * (x_grid[idx] / rhoORnsp)**(-L)
        

        if x_grid[idx] < rhoORnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jout)**2.)
        
        if x_grid[idx] < rhoRnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)
        
        if x_grid[idx] < rhocnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hin)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)

        if x_grid[idx] < rhoORne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bout)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Dout)**2.)
        
        if x_grid[idx] < rhoRne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bout)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Din)**2.)
         
        if x_grid[idx] < rhocne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bin)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Din)**2.)

       # nec_funcform_test[idx] = (nsp_funcform_test[idx] + 2. * ns2p_funcform_test[idx] + nop_funcform_test[idx]) / totchmix_for_apo[idx]

        if x_grid[idx] <= rho4069peak:
            tec[idx] = A_TE1 * (x_grid[idx] / rho4069peak)**B_te1
        elif x_grid[idx] > rhoORne_:
            tec[idx] = A_TE2 * (x_grid[idx] / rhoORne_)**B_te2
        else:
            tec[idx] = A_TE1 * (x_grid[idx] / rho4069peak)**beta_

        if nec_funcform_test[idx] < 0.000001:
            nspmixtemp = nsp_funcform_test[idx] / nec_funcform_test[idx]
            #ns2pmixtemp = ns2p_funcform_test[idx] / nec_funcform_test[idx]
            #nopmixtemp = nop_funcform_test[idx] / nec_funcform_test[idx]
            nec_funcform_test[idx] = 0.0000011
            nsp_funcform_test[idx] =  nspmixtemp * nec_funcform_test[idx]
            #ns2p_funcform_test[idx] = ns2pmixtemp *  nec_funcform_test[idx]
            #nop_funcform_test[idx] = nopmixtemp *  nec_funcform_test[idx]
            
        if ((nsp_funcform_test[idx] / nec_funcform_test[idx]) > 0.9) :
            flagsp[idx] = 1.

    model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
        tec, nec_funcform_test, nsp_funcform_test, ns2p_funcform_test, nop_funcform_test, max(x_grid), min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
    )

    idx_want = [3,4,6,7]# [0, 1, 2, 3, 5, 6, 7]
    model = model[:, idx_want].flatten()

    penalt = 1e6
    flagtotsp = np.sum(flagsp)

    if flagtotsp > 0:
        model = model +  penalt * flagtotsp / len(model)

    result = np.sum(((y - model) / err) ** 2.)# (y - model) / err
    return result


def residuals_function_simult(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d):
    
    nx_grid = len(x_grid)
    nec_tempi = p[0:nx_grid]
    nsp_tempi = p[nx_grid:2 * nx_grid]
    ns2p_tempi = np.zeros(nx_grid)#p[2 * nx_grid:3 * nx_grid]
    nop_tempi = np.zeros(nx_grid)# p[3 * nx_grid:4 * nx_grid]
    tec_tempi = p[2 * nx_grid:3 * nx_grid]

    model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
        tec_tempi, nec_tempi, nsp_tempi, ns2p_tempi, nop_tempi, max(x_grid), min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
    )
    idx_want = [3,4,6,7]#[0, 1, 2, 3, 5, 6, 7]
    model = model[:, idx_want].flatten()


    result = (y - model) / err
    
    return result


def parallel_jacobian(p, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl, n_jobs=-1):
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
        f1 = residuals_function(p1, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl)
        f2 = residuals_function(p2, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl)
        return (f1 - f2) / (2. * eps)

    jacobian_columns = Parallel(n_jobs=n_jobs)(delayed(compute_jacobian_column)(i) for i in range(len(p)))
    jacobian = np.column_stack(jacobian_columns)
    return jacobian



def parallel_jacobian_simult(p,  x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, n_jobs=-1):
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
        f1 = residuals_function_simult(p1, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d)
        f2 = residuals_function_simult(p2, x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d)
        return (f1 - f2) / (2. * eps)

    jacobian_columns = Parallel(n_jobs=n_jobs)(delayed(compute_jacobian_column_simult)(i) for i in range(len(p)))
    jacobian = np.column_stack(jacobian_columns)
    return jacobian




if __name__ == '__main__':

    # Load and preprocess data
    yptsi_in = np.loadtxt('yptsi_8_APO_LINES_optical_vary_nec_tec_lookuptable_66x58x8_CHIANTI_10.1.csv', delimiter=',')
    yptsi_in = yptsi_in.reshape((8, 58, 66)).transpose(2, 1, 0)

    x_grid = np.loadtxt('x_grid_to_plot_60_points.csv', delimiter=',')
    dawn_grid1 = np.loadtxt('dusk_grid_to_plot_60x8.csv', delimiter=',')
    dawn_grid = np.transpose(dawn_grid1)

    err_dawn_grid1 = np.loadtxt('err_dusk_grid_to_plot_60x8.csv', delimiter=',')
    err_dawn_grid = np.transpose(err_dawn_grid1)

    nel = np.concatenate(([0.000001, 0.00001, 0.0001, 0.001, 0.1], np.arange(9.) + 1., 10. * np.arange(9.) + 10., 100. * np.arange(4.) + 100., 250. * np.arange(39.) + 500.))
    Tec = np.concatenate(([0.01], 0.1 * np.arange(9.) + 0.1, 0.5 * np.arange(17.) + 1., np.arange(10.) + 10., 5. * np.arange(8.) + 20., 20. * np.arange(7.) + 60., 80. * np.arange(6.) + 200.))

    n_ne = len(nel)
    n_te = len(Tec)
    nel2dto1d = np.repeat(nel, n_te)
    tec2dto1d = np.tile(Tec, n_ne)
    yptsi_in_2dto1d = yptsi_in.reshape((n_ne * n_te, 8))

    XTITTLEE = 'Dusk $\\rho_c$ ($R_J$)'
    idx_want = [3,4,6,7]#[0, 1, 2, 3, 5, 6, 7]
    dawn_grid_temp = dawn_grid[:, idx_want]
    err_dawn_grid_temp = err_dawn_grid[:, idx_want]

    ndawn_grid = dawn_grid_temp.size
    y = dawn_grid_temp.flatten()
    err = err_dawn_grid_temp.flatten()
    x = np.ones(ndawn_grid)

    x_step_size = round(abs(x_grid[1] - x_grid[0]), 5)
    x_grid_new = round(abs(x_grid[1] - x_grid[0]), 5) * np.arange(len(x_grid)) + min(x_grid)
    x_grid = x_grid_new[::-1]

    out_fits1 = np.loadtxt('outs_and_errors_dusk_funcform_fit_60x16.csv', delimiter=',')
    out_fits = np.transpose(out_fits1)

    totchmix_for_apo = np.loadtxt('totchmix_for_apo_7pointavg_dusk_unshifted.csv', delimiter=',')

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
   # plt.savefig('initial_emiss_funcform_for_python_dusk.png')
    plt.show()
    nec_pl = 3.

    rhocne = 5.057
    rhoRne = 5.7

    rhoRns2p_ =  5.7 # tried 5.75d before Dusk
    Din2max = 0.2

    rhoORne = 5.82 #d ;Dusk tried 5.94 too thought og 5.82 better
    A = 1200.#d;1871d;1059d ; 1792d
    Bin = 0.2275#d;p_out[2]
    Bout = 0.239#d;p_out[3]
    C = 3120.#d;2994d;2783.55;3430d
    Din = 0.244#d;p_out[12]
    Dout = 0.19#d; p_out[13]
    F = nec_pl
    E = A*np.exp(-((rhoORne - rhocne)/Bout) ** 2.) + C*np.exp(-((rhoORne - rhoRne)/Dout) ** 2.)



    

    rhocnsp =    rhocne
    rhoRnsp = 5.66
    rhoORnsp = rhoORne
    G = 650.
    Hin= 0.2275
    Hout= 0.33
    I = 600.
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


   # lower_bounds = [4.9,300., 0.1, 0.1, 100., 0.05, 0.05, 100.,0.05,0.05,100.,100.,0.05,0.05, 1., 0.01, 0.01, 3. ]
    #upper_bounds = [5.1,1400., 0.5,0.6, 1000.,0.5, 1.5, 1000.,0.8,2.5,1000.,1000.,0.5,2.5,   4.,10., 10.,8.   ]
    
    lower_bounds = [4.9,300., 0.1, 0.1, 100., 0.05, 0.05,1., 0.01,0.01,3.,400.,0.1,1000.,0.1,0.1 ]
    upper_bounds = [5.1,1500., 0.5,0.6, 1000.,0.5, 1.5,4.,10.,10.,8.,2000.,0.6,5000.,0.6,0.6   ]
    # Initialize parameters (these should be aligned with the p[] indices used in the residuals function)
    #p0 = [5.27, 1298, 0.2, 0.3, 743, 0.08, 0.1, 520, 0.09, 0.91, 470, 1106, 0.1, 0.2, 3.37, 1.0, 2.0, 3.0]

    #p0 = [rhocnsp, G,Hin,Hout,I,Jin,Jout,A_TE1,B_te1,B_te2, A_TE2,A,Bout,C,Din,Dout]
    #p0 = [5.06, 625.,0.251,Hout,600.,Jin,Jout,650.,Din2,Dout2, 350., I2,Jin2,Jout2, A_TE1,B_te1 ,B_te2, A_TE2]
    p0 = [rhocnsp, G,Hin,Hout,I,Jin,Jout,A_TE1,B_te1,B_te2, A_TE2,A,Bout,C,Din,Dout]
   
    
    start_time = time.time()

    bounds = [(4.9, 5.1), (300., 1400.), (0.1, 0.5), (0.1, 0.6), (100., 1000.), (0.05, 0.5), (0.05, 1.5), (1., 4.), (0.01, 10.), (0.01, 10.), (3., 8.),(800.,2000.),(0.1,0.6),(2000.,4000.),(0.1,0.6),(0.1,0.6)]




    result = differential_evolution(
            residuals_function,
            bounds=bounds,updating='deferred',
            args=(x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d, rhoRne_, rhoORne_, rhocne_, rhoRnsp_, rho4069peak, nec_pl),
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
    
    flagsp = np.zeros(nx_grid)
    
    rhocne = popt[0]#;5.27
    rhoRne =  rhoRne_# ;p[1];5.66d;p[1];5.69
    rhoORne = rhoORne_#;5.82d;p[1];5.8 rhoRne_,rhoORne_,rhocne_
    A = popt[11]#;1960
    Bin = popt[2]#;p[12];0.2
    Bout = popt[12]#;0.2
    C = popt[13]#;3430
    Din = popt[14]#;0.2
    Dout = popt[15]#;0.2
    F = nec_pl#;5.4d;p[14];5.4d

    E = A*np.exp(-((rhoORne - rhocne)/Bout) ** 2.) + C*np.exp(-((rhoORne - rhoRne)/Dout) **2.)

    rhocnsp = popt[0]
    rhoRnsp = rhoRnsp_
    rhoORnsp = rhoORne_
    G = popt[1]
    Hin = popt[2]
    Hout = popt[3]
    I = popt[4]
    Jin = popt[5]
    Jout = popt[6]
    L = 3.37 + nec_pl  # assuming nec_pl is p[14]

    K = G * np.exp(-((rhoORnsp - rhocnsp)/Hout)**2.) + I * np.exp(-((rhoORnsp - rhoRnsp)/Jout)**2.)

   
    tec = np.zeros(nx_grid)
    A_TE1 = popt[7]
    B_te1 = popt[8]
    B_te2 = popt[9]
    A_TE2 = popt[10]

    beta_ = np.log(A_TE2 / A_TE1) / np.log(rhoORne_ / rho4069peak)

    for idx in range(nx_grid):
        nec_funcform_test[idx] = E*(x_grid[idx]/rhoORne)**(-F)
        nsp_funcform_test[idx] = K * (x_grid[idx] / rhoORnsp)**(-L)
        

        if x_grid[idx] < rhoORnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jout)**2.)
        
        if x_grid[idx] < rhoRnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hout)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)
        
        if x_grid[idx] < rhocnsp:
            nsp_funcform_test[idx] = G * np.exp(-((x_grid[idx] - rhocnsp)/Hin)**2.) + I * np.exp(-((x_grid[idx] - rhoRnsp)/Jin)**2.)

        if x_grid[idx] < rhoORne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bout)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Dout)**2.)
        
        if x_grid[idx] < rhoRne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bout)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Din)**2.)
         
        if x_grid[idx] < rhocne:
            nec_funcform_test[idx] = A * np.exp(-((x_grid[idx] - rhocne)/Bin)**2.) + C * np.exp(-((x_grid[idx] - rhoRne)/Din)**2.)

       # nec_funcform_test[idx] = (nsp_funcform_test[idx] + 2. * ns2p_funcform_test[idx] + nop_funcform_test[idx]) / totchmix_for_apo[idx]

        if x_grid[idx] <= rho4069peak:
            tec[idx] = A_TE1 * (x_grid[idx] / rho4069peak)**B_te1
        elif x_grid[idx] > rhoORne_:
            tec[idx] = A_TE2 * (x_grid[idx] / rhoORne_)**B_te2
        else:
            tec[idx] = A_TE1 * (x_grid[idx] / rho4069peak)**beta_

        if nec_funcform_test[idx] < 0.000001:
            nspmixtemp = nsp_funcform_test[idx] / nec_funcform_test[idx]
            #ns2pmixtemp = ns2p_funcform_test[idx] / nec_funcform_test[idx]
            #nopmixtemp = nop_funcform_test[idx] / nec_funcform_test[idx]
            nec_funcform_test[idx] = 0.0000011
            nsp_funcform_test[idx] =  nspmixtemp * nec_funcform_test[idx]
            #ns2p_funcform_test[idx] = ns2pmixtemp *  nec_funcform_test[idx]
            #nop_funcform_test[idx] = nopmixtemp *  nec_funcform_test[idx]
            

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
    p0 = np.concatenate([nec_fixed_found, nsp_fixed_found,tec_fixed_found])

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

    #lower_bounds_simult = np.concatenate([nec_lower_bounds, nsp_lower_bounds ,ns2p_lower_bounds ,nop_lower_bounds , tec_lower_bounds])
    #upper_bounds_simult = np.concatenate([nec_upper_bounds, nsp_upper_bounds, ns2p_upper_bounds, nop_upper_bounds, tec_upper_bounds])

    lower_bounds_simult = np.concatenate([nec_lower_bounds, nsp_lower_bounds , tec_lower_bounds])
    upper_bounds_simult = np.concatenate([nec_upper_bounds, nsp_upper_bounds, tec_upper_bounds])


    start_time = time.time()
    
    idx_want = [3,4,6,7]#[0, 1, 2, 3, 5, 6, 7]
    dawn_grid_temp = dawn_grid[:, idx_want]
    err_dawn_grid_temp = err_dawn_grid[:, idx_want]

    ndawn_grid = dawn_grid_temp.size
    y = dawn_grid_temp.flatten()
    err = err_dawn_grid_temp.flatten()
    x = np.ones(ndawn_grid)

    result_simult = least_squares(
        residuals_function_simult, p0, jac=parallel_jacobian_simult, bounds=(lower_bounds_simult, upper_bounds_simult), max_nfev=1,
        args=(x, y, err, x_grid, x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d)
    )

    popt = result_simult.x
    perr = np.sqrt(np.diag(result_simult.jac.T @ result_simult.jac))





    end_time = time.time()

    print("--- %s seconds ---" % (end_time - start_time))

    optimized_nec = popt[:len(x_grid)]
    optimized_nsp = popt[len(x_grid):2*len(x_grid)]
    optimized_ns2p = np.zeros(len(x_grid)) #popt[2 * len(x_grid):3 * len(x_grid)]
    optimized_nop = np.zeros(len(x_grid))#popt[3 * len(x_grid):4*len(x_grid)]
    optimized_tec = popt[2 * len(x_grid): 3 * len(x_grid)] #popt[len(x_grid):2*len(x_grid)]



    optimized_nec_err = perr[:len(x_grid)]

    nsp_error = perr[len(x_grid):2*len(x_grid)]#optimized_nsp * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (nsp_error / optimized_nsp) ** 2) 
    ns2p_error = np.zeros(len(x_grid))#perr[2*len(x_grid):3*len(x_grid)]#optimized_ns2p * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (ns2p_error / optimized_ns2p) ** 2) 
    nop_error = np.zeros(len(x_grid))#perr[3*len(x_grid):4*len(x_grid)]#optimized_nop * np.sqrt((optimized_nec_err / optimized_nec) ** 2 + (nop_error / optimized_nop) ** 2) 
    optimized_tec_err = perr[2*len(x_grid):3*len(x_grid)]#perr[len(x_grid):2*len(x_grid)]


    nspmix_error = ( optimized_nsp/optimized_nec) * np.sqrt((nsp_error/optimized_nsp) ** 2. + (optimized_nec_err/optimized_nec) ** 2. )
   
    output_data = np.column_stack((optimized_tec, optimized_tec_err, optimized_nec, optimized_nec_err, optimized_nsp, nsp_error,nspmix, nspmix_error))
    np.savetxt('all4lines_onlysp_nonec_ratio_diff_ev_outs_and_errors_dusk_funcform_tec3partpl_python_60x16.csv', output_data, delimiter=',', comments='')

    fig, axs = plt.subplots(4, 1, figsize=(10, 8))

    axs[0].errorbar(x_grid, optimized_tec, yerr=optimized_tec_err, fmt='o', label='Tec')
    axs[0].set_xlabel(XTITTLEE)
    axs[0].set_ylabel('T$_{ec}$ (eV)')
    axs[0].set_title('Optimized TEC Radial Profile')
    axs[0].set_ylim([0,10])
    axs[0].legend()
    axs[0].grid(True)

    axs[1].errorbar(x_grid, optimized_nec, yerr=optimized_nec_err, fmt='o', label='nec')
    axs[1].set_xlabel(XTITTLEE)
    axs[1].set_ylabel('n$_{ec}$ ($cm^{-3}$)')
    axs[1].set_title('Optimized nEC Radial Profile')
    axs[1].set_ylim([0,4000])
    axs[1].legend()
    axs[1].grid(True)

    axs[2].errorbar(x_grid, optimized_nsp, yerr=nsp_error, fmt='o', label='nsp')
    axs[2].set_xlabel(XTITTLEE)
    axs[2].set_ylabel('n$_{S^{+}}$ ($cm^{-3}$)')
    axs[2].set_title(r'Optimized n$_{S^{+}}$ Radial Profile')
    axs[2].set_ylim([0,1000])
    axs[2].legend()
    axs[2].grid(True)
    
    axs[3].errorbar(x_grid, optimized_nsp/optimized_nec, yerr=nspmix_error, fmt='o', label='nspmix')
    axs[3].set_xlabel(XTITTLEE)
    axs[3].set_ylabel('$S^{+}$ Mixing Ratio')
    axs[3].set_title(r'Optimized n$_{S^{+}}$ Radial Profile')
    axs[3].set_ylim([0,0.8])
    axs[3].legend()
    axs[3].grid(True)

   

    plt.tight_layout()
    plt.savefig('all4lines_onlysp_nonec_ratio_diff_ev_bestfit_params_python_dusk_funcform_tec3partpl.png')
    plt.show()

    model_spectra = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(
        optimized_tec, optimized_nec, optimized_nec * nspmix, optimized_nec * ns2pmix, optimized_nec * nopmix,
        np.max(x_grid), np.min(x_grid), x_step_size, yptsi_in_2dto1d, nel2dto1d, tec2dto1d
    )

    np.savetxt('all4lines_onlysp_nonec_ratio_diff_ev_model_dusk_funcform_tec3partpl_python_60x8.csv', model_spectra, delimiter=',', comments='')

    fig, axs = plt.subplots(8, 1, figsize=(10, 20), sharex=True)

    for i, ax in enumerate(axs):
        ax.errorbar(x_grid, dawn_grid[:, i], yerr=err_dawn_grid[:, i], fmt='o', label=f'Line {i} Data')
        ax.plot(x_grid, model_spectra[:, i], label=f'Line {i} Fit')
        ax.set_ylabel(f'Rayleighs (Line {i})')
        ax.legend()
        ax.grid(True)

    axs[-1].set_xlabel(XTITTLEE)
    plt.tight_layout()
    plt.savefig('all4lines_onlysp_nonec_ratio_diff_ev_bestfit_emiss_for_python_dusk_funcform_tec3partpl.png')
    plt.show()

