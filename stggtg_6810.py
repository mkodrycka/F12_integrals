"""
stggtg.py

Hard-coded fit parameters for {6,8,10}-GTG to STG fit.
"""
import psi4
import numpy as np

def stggtg(gamma, f12sq_primitive = False, cgtg = 6, eventempered=False):
    """
    Returns a psi4 Vector of coefficients and exponents for the fit.
    Parameters generated using Molpro using default weight function.
    """

    if cgtg == 6:
        coeff = np.zeros((6))
        exp  = np.zeros((6))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] =  0.270700;  exp[0] =   0.195320;
            coeff[1] =  0.305520;  exp[1] =   0.819200;
            coeff[2] =  0.182970;  exp[2] =   2.859170;
            coeff[3] =  0.109860;  exp[3] =   9.500730;
            coeff[4] =  0.068100;  exp[4] =  35.699890;
            coeff[5] =  0.042240;  exp[5] = 197.793280;
        else:
            print("We don't have coeff and exp for this gamma")

    if cgtg == 8:
        coeff = np.zeros((8))
        exp  = np.zeros((8))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] =  0.02125;  exp[0] =   0.02000;
            coeff[1] =  0.31307;  exp[1] =   0.25910;
            coeff[2] =  0.27338;  exp[2] =   1.01383;
            coeff[3] =  0.16214;  exp[3] =   3.26236;
            coeff[4] =  0.09900;  exp[4] =   10.07154;
            coeff[5] =  0.05711;  exp[5] =   31.60266;
            coeff[6] =  0.03033;  exp[6] =   99.01650;
            coeff[7] =  0.02525;  exp[7] =   310.25244;
        else:
            print("We don't have coeff and exp for this gamma")

    if cgtg == 10:
        coeff = np.zeros((10))
        exp  = np.zeros((10))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] =  -0.11505;  exp[0] =  0.02000;
            coeff[1] =  0.22241;  exp[1] =   0.04515;
            coeff[2] =  0.39087;  exp[2] =   0.48862;
            coeff[3] =  0.21771;  exp[3] =   2.05674;
            coeff[4] =  0.09624;  exp[4] =   5.70016;
            coeff[5] =  0.08163;  exp[5] =   15.16035;
            coeff[6] =  0.03626;  exp[6] =   42.93890;
            coeff[7] =  0.03583;  exp[7] =   120.84800;
            coeff[8] =  0.00438;  exp[8] =   331.71342;
            coeff[9] =  0.02058;  exp[9] =   912.32744;
            
            """
            coeff[0] = -0.10993;  exp[0] = 0.02000;
            coeff[1] = 0.21644;  exp[1] =  0.04536;
            coeff[2] = 0.39119;  exp[2] =  0.48706;
            coeff[3] = 0.21836;  exp[3] =  2.05420;
            coeff[4] = 0.09634;  exp[4] =  5.71783;
            coeff[5] = 0.08123;  exp[5] =  15.15629;
            coeff[6] = 0.03654;  exp[6] =  42.93554;
            coeff[7] = 0.03557;  exp[7] =  120.84873;
            coeff[8] = 0.00463;  exp[8] =  331.71330;
            coeff[9] = 0.02043;  exp[9] =  912.32751;
	    """

        else:
            print("We don't have coeff and exp for this gamma")

    if (eventempered == True) and (cgtg == 6):
        coeff = np.zeros((6))
        exp  = np.zeros((6))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] =  0.024095;  exp[0] = 0.064150;
            coeff[1] =  0.200202;  exp[1] = 0.192450;
            coeff[2] =  0.241131;  exp[2] = 0.577350;
            coeff[3] =  0.257333;  exp[3] = 1.732051;
            coeff[4] =  0.030452;  exp[4] = 5.196152;
            coeff[5] =  0.200003;  exp[5] = 15.588457;
        else:
            print("We don't have coeff and exp for this gamma")

    if (eventempered == True) and (cgtg == 8):
        coeff = np.zeros((8))
        exp  = np.zeros((8))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] = -0.000836;  exp[0] =   0.021383;
            coeff[1] =  0.028514;  exp[1] =   0.064150;
            coeff[2] =  0.186799;  exp[2] =   0.192450;
            coeff[3] =  0.273799;  exp[3] =   0.577350;
            coeff[4] =  0.186594;  exp[4] =   1.732051;
            coeff[5] =  0.163105;  exp[5] =   5.196152;
            coeff[6] =  0.018372;  exp[6] =   15.588457;
            coeff[7] =  0.116655;  exp[7] =   46.765372;
        else:
            print("We don't have coeff and exp for this gamma")

    if (eventempered == True) and (cgtg == 10):
        coeff = np.zeros((10))
        exp  = np.zeros((10))

        if abs(gamma - 1.0) < 1E-6:
            coeff[0] = 0.000149 ;  exp[0] =  0.007128;
            coeff[1] = -0.000984 ; exp[1] =  0.021383;
            coeff[2] = 0.027772 ;  exp[2] =  0.064150;
            coeff[3] = 0.189959 ;  exp[3] =  0.192450;
            coeff[4] = 0.265391 ;  exp[4] =  0.577350;
            coeff[5] = 0.205997 ;  exp[5] =  1.732051;
            coeff[6] = 0.121946 ;  exp[6] =  5.196152;
            coeff[7] = 0.094988 ;  exp[7] =  15.588457;
            coeff[8] = 0.012027 ;  exp[8] =  46.765372;
            coeff[9] = 0.067137 ;  exp[9] =  140.296115;
        else:
            print("We don't have coeff and exp for this gamma")

    cgtg_params = []
    if f12sq_primitive == True:
        for e_i, c_i in zip(exp, coeff):
            for e_j, c_j in zip(exp, coeff):
                cgtg_params.append((e_i + e_j, c_i*c_j))
    else: 	
    	for e, c in zip(exp, coeff):
    	    cgtg_params.append((e, c)) 	
    
    return cgtg_params	
