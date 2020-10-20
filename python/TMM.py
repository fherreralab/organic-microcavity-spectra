from dielectric_functions import sellmeier, lorentzAg, lorentzTDBC 
import numpy as np     
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import floor

def DrawCavity(cavity,fons=12,figsize=(7,7)):
    """
    cavity: Cavity list of dictionaries
    fons: Label fontsize
    figsize: Canvas size

    Outputs cavity diagram plot
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)

    cavity = cavity[::-1]
    for i,layer in enumerate(cavity):
        th = 0.1
        posText = (0.31,i*th*1.02+0.1*th)
        posLayer = (0.3, i*0.1)
        fons = 12
        if i == len(cavity)-1 and layer['name'] == 'substrate':
            rect = plt.Rectangle((0.2, i*0.1), 0.6, th, color='m', alpha=0.3)
            plt.text(*posText,f"{layer['name']} {layer['thickness']} [nm]", fontsize=fons)
        elif i == 0:
            rect = plt.Rectangle((0.2, i*0.1), 0.6, th, color='k', alpha=0.1)
            plt.text(*posText,f"{layer['name']} {layer['thickness']} [nm]", fontsize=fons)
        elif layer['name'] == 'silica':
            rect = plt.Rectangle(posLayer, 0.4, th, color='b', alpha=0.3)
            plt.text(*posText,f"{layer['name']} {layer['thickness']} [nm]", fontsize=fons)
        elif layer['name'] == 'TDBC':
            rect = plt.Rectangle(posLayer, 0.4, th, color='g', alpha=0.3)
            plt.text(*posText,f"{layer['name']} {layer['thickness']} [nm]", fontsize=fons)
        else:
            rect = plt.Rectangle(posLayer, 0.4, th, color='k', alpha=0.4)
            plt.text(*posText,f"{layer['name']} {layer['thickness']} [nm]", fontsize=fons)
        ax.add_patch(rect) 

    for i in range(4):
        arr = plt.Arrow(posLayer[0]+i*0.13,posLayer[1]*1.4,0,-posLayer[1]*0.2,width=0.1,color='y',alpha=0.6)
        ax.add_patch(arr)
    plt.axis('off')
    plt.axis('scaled')
    plt.show()

def refractiveIndices(cavity,wavelength):
    """
    cavity: Cavity list of dictionaries
    wavelength: Wavelength in nanometers

    Outputs a list of complex refractive indices (if the material has attenuation) for the given wavelength.
    """

    N = []
    nIn = 1.4587
    nSubs = nIn
    nOut = 1
    for layer in cavity:
        if layer['name'] == 'substrate':
            N.append(nSubs)
        elif layer['name'] == 'silver':
            N.append(np.sqrt(lorentzAg(wavelength)))
        elif layer['name'] == 'silica':
            N.append(np.sqrt(sellmeier(wavelength)))
        elif layer['name'] == 'TDBC':
            N.append(np.sqrt(lorentzTDBC(wavelength)))
        elif layer['name'] == 'air':
            N.append(nOut)
    return N

def angles(refractiveIdxs,incidenceAngle):
    """
    refractiveIdxs: Refractive indices list
    incidenceAngle: Incidence angle in radians

    Outputs a list of angles for each refractive index, for the given incidence angle
    """

    angles = [incidenceAngle]
    for c,idx in enumerate(refractiveIdxs[:-1]):
        angles.append(np.arcsin((idx/refractiveIdxs[c+1])*np.sin(angles[-1])))

    return angles

def phases(cavity,refractiveIdxs,angles,wavelength):
    """
    cavity: Cavity list of dictionaries
    refractiveIdxs: Refractive index list
    angles: Angle list
    wavelength: Wavelength in nanometers

    Outputs phase list
    """

    thcks = [layer['thickness'] for layer in cavity]
    return [2*np.pi*idx*th*np.cos(angle)/wavelength for idx, th, angle in zip(refractiveIdxs,thcks,angles)]

def admittances(refractiveIdxs,angles,polarization='TM'):
    """
    refractiveIdxs: Refractive index list
    angles: Angle list
    polarization: Electric polarization, 'TM' or 'TE', 'TM' by default

    Outputs tilted optical admittance list, incidence admittance and polarization
    """
    n0 = 1 # Usually the same as output
    if polarization == 'TM':
        return [idx/np.cos(angle) for idx,angle in zip(refractiveIdxs,angles)], n0/np.cos(angles[0]), 'TM'
    elif polarization == 'TE':
        return [idx*np.cos(angle) for idx,angle in zip(refractiveIdxs,angles)], n0*np.cos(angles[0]), 'TE'

def TMatrixGenerator(deltas,etas):
    """
    deltas: Phase list
    etas: Tilted optical admittance list

    Outputs cavity matrix, matrices
    """

    Matrices = np.zeros((len(deltas),2,2),dtype=np.complex_)

    Matrices = [[[np.cos(dt),complex(0,-1)*np.sin(dt)/e],[complex(0,-1)*e*np.sin(dt),np.cos(dt)]] for dt, e in zip(deltas,etas)]

    M_ensemble = np.array([[1,0],[0,1]])

    for matrix in Matrices:
        M_ensemble = M_ensemble.dot(matrix)

    return M_ensemble,Matrices

def RTA(EnsembleMatrix,admittanceIn,admittanceOut):
    """
    EnsembleMatrix: Cavity Matrix
    admittanceIn: Tilted optical admittance In
    admittanceOut: Tilted optical admittance Out

    Outputs Reflection, Transmission, Absorption
    """

    BC = EnsembleMatrix.dot([[1],[admittanceOut]])
    Y = BC[1][0]/BC[0][0]
    
    R = ((admittanceIn-Y)/(admittanceIn+Y))*np.conj((admittanceIn-Y)/(admittanceIn+Y))
    T = (4*admittanceIn*np.real(admittanceOut))/((admittanceIn*BC[0][0]+BC[1][0])*np.conj(admittanceIn*BC[0][0]+BC[1][0]))
    A = (4*admittanceIn*np.real(BC[0][0]*np.conj(BC[1][0])-admittanceOut))/((admittanceIn*BC[0][0]+BC[1][0])*np.conj(admittanceIn*BC[0][0]+BC[1][0]))

    return R, T, A

def plotRTA(wavelength, RTAdata, icdAngle = 0, maxAngle = 90, toPlot='RTA', polarization = 'TM', figsize=(7,7), scaleFixed=True): 
    """
    wavelength: Wavelength list
    RTAdata: Data dictionary to plot 
    icdAngle: Incidence Angle, default 0°
    maxAngle: Maximum calculated angle, default 90°
    toPlot: RTA data to plot
    polarization: Light Polarization
    figsize: Canvas size
    scaleFixed: Fix scale
    incidenceAngles: Incidence Angle List

    Outputs a plot with the selected matrices on given incidence angle
    """
    
    a = icdAngle*len(RTAdata['R'][0,:])/maxAngle

    fig = plt.figure(figsize=figsize)
    for p in list(toPlot):
        if p == 'R':
            plt.plot(wavelength,RTAdata[p][floor(a),:],'r',label='Reflection')
        if p == 'T':
            plt.plot(wavelength,RTAdata[p][floor(a),:],'b',label='Transmission')
        if p == 'A':
            plt.plot(wavelength,RTAdata[p][floor(a),:],'g',label = 'Absorption')

    plt.legend()
    plt.ylabel('Intensity [a.u]')
    plt.xlabel('Wavelength [nm]')
    if scaleFixed:
        plt.xlim([380,830])
        plt.ylim([0,1])
    plt.title(f'Intensity spectrum at ~{icdAngle:.0f}°, {polarization} polarization', fontsize=20)
    plt.show()
