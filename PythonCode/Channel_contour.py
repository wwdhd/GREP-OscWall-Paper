"""
       FILE: channel.py
DESCRIPTION: Post processes the channel flow case, displaying contours and plots
"""
import csv
import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.misc import avg_over_axis
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker

mpl.use('TkAgg')

plt.ioff()

REFPATH="./"
INPUT_FILE="input.json"

#### TIME INPUT ######
t = "140000" 




def main():
    # Gathering the data
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]
    postprocess.load(time=[t])
    
    umean = postprocess.fields["umean"].data[t]
    vmean = postprocess.fields["vmean"].data[t]
    wmean = postprocess.fields["wmean"].data[t]  
    
    pmean = postprocess.fields["pmean"].data[t]
    
    uumean = postprocess.fields["uumean"].data[t]
    vvmean = postprocess.fields["vvmean"].data[t]
    wwmean = postprocess.fields["wwmean"].data[t]  
    
    uvmean = postprocess.fields["uvmean"].data[t]
    uwmean = postprocess.fields["uwmean"].data[t]
    vwmean = postprocess.fields["vwmean"].data[t]  


    print("========================================================================")
    #Calculating flow properties
    print("FLOW PROPERTIES")
    retau = 180
    reb = (np.exp(np.log(retau/0.09) / 0.88) / 2) 
    print(f"Friction Reynolds Number (Re_tau) {retau:.4f}")
    print(f"Bulk Reynolds Number (Re_b) {reb:.4f}")
    
    rho = 1  # density
    nu = 1/reb  # kinematic viscosity
    
    print(f"Kinematic Viscosity (nu)", nu)

    #Printing Mesh Size for Double Check
    print("Mesh Size")
    print(f"z elements = {len(umean[0][0])}") # Z
    print(f"y elements = {len(umean[0])}") # Y
    print(f"x elements = {len(umean)}") # X
    
    print("========================================================================")
    
    #Slice Plane
    x_coordinates = []
    for f in range(0,len(umean)):
        x_float = f*mesh.dx ; x_coordinates.append(x_float)
    
    y_coordinates = []
    for f in range(0,len(umean[0])):
        y_float = f*mesh.dy ; y_coordinates.append(y_float)
    
    z_coordinates = []
    for f in range(0,len(umean[0][0])):
        z_float = f*mesh.dz ; z_coordinates.append(z_float)
        
    X, Y = np.meshgrid(x_coordinates, y_coordinates, indexing='ij')  
    Y2, Z = np.meshgrid(y_coordinates, z_coordinates, indexing='ij')  
    dim = mesh.Lz/2
    dim1 = dim/mesh.dz

    ub = np.mean(umean, axis=(0, 2))  # Averaging over x and z
    ub = np.trapezoid(ub, y_coordinates) / (y_coordinates[-1] - y_coordinates[0])

    print(f"Bulk velocity (U_b): {ub}")

    tau_w = tau_wall_calculation(pmean, mesh, rho, ub)

    #print(tau_w)

    #print(mesh.Lx)
    
    
def tau_wall_calculation(pmean, mesh, rho, ub):
    # Compute mean pressure at inlet and outlet
    inlet_avg = np.mean(pmean[0, :, :])  # Mean pressure at inlet
    outlet_avg = np.mean(pmean[-1, :, :])  # Mean pressure at outlet

    pressure_scale = 0.5 * rho * ub**2
    scaled_dp_dx = (outlet_avg - inlet_avg) / mesh.Lx * pressure_scale
    tau_w = -scaled_dp_dx


    print(pressure_scale)
    print(f"Scaleddpdx: {scaled_dp_dx}")
    print(f"Inlet pressure avg: {inlet_avg}")
    print(f"Outlet pressure avg: {outlet_avg}")

    return tau_w

    
if __name__ == "__main__":
    main()
