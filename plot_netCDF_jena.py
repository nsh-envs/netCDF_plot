# Author: Niels S Hvidberg
# Description: Python program, converted from
#              a jupyter notebook.

# Import packages
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import datetime as dt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Printing text in colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def show(xr_dataset):
    """
    This function gives you a way to display xarray datasets per command, instead of implicitly
    by calling the dataset instance as the last thing in a cell.
    """
    from IPython.display import display, HTML
    display(HTML(xr.core.formatting_html.dataset_repr(xr_dataset)))

def xr_open_old(data_dir, days=[1, 31], hours=[0, 23]):
    infile = []
    dataset = []
    for i in range(days[0], days[1]+1):
        if (i == 1):
            z = 1
        else:
            z = 0
    
        for j in range(z+hours[0], hours[1]+1):
            name = f'201101{i:02d}{j:02d}.nc'
            infile.append(data_dir+'/'+name)
            dataset.append(xr.open_dataset(infile[-1]))
    
    return dataset

def xr_open(data_dir, days=[1, 31], hours=[0, 23]):
    infile = []
    dataset = []
    for i in range(days[0], days[1]+1):
        if (i == 1):
            z = 1
        else:
            z = 0
    
        for j in range(z+hours[0], hours[1]+1):
            name = f'1_201101{i:02d}{j:02d}.nc'
            infile.append(data_dir+'/'+name)
            dataset.append(xr.open_dataset(infile[-1]))
    
    return dataset

def JENA_plotter(dataset, species, isel={'mtime':0}, sel={}, close=False, central_longitude=-32, output_dir='plots', figspec='pCO2'):
    """
    Docstring for "JENA_plotter":
        This function is used to plot netCDF files on a polar
        stereographic map.

        dataset  : A list of xarray-dataset objects with 3D field data.
        species  : The kind to plot. For DEHM CO2 the tracers are called CO2.ntrac 
        isel/sel : A dictionary containing the variable and value or slices to choose
                   in the form {var:val/slice/index}. "isel" goes
                   by index and "sel" goes by value.
        meandim  : Calculate the mean over this dimension. Most often 'time' or 'z'.

        The function creates a plot and optionally saves and closes the plot again. The
        plot is saved under the name:
        
            f'{output_dir}/{figname}.png'
    """
    # Select values by index/slice/value.
    if isel:
        da = dataset[species].isel(isel)
    if sel:
        da = dataset[species].sel(sel)
    else:
        da = dataset[species].isel(mtime=0)

    # Set latitude and longitude values
    lat, lon = da['lat'].values, da['lon'].values
    
    # plot the first layer
    pdata = da.values
    
    # This is the map projection we want to plot *onto*
    map_proj = ccrs.NorthPolarStereo(central_longitude=central_longitude)
    fig, ax = plt.subplots(subplot_kw={'projection': map_proj}, figsize=(7,7))

    # Set appearance for coastlines
    ax.coastlines(color='black',zorder=2,alpha=1,linewidth=1.5)

    # Set appearance for gridlines
    gl = ax.gridlines(draw_labels=True, zorder=1, alpha=0.7, linewidth=1, crs=ccrs.PlateCarree())

    # Set label size
    gl.xlabel_style = {'size':7}
    gl.ylabel_style = {'size':7}

    # Set colormap
    cmap = plt.get_cmap('jet')
    cmap.set_bad('w',1.0)
    
    # Plot image
    im = ax.pcolormesh(lon,
                       lat,
                       pdata,
                       cmap=cmap,
#                      vmin=-0.1, # Min value for colormap (especially used for videos)
#                      vmax=14,   # Max value for colormap
                       transform=ccrs.PlateCarree(),
                       zorder=0,  # ?
                       shading='auto')

    # Set Colorbar
    fig.colorbar(im, ax=ax)

    # Check if plot should be closed and saved instead of shown (Used for creating frames)
    if close:
        figname = f'plot_{figspec}.png'
        plt.savefig(f'{output_dir}/{figname}.png')
        plt.close(fig)

## END OF DEFINITIONS #########################################################################
# Choose case and scroll down to the relevant section
case = 11

# DEHM data
if case == 1 :
    dir1 = 'netcdf1'
    dir2 = 'netcdf2'
    dir3 = 'netcdf3'
    dir4 = 'netcdf4'
    
    ds1 = xr_open_old(dir1)
    ds2 = xr_open_old(dir2)
    ds3 = xr_open_old(dir3)
    ds4 = xr_open_old(dir4)
    ds = ds1
    show(ds[0])

if case == 11 :
    dir1 = 'only_cams_glob/2011/netcdf'
    dir2 = 'only_jena_ocn/2011/netcdf'
    
    ds1 = xr_open(dir1)
    ds2 = xr_open(dir1)
    ds = ds1
    show(ds[0])

# CAMS data
if case == 2 :
    infile_cams       = '/home/niels/Documents/Data/CAMS/CAMS-GLOB/CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_co2_excl_short-cycle_org_C_v5.3_monthly.nc'
    ds_cams  = xr.open_dataset(infile_cams)
    ds       = ds_cams
    show(ds)

# Jena data
if case == 3 : 
    infile_jena = '/home/niels/Documents/Data/JENA/oc_v2023.pCO2.nc'
    ds_jena = xr.open_dataset(infile_jena)
    ds = ds_jena
    show(ds)

## END OF DATA SELECTION, BEGINNING OF SPECIFIC CODE ##########################################

# Find the index for a specific time:
# 2011 = 19723
sel_time = '2011-01'
_ = np.argwhere(ds.mtime.values==ds.mtime.sel(mtime=sel_time,method='nearest').values)[0]
del _
# Example of using slice to select slices of the data. 
# This can both be done with isel (using index) or sel (using values)
_ = ds['mtime'].isel({'mtime':slice(19723,19754)})
del _

species = 'pCO2'
sel={'mtime':ds['mtime'].values[0], 'lat':slice(55.0,57.0), 'lon':slice(8.0,16.0)}
central_longitude=-32
# print(ds[species].sel(sel)) # Testing settings

JENA_plotter(dataset=ds,
             species=species,
             sel=sel,
             central_longitude=central_longitude,
             output_dir='plots',
             figspec='pCO2',
            )