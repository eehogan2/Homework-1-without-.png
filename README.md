Homework-1-without-.png
=======================
{
 "metadata": {
  "name": "",
  "signature": "sha256:6898369147f4cf9dd5d4a4075ef3499fb8187a36ee6f561d9cbe669bd0dcf27b"
 },
 "nbformat": 3,
 "nbformat_minor": 0,
 "worksheets": [
  {
   "cells": [
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "%matplotlib inline \n",
      "#Plot all figures inside the notebook. Have to have this as can't plot as other popups.\n",
      "import netCDF4\n",
      "from mpl_toolkits.basemap import Basemap, shiftgrid\n",
      "from netCDF4 import Dataset, date2index\n",
      "import numpy as np\n",
      "import matplotlib.pyplot as plt\n",
      "from datetime import datetime"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 2
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "#Question 1\n",
      "\n",
      "#For January\n",
      "\n",
      "# specify an url, the JARKUS dataset in this case\n",
      "url = 'http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_pressure_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_pressure_omega.jnl'# create a dataset object\n",
      "dataset = netCDF4.Dataset(url)\n",
      "\n",
      "url2 = 'http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_pressure_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_pressure_uwnd.jnl'# create a dataset object\n",
      "dataset2 = netCDF4.Dataset(url2)\n",
      "\n",
      "url3 = 'http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_pressure_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_pressure_vwnd.jnl'# create a dataset object\n",
      "dataset3 = netCDF4.Dataset(url3)\n",
      "\n",
      "\n",
      "timevar = dataset.variables['TIME']\n",
      "levvar = dataset2.variables['LEV']\n",
      "\n",
      "print timevar[:]\n",
      "print levvar[:]\n",
      "print dataset2.variables['uwnd']\n",
      "\n",
      "timeindex = 0 # find time index for desired date. This is being used as the index for the pressure.\n",
      "presindex = 6 # this is the index for our month\n",
      "# read sst.  Will automatically create a masked array using\n",
      "# missing_value variable attribute. 'squeeze out' singleton dimensions.\n",
      "omega = dataset.variables['omega'][timeindex,presindex,:,:].squeeze()\n",
      "u = dataset2.variables['uwnd'][timeindex,11,:,:].squeeze()\n",
      "v = dataset3.variables['vwnd'][timeindex,11,:,:].squeeze()\n"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "[   0.   31.   59.   90.  120.  151.  181.  212.  243.  273.  304.  334.]\n",
        "[   10.    20.    30.    50.    70.   100.   150.   200.   250.   300.\n",
        "   400.   500.   600.   700.   850.   925.  1000.]\n",
        "<type 'netCDF4.Variable'>\n",
        "float32 uwnd(TIME, LEV, LAT, LON)\n",
        "    direction: IJKL\n",
        "    missing_value: -9.96921e+36\n",
        "    ferret_datatype: FLOAT\n",
        "    infile_datatype: FLOAT\n",
        "    _FillValue: -9.96921e+36\n",
        "    long_name: monthly long term mean u wind\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/pressure/uwnd\n",
        "    dataset_index: 1\n",
        "unlimited dimensions: \n",
        "current shape = (12, 17, 73, 144)\n",
        "filling off\n",
        "\n"
       ]
      }
     ],
     "prompt_number": 4
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset.variables['LAT'][:]\n",
      "lons = dataset.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,omega,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Omega and wind analysis for 500 hPa January')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f978bcd1b10>"
       ]
      }
     ],
     "prompt_number": 6
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "np.shape(omega)"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "pyout",
       "prompt_number": 5,
       "text": [
        "(73, 144)"
       ]
      }
     ],
     "prompt_number": 5
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "\n",
      "#Question 1\n",
      "\n",
      "#For July\n",
      "\n",
      "timeindex = 6 # find time index for desired date. This is being used as the index for the pressure.\n",
      "presindex = 6 # this is the index for our month\n",
      "\n",
      "omega = dataset.variables['omega'][timeindex,presindex,:,:].squeeze()\n",
      "u = dataset2.variables['uwnd'][timeindex,11,:,:].squeeze()\n",
      "v = dataset3.variables['vwnd'][timeindex,11,:,:].squeeze()\n"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 7
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset.variables['LAT'][:]\n",
      "lons = dataset.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,omega,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Omega and wind analysis for 500 hPa July')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png":  PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f978c564ad0>"
       ]
      }
     ],
     "prompt_number": 8
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Question 2\n",
      "\n",
      "# 850hPa relative humidity and winds for January\n",
      "\n",
      "url4='http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_pressure_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_pressure_rhum.jnl'\n",
      "dataset4 = netCDF4.Dataset(url4)\n",
      "print dataset4.variables['rhum']\n",
      "levs=dataset4.variables['LEV']\n",
      "print levs[:]\n",
      "time=dataset4.variables['TIME']\n",
      "print time[:]\n",
      "timeindex=0\n",
      "presindex=14\n",
      "\n",
      "rhum = dataset4.variables['rhum'][timeindex,5,:,:].squeeze()\n",
      "u = dataset2.variables['uwnd'][timeindex,presindex,:,:].squeeze()\n",
      "v = dataset3.variables['vwnd'][timeindex,presindex,:,:].squeeze()\n",
      "\n",
      "\n"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "<type 'netCDF4.Variable'>\n",
        "float32 rhum(TIME, LEV, LAT, LON)\n",
        "    direction: IJKL\n",
        "    missing_value: -9.96921e+36\n",
        "    ferret_datatype: FLOAT\n",
        "    infile_datatype: FLOAT\n",
        "    _FillValue: -9.96921e+36\n",
        "    long_name: monthly long term mean of relative humidity\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/pressure/rhum\n",
        "    dataset_index: 1\n",
        "unlimited dimensions: \n",
        "current shape = (12, 8, 73, 144)\n",
        "filling off\n",
        "\n",
        "[  300.   400.   500.   600.   700.   850.   925.  1000.]\n",
        "[   0.   31.   59.   90.  120.  151.  181.  212.  243.  273.  304.  334.]\n"
       ]
      }
     ],
     "prompt_number": 9
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset.variables['LAT'][:]\n",
      "lons = dataset.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,rhum,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Relative Humidity and wind analysis for 850 hPa January')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png":  PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f978bcee790>"
       ]
      }
     ],
     "prompt_number": 10
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Question 2\n",
      "\n",
      "# 850hPa relative humidity and winds for July\n",
      "\n",
      "timeindex=6\n",
      "presindex=14\n",
      "\n",
      "rhum = dataset4.variables['rhum'][timeindex,5,:,:].squeeze()\n",
      "u = dataset2.variables['uwnd'][timeindex,presindex,:,:].squeeze()\n",
      "v = dataset3.variables['vwnd'][timeindex,presindex,:,:].squeeze()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 11
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset.variables['LAT'][:]\n",
      "lons = dataset.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,rhum,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Relative Humidity and wind analysis for 850 hPa July')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f978be7f110>"
       ]
      }
     ],
     "prompt_number": 12
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Question 3\n",
      "\n",
      "# Net Surface sensible heat flux for January\n",
      "\n",
      "url5='http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_surface_gauss_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_surface_gauss_shtfl.jnl'\n",
      "dataset5 = netCDF4.Dataset(url5)\n",
      "\n",
      "time=dataset5.variables['TIME']\n",
      "print dataset5.variables['shtfl']\n",
      "print dataset5.variables['LAT']\n",
      "print time[:]\n",
      "\n",
      "timeindex=0\n",
      "\n",
      "shtfl = dataset5.variables['shtfl'][timeindex,:,:].squeeze()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "<type 'netCDF4.Variable'>\n",
        "float32 shtfl(TIME, LAT, LON)\n",
        "    direction: IJL\n",
        "    missing_value: -9.96921e+36\n",
        "    ferret_datatype: FLOAT\n",
        "    infile_datatype: FLOAT\n",
        "    _FillValue: -9.96921e+36\n",
        "    long_name: monthly long term mean of sensible heat net flux\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/shtfl\n",
        "    dataset_index: 1\n",
        "unlimited dimensions: \n",
        "current shape = (12, 94, 192)\n",
        "filling off\n",
        "\n",
        "<type 'netCDF4.Variable'>\n",
        "float64 LAT(LAT)\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/shtfl\n",
        "    grads_size: 94\n",
        "    grads_mapping: levels\n",
        "    minimum: -88.542\n",
        "    direction: J\n",
        "    grads_dim: y\n",
        "    resolution: 1.90413\n",
        "    point_spacing: uneven\n",
        "    units: degrees_north\n",
        "    bounds: lat_bnds\n",
        "    infile_datatype: DOUBLE\n",
        "    long_name: latitude\n",
        "    start: -88.542\n",
        "    maximum: 88.542\n",
        "    length: 94\n",
        "    orig_file_axname: lat\n",
        "    modulo: no\n",
        "    end: 88.542\n",
        "unlimited dimensions: \n",
        "current shape = (94,)\n",
        "filling off\n",
        "\n",
        "[   0.   31.   59.   90.  120.  151.  181.  212.  243.  273.  304.  334.]\n"
       ]
      }
     ],
     "prompt_number": 13
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset5.variables['LAT'][:]\n",
      "lons = dataset5.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,shtfl,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "#urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "#Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "#qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Net surface sensible heat flux for January')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f7aacf66210>"
       ]
      }
     ],
     "prompt_number": 13
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "#Question 3\n",
      "\n",
      "# Net Surface Sensible heat flux for July\n",
      "\n",
      "\n",
      "timeindex=6\n",
      "\n",
      "shtfl = dataset5.variables['shtfl'][timeindex,:,:].squeeze()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [],
     "prompt_number": 14
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset5.variables['LAT'][:]\n",
      "lons = dataset5.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,shtfl,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "#urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "#Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "#qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Net surface sensible heat flux for July')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f7aacf66f10>"
       ]
      }
     ],
     "prompt_number": 15
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [],
     "language": "python",
     "metadata": {},
     "outputs": []
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Question 4\n",
      "\n",
      "# Net Surface Latent heat flux for January.\n",
      "\n",
      "url7='http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_surface_gauss_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_surface_gauss_lhtfl.jnl'\n",
      "dataset7 = netCDF4.Dataset(url7)\n",
      "\n",
      "time=dataset7.variables['TIME']\n",
      "print dataset7.variables['lhtfl']\n",
      "print dataset7.variables['LAT']\n",
      "print time[:]\n",
      "\n",
      "timeindex=0\n",
      "\n",
      "lhtfl = dataset7.variables['lhtfl'][timeindex,:,:].squeeze()\n"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "<type 'netCDF4.Variable'>\n",
        "float32 lhtfl(TIME, LAT, LON)\n",
        "    direction: IJL\n",
        "    missing_value: -9.96921e+36\n",
        "    ferret_datatype: FLOAT\n",
        "    infile_datatype: FLOAT\n",
        "    _FillValue: -9.96921e+36\n",
        "    long_name: monthly long term mean of latent heat net flux\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/lhtfl\n",
        "    dataset_index: 1\n",
        "unlimited dimensions: \n",
        "current shape = (12, 94, 192)\n",
        "filling off\n",
        "\n",
        "<type 'netCDF4.Variable'>\n",
        "float64 LAT(LAT)\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/lhtfl\n",
        "    grads_size: 94\n",
        "    grads_mapping: levels\n",
        "    minimum: -88.542\n",
        "    direction: J\n",
        "    grads_dim: y\n",
        "    resolution: 1.90413\n",
        "    point_spacing: uneven\n",
        "    units: degrees_north\n",
        "    bounds: lat_bnds\n",
        "    infile_datatype: DOUBLE\n",
        "    long_name: latitude\n",
        "    start: -88.542\n",
        "    maximum: 88.542\n",
        "    length: 94\n",
        "    orig_file_axname: lat\n",
        "    modulo: no\n",
        "    end: 88.542\n",
        "unlimited dimensions: \n",
        "current shape = (94,)\n",
        "filling off\n",
        "\n",
        "[   0.   31.   59.   90.  120.  151.  181.  212.  243.  273.  304.  334.]\n"
       ]
      }
     ],
     "prompt_number": 22
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset7.variables['LAT'][:]\n",
      "lons = dataset7.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,lhtfl,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "#urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "#Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "#qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Net surface latent heat flux for January')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED
       "text": [
        "<matplotlib.figure.Figure at 0x7f7aad07aad0>"
       ]
      }
     ],
     "prompt_number": 23
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "# Question 4\n",
      "\n",
      "# Net Surface Latent heat flux for July.\n",
      "\n",
      "url7='http://apdrc.soest.hawaii.edu/thredds/dodsC/las/ncep_mon_surface_gauss_clima/data_apdrc.soest.hawaii.edu_dods_public_data_Reanalysis_Data_NCEP_NCEP_clima_surface_gauss_lhtfl.jnl'\n",
      "dataset7 = netCDF4.Dataset(url6)\n",
      "\n",
      "time=dataset7.variables['TIME']\n",
      "print dataset7.variables['lhtfl']\n",
      "print dataset7.variables['LAT']\n",
      "print time[:]\n",
      "\n",
      "timeindex=6\n",
      "\n",
      "lhtfl = dataset7.variables['lhtfl'][timeindex,:,:].squeeze()\n"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "output_type": "stream",
       "stream": "stdout",
       "text": [
        "<type 'netCDF4.Variable'>\n",
        "float32 lhtfl(TIME, LAT, LON)\n",
        "    direction: IJL\n",
        "    missing_value: -9.96921e+36\n",
        "    ferret_datatype: FLOAT\n",
        "    infile_datatype: FLOAT\n",
        "    _FillValue: -9.96921e+36\n",
        "    long_name: monthly long term mean of latent heat net flux\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/lhtfl\n",
        "    dataset_index: 1\n",
        "unlimited dimensions: \n",
        "current shape = (12, 94, 192)\n",
        "filling off\n",
        "\n",
        "<type 'netCDF4.Variable'>\n",
        "float64 LAT(LAT)\n",
        "    dataset: http://apdrc.soest.hawaii.edu/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/surface_gauss/lhtfl\n",
        "    grads_size: 94\n",
        "    grads_mapping: levels\n",
        "    minimum: -88.542\n",
        "    direction: J\n",
        "    grads_dim: y\n",
        "    resolution: 1.90413\n",
        "    point_spacing: uneven\n",
        "    units: degrees_north\n",
        "    bounds: lat_bnds\n",
        "    infile_datatype: DOUBLE\n",
        "    long_name: latitude\n",
        "    start: -88.542\n",
        "    maximum: 88.542\n",
        "    length: 94\n",
        "    orig_file_axname: lat\n",
        "    modulo: no\n",
        "    end: 88.542\n",
        "unlimited dimensions: \n",
        "current shape = (94,)\n",
        "filling off\n",
        "\n",
        "[   0.   31.   59.   90.  120.  151.  181.  212.  243.  273.  304.  334.]\n"
       ]
      }
     ],
     "prompt_number": 24
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [
      "lats = dataset7.variables['LAT'][:]\n",
      "lons = dataset7.variables['LON'][:]\n",
      "lons, lats = np.meshgrid(lons,lats)\n",
      "# create figure, axes instances.\n",
      "fig = plt.figure(figsize=(11,8.5))\n",
      "ax = fig.add_axes([0.05,0.05,0.9,0.9])\n",
      "# create Basemap instance.\n",
      "# coastlines not used, so resolution set to None to skip\n",
      "# continent processing (this speeds things up a bit)\n",
      "m =\\\n",
      "Basemap(llcrnrlon=0,llcrnrlat=-50,urcrnrlon=360,urcrnrlat=50,projection='mill')\n",
      "# draw line around map projection limb.\n",
      "# color background of map projection region.\n",
      "# missing values over land will show up this color.\n",
      "m.drawmapboundary(fill_color='0.3')\n",
      "# plot sst, then ice with pcolor\n",
      "im1 = m.pcolormesh(lons,lats,lhtfl,shading='flat',cmap=plt.cm.jet,latlon=True)\n",
      "# plot wind vectors on projection grid.\n",
      "# first, shift grid so it goes from -180 to 180 (instead of 0 to 360\n",
      "# in longitude).  Otherwise, interpolation is messed up.\n",
      "#ugrid,newlons = shiftgrid(180.,u,lons,start=False)\n",
      "#vgrid,newlons = shiftgrid(180.,v,lons,start=False)\n",
      "# transform vectors to projection grid.\n",
      "#urot,vrot,x,y = m.rotate_vector(u[::2,::2],v[::2,::2],lons[::2,::2],lats[::2,::2],returnxy=True)\n",
      "# now plot.\n",
      "#Q = m.quiver(x,y,urot,vrot) #or specify, e.g., width=0.003, scale=400)\n",
      "#qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')\n",
      "\n",
      "\n",
      "m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])\n",
      "m.drawmeridians(np.arange(-180.,180.,60.),labels=[1,0,0,1])\n",
      "# add colorbar\n",
      "cb = m.colorbar(im1,\"bottom\", size=\"5%\", pad=\"17%\")\n",
      "# add a title.\n",
      "ax.set_title('Net surface latent heat flux for July')\n",
      "m.drawcoastlines()\n",
      "plt.show()"
     ],
     "language": "python",
     "metadata": {},
     "outputs": [
      {
       "metadata": {},
       "output_type": "display_data",
       "png": PNG HAS BEEN REMOVED.
       "text": [
        "<matplotlib.figure.Figure at 0x7f7aac5b4e50>"
       ]
      }
     ],
     "prompt_number": 26
    },
    {
     "cell_type": "raw",
     "metadata": {},
     "source": [
      "Question 5 - Is there a connection between the 500hPa winds and 850hPa relative humidity? Why?\n",
      "\n",
      "Generally, where we have converging winds at 500hPa we would expect to see lifting (cooling air, leading to condensation and cloud) and therefore high relative humidity levels. The movement of winds are closely connected to the hadley cells in a climatology. Lower wind speeds occur where air is lifting/sinking along the edges of these cells. High relative humidity occurs near the equator, where air is being lifted, and low relative humidity occurs at higher latitudes, associated with the downward movement of air (e.g. at the Sahara). We can see lifting through the negative values of omega at 500hPa and these patterns corresponds to high relative humidity values.\n",
      "We also see low relative humidity values where we see positive omega values.\n",
      "\n",
      "There are slight converging winds at 500hPa over the Maritime continent and East-Central Pacific (in January) that may be responsible for these areas of high relative humidity. We can also see throughout the year high relative humidity values that are associated the upward movement of air due to elevated land such as the Andes/Rockies."
     ]
    },
    {
     "cell_type": "raw",
     "metadata": {},
     "source": [
      "Question 6 - How are the sensible and latent heat fluxes related to 500hPa omega and 850hPa RH? Why?\n",
      "\n",
      "We would expect a high sensible heat flux over land experiencing summer (Northern Hemisphere in July, Southern Hemisphere in January) as land has a low heat capacity and is quickly warmed. Heat therefore will be lost around these areas more readily than the ocean, or land that is experiencing winter. The high sensible heat flux over Australia in January is coupled with low relative humidity values. The same occurs over the Sahara in July. These areas are associated with little cloud (low relative humidity as the air is sinking - positive omega values), so the Sun's rays are not blocked by cloud, and the sensible heat flux is higher.\n",
      "\n",
      "\n",
      "Latent heat flux is associated with the loss of water from the surface into the atmosphere via evaporation.\n",
      "The evaporation of water (latent heat) that leads to cloud over the equator in January (seen at the equator through high relative humidity values) is occuring north and south of the equator, as the trade winds are evaporating water from the ocean and advecting this moist air to the equator. The converging winds rise (seen through the negative omega values) and cloud forms. The cloud that is produced at the equator via this process blocks the Sun from the surface and prevents evaporation at the equator, so the latent heat flux here is close to 0.\n",
      "\n",
      "In July, there is a large flux of latent heat from the West Indian Ocean which advects over India and the Himalayas, which can be seen in the high relative humidity and negative omega values around the region. This pattern is associated with the Indian monsoon."
     ]
    },
    {
     "cell_type": "code",
     "collapsed": false,
     "input": [],
     "language": "python",
     "metadata": {},
     "outputs": []
    }
   ],
   "metadata": {}
  }
 ]
}
