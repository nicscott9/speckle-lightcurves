import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import numpy as np
from astropy.time import Time
import sys
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import math 
import tkinter.filedialog

fitsfile = tkinter.filedialog.askopenfilename() #click the fits file you want to generate timing info for
path = tkinter.filedialog.askdirectory()+'/' #doubleclick the folder you want to save the files in
hdu_list=fits.open(fitsfile)
fits.info(fitsfile)
image_data = fits.getdata(fitsfile,ext=0)
#end_of_frame_UTC = Table(hdu_list[1].data)
end_of_frame_UTC = Table.read(fitsfile, hdu=1)

#fits keywords to print
keywords=['ACT','EXPTIME','OBSTIME','EXPENDTM','KCT','EXPOSURE','MJD-OBS','JD','UTC','SPKLEDAT','OBJECT','RA','DEC']
for i in keywords:
    print(i+' = '+str(hdu_list[0].header[i]))  #print FITS header keyword value

frames=hdu_list[0].header['NAXIS3'] #number of frames (counts from 0)
name = str(hdu_list[0].header['UTC'])+str(hdu_list[0].header['OBJECT'])

#save the time table to a txt (unix time) and convert it to mjd and save a txt
np.savetxt(path+name+'_timing.txt', end_of_frame_UTC, delimiter=" ") 
t=np.loadtxt(path+name+'_timing.txt', dtype='str')
unixtime=Time(t, format='unix')
mjdtime=unixtime.to_value('mjd','long')
np.savetxt(path+name+'_mjdtiming.txt', mjdtime, delimiter=" ") 

#plot the image
fig, ax = plt.subplots()
plt.figure(1)
plt.imshow(image_data[0,:,:])
plt.savefig(path+'image'+name+'.png',dpi=200)

#plot the mjd
plt.figure(2)
plt.plot(mjdtime[0:frames-1], color='k', marker='o', linestyle='solid', linewidth=1, markersize=2)
plt.xlabel('frame')
plt.ylabel('MJD (s)')
plt.savefig(path+'timing'+name+'.png',dpi=200)

#work out the bjd
#https://mail.python.org/pipermail/astropy/2014-April/002844.html
#https://docs.astropy.org/en/stable/time/ 
#http://astroutils.astronomy.ohio-state.edu/time/

#time to barycenter correction -needs target ra, dec, and observatory location. it pulls these from the fits header
#enter target in qoutes below if you don't want to pull the object from the header
#target="ip_peg"
#obj=coord.SkyCoord.from_name(target)  
ra = hdu_list[0].header['RA']
dec = hdu_list[0].header['DEC']

obj = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='fk5')
#l=coord.EarthLocation.get_site_names()
telescopedict = {
  "Gemini-North": "gemini_north",
  "Gemini-South": "gemini_south",
  "WIYN": "kitt peak"
}
telescope = telescopedict[hdu_list[0].header['OBSERVAT']]   
times = Time(mjdtime, format='mjd', scale='utc', location=coord.EarthLocation.of_site(telescope))  
ltt_bary = times.light_travel_time(obj)  

np.savetxt(path+name+'_ltt_barytiming.txt', ltt_bary, delimiter=" ", fmt='%s') 
bjdtime=mjdtime+ltt_bary
np.savetxt(path+name+'bjdtiming.txt', bjdtime, delimiter=" ", fmt='%s') 

plt.show()