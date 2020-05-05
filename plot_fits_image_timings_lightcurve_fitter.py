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
from astropy import modeling
from scipy.optimize import curve_fit
import math 

# create a gaussian function to model and fit to data
def func(x, a, x0, sigma, yoffset):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+yoffset

#set this to the path that has the files below
path = "/Users/njscott/Documents/pulsating_wd/"

#put the file names here and a guess at the eclipse center
date20200216={"filename":"N20200216A0230b.fits","lc_filename":"ztf-16.lc","eclipse_center_guess":156,"pcolor":'b',"fitcolor":'b'} #these are all blue files
date20200217a={"filename":"N20200217A0319b.fits","lc_filename":"ztf-17a.lc","eclipse_center_guess":140,"pcolor":'b',"fitcolor":'b'}
date20200217b={"filename":"N20200217A0320b.fits","lc_filename":"ztf-17b.lc","eclipse_center_guess":80,"pcolor":'b',"fitcolor":'b'}
date20200218={"filename":"N20200218A0536b.fits","lc_filename":"ztf-18.lc","eclipse_center_guess":57,"pcolor":'b',"fitcolor":'b'}
date20200216r={"filename":"N20200216A0230r.fits","lc_filename":"ztf-16r.lc","eclipse_center_guess":152,"pcolor":'r',"fitcolor":'r'} #these are all red files
date20200217ar={"filename":"N20200217A0319r.fits","lc_filename":"ztf-17ar.lc","eclipse_center_guess":140,"pcolor":'r',"fitcolor":'r'}
date20200217br={"filename":"N20200217A0320r.fits","lc_filename":"ztf-17br.lc","eclipse_center_guess":80,"pcolor":'r',"fitcolor":'r'}
date20200218r={"filename":"N20200218A0536r.fits","lc_filename":"ztf-18r.lc","eclipse_center_guess":57,"pcolor":'r',"fitcolor":'r'}

#select which observation here:
#file=date20200216 #these are all blue files
#file=date20200217a
#file=date20200217b
#file=date20200218
file=date20200216r #these are all red files

#don't have these yet
#file=date20200217ar
#file=date20200217br
#file=date20200218r

fitsfile = path+file["filename"]
#file = askopenfilename()
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

#save the time table to a txt (unix time) and convert it to mjd and save a txt
np.savetxt(path+str(file["filename"])+'_timing.txt', end_of_frame_UTC, delimiter=" ") 
t=np.loadtxt(path+str(file["filename"])+'_timing.txt', dtype='str')
unixtime=Time(t, format='unix')
mjdtime=unixtime.to_value('mjd','long')
np.savetxt(path+str(file["filename"])+'_mjdtiming.txt', mjdtime, delimiter=" ") 

#plot the image
fig, ax = plt.subplots()
plt.figure(1)
plt.imshow(image_data[0,:,:])
plt.savefig(path+'image'+str(file["filename"])+'.png',dpi=200)

#plot the mjd
plt.figure(2)
plt.plot(mjdtime[0:frames-1], color=str(file["pcolor"]), marker='o', linestyle='solid', linewidth=1, markersize=2)
plt.xlabel('frame')
plt.ylabel('MJD (s)')
plt.savefig(path+'timing'+str(file["filename"])+'.png',dpi=200)

#work out the bjd
#https://mail.python.org/pipermail/astropy/2014-April/002844.html
#https://docs.astropy.org/en/stable/time/ 
#http://astroutils.astronomy.ohio-state.edu/time/

#time to barycenter correction -needs target ra, dec, and observatory location. it pulls these from the fits header
#enter target in qoutes below
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

np.savetxt(path+str(file["filename"])+'_ltt_barytiming.txt', ltt_bary, delimiter=" ", fmt='%s') 
bjdtime=mjdtime+ltt_bary
np.savetxt(path+str(file["filename"])+'bjdtiming.txt', bjdtime, delimiter=" ", fmt='%s') 

#read the data in using pandas
dat_file = path+file["lc_filename"]
df = pd.read_csv(dat_file, index_col=None, delim_whitespace=True, names=['frame','mag'])

#select region of interest (optional)
framestart=0 #starting from 0
frameend=frames
numframes=frameend-framestart
frames=df['frame'].values[framestart:frameend]
mag=df['mag'].values[framestart:frameend]

x = np.linspace(framestart, frameend, 1000) #array for model fit

#https://python4esac.github.io/fitting/examples1d.html
popt, pcov = curve_fit(func, frames, mag, p0=[-1,file["eclipse_center_guess"],1,9])
#popt returns the best fit values for parameters of the given model (func)
print(popt)
ym = func(frames, popt[0], popt[1], popt[2], popt[3])

#plots
fig, ax = plt.subplots()
#plt.plot(frames,mag, color='blue', marker='o', linestyle='solid', linewidth=1, markersize=2)
plt.scatter(frames,mag, c=file["pcolor"], s=10)
ax.set_xlim(framestart, frameend)
#set y range 
maglow=8.5
magup=13
ax.set_ylim(magup, maglow)
ax.set_xlabel('frame')
ax.set_ylabel('instrument magnitude')
ax.set_title('ecliping WD binary')

plt.plot(frames, ym, c=file["fitcolor"], label='Best fit')
plt.title(str(file["lc_filename"]))
ax.legend()

eclipse_center = popt[1]
print('eclipse center frame= '+str(eclipse_center))
frame=int(round(popt[1]))

#interpolate timing points for the eclipse center
eclipse_frame1 = math.floor(eclipse_center)
eclipse_frame2 = math.ceil(eclipse_center)
def findYPoint(xa,xb,ya,yb,xc):
    m = (ya - yb) / (xa - xb)
    yc = (xc - xb) * m + yb
    return yc
unix_fit = findYPoint(eclipse_frame1,eclipse_frame2,unixtime[eclipse_frame1],unixtime[eclipse_frame2],eclipse_center)
mjd_fit = findYPoint(eclipse_frame1,eclipse_frame2,mjdtime[eclipse_frame1],mjdtime[eclipse_frame2],eclipse_center)
ltt_bary_fit = findYPoint(eclipse_frame1,eclipse_frame2,ltt_bary[eclipse_frame1],ltt_bary[eclipse_frame2],eclipse_center)
bjd_fit = findYPoint(eclipse_frame1,eclipse_frame2,bjdtime[eclipse_frame1],bjdtime[eclipse_frame2],eclipse_center)

print(ax)
label='fit (frame) = '+str(eclipse_center) #str(frame)
labeltime='BJD = '+bjd_fit.to_value('jd', 'str')
plt.annotate(label, (frame,magup-0.25), textcoords="offset points", xytext=(0,0), ha='center')
plt.annotate(labeltime, (frame,magup-0.09), textcoords="offset points", xytext=(0,0), ha='center')

plt.savefig(path+'lc'+str(file["lc_filename"])+'.png',dpi=200)

print('eclipse frame (rounded)= '+str(frame))
#print ('eclipse fit (unix) = '+str(unixtime[frame]))
#print ('eclipse fit (mjd) = '+str(mjdtime[frame]))
#print ('eclipse fit (ltt_bary) = '+str(ltt_bary[frame]))
#print ('eclipse fit (bjd) = '+str(bjdtime[frame]))

print('eclipse fit time (linear interpolated)')
print ('interp eclipse fit (unix) = '+unix_fit.to_value('unix', 'str'))
print ('interp eclipse fit (mjd) = '+str(mjd_fit))
print ('interp eclipse fit (ltt) = '+ltt_bary_fit.to_value('jd', 'str'))
print ('interp eclipse fit (bjd) = '+bjd_fit.to_value('jd', 'str'))

fits=np.array([unix_fit,mjd_fit,ltt_bary_fit,bjd_fit])
headers = np.array(["unix","mjd","ltt","bjd"])
data = np.array([headers,fits])
data = data.T  #here you transpose your data, so to have it in two columns

np.savetxt(path+str(file["filename"])+"interpeclipsefits.txt", data, fmt=['%s','%s'])

plt.show()