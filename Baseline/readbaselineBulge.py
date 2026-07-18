import pandas as pd 
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams 
import scipy.special as ss
import warnings
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

from numpy import matrix
from matplotlib import colors

cm = colors.ListedColormap(['purple', 'blue', 'darkgreen','yellowgreen', 'orange', 'red'])
###############################################################################

nam0 = ['observationId','fieldRA','fieldDec','observationStartMJD','flush_by_mjd','visitExposureTime','band','filter','rotSkyPos',
'rotSkyPos_desired','numExposures','airmass','seeingFwhm500','seeingFwhmEff','seeingFwhmGeom',
'skyBrightness','night','slewTime','visitTime','slewDistance','fiveSigmaDepth','altitude','azimuth','paraAngle',
'pseudoParaAngle','cloud','moonAlt','sunAlt','scheduler_note','target_name','target_id','observationStartLST',
'rotTelPos','rotTelPos_backup','moonAz','sunAz','sunRA','sunDec','moonRA','moonDec','moonDistance','solarElong',
'moonPhase','cummTelAz','observation_reason','science_program','cloud_extinction', 'test1','test2']##49

nam1 = ['ID','RA', 'DEC', 't', 'texp', 'filter', 'air', 'see', 'skyB', 'Tv', 'sig5', 'target']#12

idx = [0, 1, 2, 3, 5, 7, 11, 13, 15, 18, 20, 29]

assert len(idx) == len(nam1), "idx and new_names must match"
assert max(idx) < len(nam0), "idx out of range"


num, lon, lat = np.loadtxt('./layout_7f_3.centers', unpack=True) ## Bulge
# Numbers below may change based on the Roman Bulge survey details
l0 = -0.219 - 0.2 - 3.5 / 2
l1 = 1.4134 + 0.2 + 3.5 / 2
l2 = 1.0053 - 0.2 - 3.5 / 2

b0 = -1.64  - 0.2 - 3.5 / 2
b1 = -0.85  + 0.2 + 3.5 / 2
b2 = -1.64  + 0.2 + 3.5 / 2

#RA0, RA1, DEC0, DEC1 = float(75.0 - 3.5 / 2.0), float(90.0 + 3.5 / 2.0),  float(-75.0 -3.5 / 2.0), float(-60.0 + 3.5 / 2.0)##Bulge
#fil = open("./Bulgebaseline.dat", "w")
#fil.close()

#################################################################################

selected_cols = [nam0[i] for i in idx]
query = f"""
SELECT {', '.join(selected_cols)}
FROM observations
"""

conn = sqlite3.connect("baseline_v5.1.0_10yrs.db")
df = pd.read_sql_query(query, conn)
conn.close()

#print(df.keys())
df.columns = nam1
#print(df.keys())

nm  = int(len(df))
ra  = np.zeros((nm))
dec = np.zeros((nm))
#dis=np.sqrt( (ra-RA0)**2.0 + (dec-DEC0)**2.0)## degree

# Create SkyCoord object
coords = SkyCoord(
    ra  = df["RA"].values * u.deg,
    dec = df["DEC"].values * u.deg,
    frame = "icrs"
)

# Add Galactic coordinates to DataFrame
l = coords.galactic.l.deg
df["l"] = Angle(l * u.deg).wrap_at(180 * u.deg).degree
df["b"] = coords.galactic.b.deg

l = df["l"].values
b = df["b"].values

nvis = int(20000)
tstA = np.zeros((nvis, 13))
nfil = 0
nr   = 0

for i in range(nm):   
    if (l[i] >= l0 and l[i] <= l1 and b[i] >= b0 and b[i] <= b1) and not (l[i] < l2 and b[i] > b2): 
        #print(df["RA"][i], df["DEC"][i], df['l'][i], df['b'][i], df['t'][i], df['filter'][i]) 
        #print(df['air'][i], df['see'][i], df['skyB'][i] )
        #print(df['Tv'][i],  df['sig5'][i], df['target'][i],   df['ID'][i] , df['texp'][i] )
        #print("**********************************************************")
        #input("Enter a number")
        if(df['filter'][i]=='u'):  nfil=0
        if(df['filter'][i]=='g'):  nfil=1
        if(df['filter'][i]=='r'):  nfil=2
        if(df['filter'][i]=='i'):  nfil=3
        if(df['filter'][i]=='z'):  nfil=4
        if(df['filter'][i]=='y'):  nfil=5
        
        tstA[nr,:] = np.array([df['ID'][i], df['RA'][i], df['DEC'][i], df['l'][i], df['b'][i], df['t'][i], nfil, df['air'][i], df['see'][i], df['skyB'][i], df['Tv'][i], df['sig5'][i], df['texp'][i]])##13
        
        nr+=1

print(nr)
print("No. year of observations:" , float(np.max(tstA[:nr,5]) - np.min(tstA[:nr,5])) / 365.2425)

###############################################################################
dist  = np.zeros((nr))
fil   = open("./BulgeBaseline.dat", "a")
fil.write("#ID  RA  Dec  l  b  time  filter  airmass  seeing  skyBrightness visittime sigma5 targetname distance\n")

tst   = np.zeros((nr, 14))
idx   = np.argsort(tstA[:nr, 5])

time0 = float(tstA[int(idx[0]), 5])

for i in range(nr):
    tst[i,:-1] = tstA[int(idx[i]), :]
    tst[i,5] = tst[i, 5] - time0

    if(i > 0): 
        dist[i] = np.sqrt((tst[i, 1] - tst[i-1, 1]) ** 2.0 + (tst[i, 2] - tst[i-1, 2]) ** 2.0)

    tst[i,-1] = dist[i]

    np.savetxt(fil, tst[i,:].reshape((-1, 14)),
               fmt ="%d  %.6f  %.6f  %.6f  %.6f  %.6f  %d  %.6f  %.6f  %.6f  %.1f  %.6f  %.1f  %.6f") ##14 
            
fil.close()
print("Distance:  ", np.mean(dist[1:nr]), np.min(dist[1:nr]), np.max(dist[1:nr]))

###############################################################################

namm = [r"$ID$", r"$\rm{RA}$", r"$\rm{DEC}$", "l", "b", r"$\rm{time}$", r"$\rm{Filter}$", r"$\rm{airmass}$", r"$\rm{seeing}$", r"$\rm{Sky_{Brightness}}$", r"$\rm{Time}-\rm{visit}$", r"$\rm{Sigma5}$", r"$t_{EXP}$"]
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
ax1=fig.add_subplot(111)

for i in range(nr): 
    nc = abs(int(tst[i,6]))
    plt.plot(tst[i,3], tst[i,4], ".", markersize=7.2, color=cm(nc))##"o", markersize=3, color=col[tst[i,4]])
#plt.scatter(tst[:,3], tst[:,4], color=cm(abs(tst[:,6])))
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
ax1.set_aspect('equal', adjustable='box')
plt.xlim(-3 - 3.5 / 2 , 3 + 3.5 / 2)
plt.ylim(-3 - 3.5 / 2 , 3 + 3.5 / 2)
plt.xlabel(r"$l~[\rm{deg}]$", fontsize=20, labelpad=0.05)
plt.ylabel(r"$b~[\rm{deg}]$", fontsize=20, labelpad=0.05)
ax1.invert_xaxis()
fig=plt.gcf()
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./lbBulge.jpg" , dpi=200)

###############################################################################
for i in range(13):
    if i == 3 or i == 4:
        continue
    plt.clf()
    plt.cla()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(tst[:nr,i],30,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.5,rwidth=1.5)
    y_vals = ax.get_yticks()
    ax.set_yticks(y_vals)
    ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/nr))) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    ax.set_xlabel(str(namm[i]),fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.legend(prop={"size":12.5})
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./jpg/histBulge{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted *****************************" )   


################################################################################


## Cadence Histogram
plt.clf()
plt.cla()
fig= plt.figure(figsize=(8,6))
ax= plt.gca()              
plt.hist(dist,100,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.5,rwidth=1.5)
#y_vals = ax.get_yticks()
#ax.set_yticks(y_vals)
#ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/nr))) for x in y_vals]) 
#y_vals = ax.get_yticks()
#plt.ylim([np.min(y_vals), np.max(y_vals)])
ax.set_ylabel(r"$\rm{Distribution}$",fontsize=19,labelpad=0.1)
ax.set_xlabel(r"$\rm{distance}(\rm{degree})$",fontsize=19,labelpad=0.1)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlim(0.0,50.0)
plt.legend(prop={"size":12.5})
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.savefig("./jpg/DistanceBulge.jpg" , dpi=200)

################################################################################




















