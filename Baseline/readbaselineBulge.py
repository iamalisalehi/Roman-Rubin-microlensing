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

nam1 = ['ID','RA', 'DEC', 't', 'texp', 'filter', 'air', 'see', 'skyB', 'Tv', 'sig5', 'target', 'rot']#13

idx = [0, 1, 2, 3, 5, 7, 11, 13, 15, 18, 20, 29, 8]   # 8 = rotSkyPos (Deviation 80)

assert len(idx) == len(nam1), "idx and new_names must match"
assert max(idx) < len(nam0), "idx out of range"


# WHICH POINTINGS (Deviation 69). The simulator scans every sky point within
# gbtds_geometry.scan_reach() (= 2 x 1.75 deg + the field reach) of a Roman field centre,
# spring or autumn roll. A Rubin visit images a sightline if its pointing centre is within
# FoV = 1.75 deg of it, so EVERY pointing centred within scan_reach() + FoV of a field centre
# can image some scanned sightline and must be in the list -- otherwise sightlines near the
# region's edge silently undercount their Rubin visits. Until Deviation 69 the list was cut by
# a box (l0..l1, b0..b1, minus an l2/b2 corner) on the pointing CENTRES around an older layout:
# 3,686 visits, against 12,348 by this rule (every one of the old 3,686 is among them).
# src/run/sightlines.cpp refuses a list with a pointing outside this reach.
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
import gbtds_geometry as G          # noqa: E402
from cparams import P               # noqa: E402
# Deviation 80: Rubin's coverage is now LSSTCam's active silicon (Baseline/lsstcam_fov), which
# reaches 1.94 deg from the boresight (the radius rubin_scheduler crops it at), not 1.75.
RUBIN_MAX_RADIUS = 1.94
REACH = G.scan_reach() + RUBIN_MAX_RADIUS

# DAY 0 OF THE SIMULATION CLOCK, PINNED. It used to be "the earliest selected visit", which
# made the clock depend on the selection: the rule above reaches pointings observed from MJD
# 60981.0 (2025-11-01), which would have moved day 0 back by 160 d and silently shifted
# Roman's placement on the clock with it. Day 0 stays MJD 61141.312002288 = 2026-04-11, the
# first visit of the original selection (observationId 74334), on which every run, document
# and the generator's MISSION_START_DAY are defined. Visits before it, or after Tobs, fall
# outside the simulated window and are dropped (40 before day 0, none after, on v5.1.0).
TIME0_MJD = 61141.312002288
TOBS_DAYS = P.Tobs                   # config/parameters.h Tobs

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

nvis = int(len(df))
tstA = np.zeros((nvis, 14))
nfil = 0
nr   = 0

reach_dist = np.min([np.hypot(l - fl, b - fb) for fl, fb, _ in G.placements()], axis=0)
t_sim = df["t"].values - TIME0_MJD
in_reach = reach_dist <= REACH
in_clock = (t_sim >= 0.0) & (t_sim <= TOBS_DAYS)
print(f"pointings within {REACH:.4f} deg of a Roman field: {in_reach.sum()}; "
      f"dropped outside the clock [0, {TOBS_DAYS:.1f}] d: {(in_reach & ~in_clock).sum()}")
for i in range(nm):   
    if in_reach[i] and in_clock[i]: 
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
        
        tstA[nr,:] = np.array([df['ID'][i], df['RA'][i], df['DEC'][i], df['l'][i], df['b'][i], df['t'][i], nfil, df['air'][i], df['see'][i], df['skyB'][i], df['Tv'][i], df['sig5'][i], df['texp'][i], df['rot'][i]])##14
        
        nr+=1

print(nr)
print("No. year of observations:" , float(np.max(tstA[:nr,5]) - np.min(tstA[:nr,5])) / 365.2425)

###############################################################################
dist  = np.zeros((nr))
# "w", NOT "a": append mode silently doubles the file if the script is run twice,
# which embeds a second header mid-file. The C++ reader's operator>> then fails on
# that header and zero-fills every remaining record without any CHECK firing.
fil   = open("./BulgeBaseline.dat", "w")
fil.write("#ID  RA  Dec  l  b  time  filter  airmass  seeing  skyBrightness visittime sigma5 texp distance rotSkyPos\n")
print(f"BulgeBaseline.dat: writing {nr} visit rows -- afterwards run `python3 tools/sync_data_products.py` from the repo root to update Nl in config/data_products.h")

tst   = np.zeros((nr, 15))
idx   = np.argsort(tstA[:nr, 5])

time0 = TIME0_MJD

for i in range(nr):
    row = tstA[int(idx[i]), :]
    tst[i, :13] = row[:13]
    tst[i, 5] = tst[i, 5] - time0

    if(i > 0): 
        dist[i] = np.sqrt((tst[i, 1] - tst[i-1, 1]) ** 2.0 + (tst[i, 2] - tst[i-1, 2]) ** 2.0)

    tst[i, 13] = dist[i]
    tst[i, 14] = row[13]          # rotSkyPos [deg] (Deviation 80)

    np.savetxt(fil, tst[i,:].reshape((-1, 15)),
               fmt ="%d  %.6f  %.6f  %.6f  %.6f  %.6f  %d  %.6f  %.6f  %.6f  %.1f  %.6f  %.1f  %.6f  %.6f") ##15 
            
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
fig.savefig("../pics/lbBulge.jpg" , dpi=200)

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
    fig.savefig("../pics/histBulge{0:d}.jpg".format(i),dpi=200)
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
fig.savefig("../pics/DistanceBulge.jpg" , dpi=200)

################################################################################




















