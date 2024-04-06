"""
The script has been developed by Trent Turbyfill with the oversight of Dr. Dennis Wagenaar. 
This script is the master program complementary to holland_model.py. 

The script is an derived from: 
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. 
Generation of a global synthetic tropical cyclone hazard dataset using STORM. 
Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

The methodology is heavily inspired by 
Lin, N., and Chavas, D. ( 2012), On hurricane parametric wind and applications in storm surge modeling, 
J. Geophys. Res., 117, D09120, doi:10.1029/2011JD017126.

The original script which this is based on can be accessed here:
https://github.com/NBloemendaal/STORM-return-periods/blob/master/masterprogram.py
"""

# Import required libraries
import pandas as pd
import numpy as np
import math
from scipy import spatial
import holland_model as hm

################################################
#               PRE-PROCESSING                 #
################################################

# Function definitions
def haversine(lon1, lat1, lon2, lat2):
    # The radius of the Earth in kilometers and conversion from degrees to radians
    R, rad = 6371.0, math.radians
    dlon, dlat = rad(lon2 - lon1), rad(lat2 - lat1)
    a = math.sin(dlat / 2)**2 + math.cos(rad(lat1)) * math.cos(rad(lat2)) * math.sin(dlon / 2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

def calculate_speed(group):
    group.sort_values('ISO_TIME', inplace=True)
    time_diff = group['timeall'].diff().dt.total_seconds() / 3600
    distances = [0] + [haversine(lon1, lat1, lon2, lat2) for lon1, lat1, lon2, lat2 in zip(group['lonall'][:-1], group['latall'][:-1], group['lonall'][1:], group['latall'][1:])]
    group['STORM_MOVE_SPEED'] = np.divide(distances, time_diff / 24, out=np.zeros_like(distances), where=time_diff!=0)
    group['STORM_MOVE_SPEED'] = group['STORM_MOVE_SPEED'].fillna(0).abs()  # fill NaNs with 0 and take absolute value
    return group

# Data loading and preprocessing
df = pd.read_csv("ibtracs.NA.list.v04r00.csv", usecols=['SID','ISO_TIME', 'SEASON', 'LAT', 'LON', 'USA_LAT', 'USA_LON', 'WMO_WIND', 'USA_WIND', 'WMO_PRES', 'USA_PRES', 'USA_RMW', 'NATURE', 'STORM_SPEED', 'STORM_DIR'])
df = df[df['SEASON'] != 'Year']
df = df[df['SEASON'].astype(int) >= 1975]
df['yearall'] = df['SEASON'].astype(int)
df['ISO_TIME'] = pd.to_datetime(df['ISO_TIME'])
df['timeall'] = df['ISO_TIME']

# Handling missing values and data conversions
numeric_cols = ['LAT', 'LON', 'WMO_WIND', 'WMO_PRES', 'USA_LAT', 'USA_LON', 'USA_WIND', 'USA_PRES', 'USA_RMW', 'STORM_SPEED', 'STORM_DIR']
df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
placeholder = -9999.9
df['latall'] = df['LAT'].combine_first(df['USA_LAT']).fillna(placeholder).round(1).replace(placeholder, np.nan)
df['lonall'] = df['LON'].combine_first(df['USA_LON']).fillna(placeholder).round(1).replace(placeholder, np.nan)
df['windall'] = df['WMO_WIND'].combine_first(df['USA_WIND'])
df['presall'] = df['WMO_PRES'].combine_first(df['USA_PRES'])
df['rmaxall'] = df['USA_RMW']

# Storm movement speed calculation
df = df.groupby('SID').apply(calculate_speed)
df.reset_index(drop=True, inplace=True)
df.dropna(subset=['latall', 'lonall', 'windall', 'presall', 'rmaxall', 'STORM_SPEED', 'STORM_DIR', 'STORM_MOVE_SPEED'], inplace=True)

# Adjusting longitude values
df.loc[df['lonall'] < 0, 'lonall'] += 360

# Initialize a dictionary to store DataFrames for each SID
sid_dfs = {}

# Split subset_df3 into separate DataFrames for each SID
unique_sids = df['SID'].unique()
for SID in unique_sids:
    sid_dfs[SID] = df[df['SID'] == SID]

################################################
#               HOLLAND MODEL                  #
################################################

# Function definitions
def Basins_WMO(basin):
    if basin=='EP': #Eastern Pacific
        lat0,lat1,lon0,lon1=5,60,180,285
    if basin=='NA': #North Atlantic
        lat0,lat1,lon0,lon1=5,60,255,359
    if basin=='NI': #North Indian
        lat0,lat1,lon0,lon1=5,60,30,100
    if basin=='SI': #South Indian
        lat0,lat1,lon0,lon1=-60,-5,10,135
    if basin=='SP': #South Pacific
        lat0,lat1,lon0,lon1=-60,-5,135,240
    if basin=='WP': #Western Pacific
        lat0,lat1,lon0,lon1=5,60,100,180

    return lat0,lat1,lon0,lon1

def reset_variables(basin):
    lat0, lat1, lon0, lon1 = Basins_WMO(basin)
    res, alpha, beta_bg, SWRF, CF, tc_radius, Patm = 0.1, 0.55, 20., 0.85, 0.915, 1000., 101325.
    max_distance = tc_radius / 110.
    n_cols, n_rows = 36, 1000
    latspace = np.arange(lat0 + res / 2., lat1 + res / 2., res) if lat0 > 0 else np.arange(lat0 - res / 2., lat1 - res / 2., res)
    lonspace = np.arange(lon0 + res / 2., lon1 + res / 2., res)
    points = [(i, j) for i in latspace for j in lonspace]
    wind_field = {i: [] for i in range(len(points))}
    shadowlist = {kk: [] for kk in range(len(points))}
    tree = spatial.cKDTree(points)
    return points, wind_field, shadowlist


# Wind field calculation and data structures
basin = 'NA'
lat0, lat1, lon0, lon1 = Basins_WMO(basin)
res, alpha, beta_bg, SWRF, CF, tc_radius, Patm = 0.1, 0.55, 20., 0.85, 0.915, 1000., 101325.
max_distance = tc_radius / 110.
n_cols, n_rows = 36, 1000
latspace = np.arange(lat0 + res / 2., lat1 + res / 2., res) if lat0 > 0 else np.arange(lat0 - res / 2., lat1 - res / 2., res)
lonspace = np.arange(lon0 + res / 2., lon1 + res / 2., res)


# Initialize a list to store the final wind speed DataFrames for concatenation
final_dfs = []

# Wind field analysis loop
for SID, storm_df in sid_dfs.items():

    # Reset variables for each SID
    points, wind_field, shadowlist = reset_variables(basin='NA')  # Assuming reset_variables() is defined correctly

    tree = spatial.cKDTree(points)
    # Initialize a list to store the data
    wind_data_records = []

#loop over the different TCs
    latslice = storm_df["latall"].values.copy()
    lonslice = storm_df["lonall"].values.copy()
    windslice = storm_df["windall"].values.copy() * 0.5144  # 1 knot = 0.5144 m/s
    presslice = storm_df["presall"].values.copy()
    timeslice = storm_df["timeall"].values.copy()
    rmaxslice = storm_df["rmaxall"].values.copy() * 1.852  # 1 nautical mile = 1.852 km


    for j in range(1,len(latslice)):
      lat0,lat1,lon0,lon1,t0,t1=latslice[j-1],latslice[j],lonslice[j-1],lonslice[j],timeslice[j-1],timeslice[j]
      U10,Rmax,P=windslice[j],rmaxslice[j],presslice[j]
      dt_seconds = ((t1 - t0) / np.timedelta64(1, 's'))

      #Generate the seperate list of coastal points that are in the spyderweb
      distances, indices = tree.query((lat1,lon1),k=len(points), p=2,distance_upper_bound=max_distance)
      points_to_save=[points[indices[k]] for k in range(len(distances)) if distances[k]<max_distance]

      #Spyderweb step 1: Generate the spyderweb mesh --> predefined function!

      rlist,thetalist,xlist,ylist=hm.Generate_Spyderweb_mesh(n_cols,n_rows,tc_radius,lat0)

      latlist,lonlist=hm.Generate_Spyderweb_lonlat_coordinates(xlist,ylist,lat1,lon1)

      #Spyderweb step 2: Calculate the background wind --> predefined function!

      [bg,ubg,vbg]=hm.Compute_background_flow(lon0,lat0,lon1,lat1,dt_seconds)

      #Spyderweb step 3: Subtract the background flow from U10 (tropical cyclone's 10-meter wind speed)
      #For this, first convert U10 to surface level using the SWRF-constant
      #next, subtract a fraction alpha of the background flow.

      Usurf=(U10/SWRF)-(bg*alpha) #1-minute maximum sustained surface winds
      P_mesh=np.zeros((xlist.shape))
      Pdrop_mesh=np.zeros((xlist.shape))

      up=np.zeros((xlist.shape))
      vp=np.zeros((xlist.shape))

      #Spyderweb step 4: Calculate wind and pressure profile using the Holland model
      for l in range(1,n_rows):
          r=rlist[0][l]
          Vs,Ps=hm.Holland_model(lat1,P,Usurf,Rmax,r)
          Vs=Vs*SWRF      #Convert back to 10-min wind speed

          P_mesh[:,l].fill(Ps/100.)    #in Pa
          Pdrop_mesh[:,l].fill((Patm-Ps)/100.)   #in Pa

          beta=hm.Inflowangle(r,Rmax,lat0)

          for k in range(0,n_cols):
            ubp=alpha*(ubg*math.cos(math.radians(beta_bg))-np.sign(lat0)*vbg*math.sin(math.radians(beta_bg)))
            vbp=alpha*(vbg*math.cos(math.radians(beta_bg))+np.sign(lat0)*ubg*math.sin(math.radians(beta_bg)))

            up[k,l]=-Vs*math.sin(thetalist[:,0][k]+beta)+ubp
            vp[k,l]=-Vs*math.cos(thetalist[:,0][k]+beta)+vbp

      u10=CF*up
      v10=CF*vp
      windfield=np.sqrt(u10**2.+v10**2.)

      spy_points=[]
      wind_points=[]

      for k in range(n_cols):
        for l in range(n_rows):
          spy_points.append((latlist[k,l],lonlist[k,l]))
          wind_points.append(windfield[k,l])

      tree2=spatial.cKDTree(spy_points)
      current_time = storm_df["timeall"].iloc[j-1]

      #overlay the spyderweb grid with the regular grid
      for (lat,lon),idx in zip(points_to_save,range(len(points_to_save))):
        local_dist, local_ind = tree2.query((lat,lon),k=1, p=2,distance_upper_bound=max_distance)
        shadowlist[indices[idx]].append((wind_points[local_ind], current_time))

    for m in range(len(shadowlist)):
        if shadowlist[m] and np.max([wind_speed for wind_speed, _ in shadowlist[m]]) >= 18:
            lat, lon = points[m]
            wind_speeds = [wind_speed for wind_speed, _ in shadowlist[m]]
            times = [time for _, time in shadowlist[m]]

            # Append the data to the list
            wind_data_records.append({'SID': SID, 'lat': lat, 'lon': lon, 'time': times, 'wind_speeds': wind_speeds})

    # Create DataFrame from the list
    wind_data_df = pd.DataFrame(wind_data_records)

    # Add the current SID's DataFrame to the list for final concatenation
    final_dfs.append(wind_data_df)

# Concatenate all SID DataFrames to get the final DataFrame
final_wind_data_df = pd.concat(final_dfs, ignore_index=True)
