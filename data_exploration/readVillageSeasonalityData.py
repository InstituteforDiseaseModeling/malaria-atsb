import math
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
datadir = 'C:\\Users\\jkurlander\\Dropbox (IDM)'

#dataName = os.path.join(datadir,
#                      'Malaria Team Folder\\data\\Novel VC\\ATSB\\phase 2 data\\merged_2018_07\\CDC_Merged.csv')
dataName = os.path.join(datadir,
                      'Malaria Team Folder\\projects\\atsb\\Jake\'s Notes\\HLC_Merged_Plus_Intervention_Status.csv')

df = pd.read_csv(dataName)
df = df[df['Year'] == 2016]
#df = df[pd.notnull(df['Date'])]
monthNames = ['April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
conVillages = ['Krekrelo', 'Sirakele', 'Trekrou', 'Farabale', 'Kignele', 'Tiko', 'Sambadani']
expVillages = ['Madina', 'Nianguanabougou', 'Korea', 'Balala', 'Cissebougou', 'Balandougou', 'Trekrouba']

#a = df.loc[df['column'] == value]
populations = {}
for month in range(4, 13):
    months = {}
    for village in conVillages:
        months[village] = 0
        for r in df.loc[(df['Vilage'] == village) & df['Month'].isin([monthNames[month-4]])].iterrows():
            #series = r[1].to_dict()
            #for j in range(10):
            #    strName = str(j+1) + ' Station female'
            #    months[village] += int(series[strName])
            months[village]+=1
    populations[month] = months
print(str(populations))
min = [1] * 9
max = [0] * 9
avg = [0] * 9
for village in conVillages:
    pop = [0]*9
    for month in range(4, 13):
        pop[month-4] = populations[month][village]
    j = sum(pop)
    pop = [x/j for x in pop]
    for month in range (4, 13):
        if pop[month-4] < min[month-4]:
            min[month-4] = pop[month-4]
        if pop[month-4] > max[month-4]:
            max[month-4] = pop[month-4]
        avg [month-4] += pop[month-4]

    #plt.plot(range(4, 13), pop)
avg = [x/float(len(setOfVillages)) for x in avg]
plt.plot(range(4, 13), avg, label='Average')
plt.fill_between(range(4, 13), min, max, alpha=0.4, linewidth = 0, label = 'Total Control Range')











populationsExp = {}
for month in range(4, 13):
    months = {}
    for village in expVillages:
        months[village] = 0
        for r in df.loc[(df['Vilage'] == village) & df['Month'].isin([monthNames[month-4]])].iterrows():
            #series = r[1].to_dict()
            #for j in range(10):
            #    strName = str(j+1) + ' Station female'
            #    months[village] += int(series[strName])
            months[village]+=1
    populationsExp[month] = months
minExp = [1] * 9
maxExp = [0] * 9
avgExp = [0] * 9
for village in expVillages:
    pop = [0]*9
    for month in range(4, 13):
        pop[month-4] = populationsExp[month][village]
    j = sum(pop)
    pop = [x/j for x in pop]
    for month in range (4, 13):
        if pop[month-4] < min[month-4]:
            min[month-4] = pop[month-4]
        if pop[month-4] > max[month-4]:
            max[month-4] = pop[month-4]
        avg[month-4] += pop[month-4]

    #plt.plot(range(4, 13), pop)
avgExp = [x/float(len(setOfVillages)) for x in avgExp]
plt.plot(range(4, 13), avgExp, label='Average')
plt.fill_between(range(4, 13), minExp, maxExp, alpha=0.4, linewidth = 0, label = 'Total Range Exp')









plt.legend()
#plt.legend(df.Vilage.unique())
plt.xlabel('Month')
plt.ylabel('% of total population found in this month')
plt.title('Seasonality of Control Village Mosquito Populations (HLC)')
plt.show()