import math
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
datadir = 'C:\\Users\\jkurlander\\Dropbox (IDM)'

dataName = os.path.join(datadir,
                      'Malaria Team Folder\\data\\Novel VC\\ATSB\\phase 2 data\\merged_2018_07\\CDC_Merged.csv')
df = pd.read_csv(dataName)
df = df[pd.notnull(df['Date'])]
#a = df.loc[df['column'] == value]
conVillages = ['Krekrelo', 'Sirakele', 'Trekrou', 'Farabale', 'Kignele', 'Tiko', 'Sambadani']
expVillages = ['Madina', 'Nianguanabougou', 'Korea', 'Balala', 'Cissebougou', 'Balandougou', 'Trekrouba', 'Krekrelo']
populationsCon = {}
for month in range(3, 12):
    months = {}
    for village in conVillages:
        months[village] = 0
        for r in df.loc[(df['Vilage'] == village) & df['Month'].isin([month]) & df['Year'].isin([17])].iterrows():
            series = r[1].to_dict()
            for j in range(10):
                strName = str(j+1) + ' Station female'
                months[village] += int(series[strName])
    populationsCon[month] = months
minCon = [1] * 10
maxCon = [0] * 10
avgCon = [0] * 10
for village in conVillages:
    pop = [0]*10
    for month in range (3, 12):
        pop[month-3] += populationsCon[month][village]
    j = sum(pop)
    pop = [x/j for x in pop]
    for month in range (3, 13):
        if pop[month-3] < minCon[month-3]:
            minCon[month-3] = pop[month-3]
        if pop[month-3] > maxCon[month-3]:
            maxCon[month-3] = pop[month-3]
        avgCon[month-3] += pop[month-3]
    #plt.plot(range(3, 13), pop)
avgCon = [x/len(conVillages) for x in avgCon]


populationsExp = {}
for month in range(3, 12):
    months = {}
    for village in expVillages:
        months[village] = 0
        for r in df.loc[(df['Vilage'] == village) & df['Month'].isin([month]) & df['Year'].isin([17])].iterrows():
            series = r[1].to_dict()
            for j in range(10):
                strName = str(j+1) + ' Station female'
                months[village] += int(series[strName])
    populationsExp[month] = months
minExp = [1] * 10
maxExp = [0] * 10
avgExp = [0] * 10
for village in expVillages:
    pop = [0]*10
    for month in range (3, 12):
        pop[month-3] += populationsExp[month][village]
    j = sum(pop)
    pop = [x/j for x in pop]
    for month in range (3, 13):
        if pop[month-3] < minExp[month-3]:
            minExp[month-3] = pop[month-3]
        if pop[month-3] > maxExp[month-3]:
            maxExp[month-3] = pop[month-3]
        avgExp[month-3] += pop[month-3]
    #plt.plot(range(3, 13), pop)
avgExp = [x/len(expVillages) for x in avgExp]
plt.plot(range(3, 13), avgExp, label='Average')
plt.plot(range(3, 13), avgCon, label='Average')
plt.fill_between(range(3, 13), minExp, maxExp, alpha=0.5, linewidth = 0, label = 'Total Experimental Village Range')
plt.fill_between(range(3, 13), minCon, maxCon, alpha=0.3, linewidth = 0, label = 'Total Control Village Range')





plt.title('2016 Seasonality of Village Mosquito Populations (CDC)')
plt.xticks(range(3, 13), ['March', 'April', 'May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
plt.ylabel('HBR')
plt.ylabel('Proportion of total population found in this month')
plt.xlabel('Month')
plt.legend()
plt.show()