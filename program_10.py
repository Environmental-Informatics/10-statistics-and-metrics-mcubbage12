#!/bin/env python
# Created on March 25, 2020
#  by Keith Cherkauer
#edited by Marissa Cubbage 4/10/2019
#
# This script servesa as the solution set for assignment-10 on descriptive
# statistics and environmental informatics.  See the assignment documention 
# and repository at:
# https://github.com/Environmental-Informatics/assignment-10.git for more
# details about the assignment.


#import packages
import pandas as pd
import scipy.stats as stats
import numpy as np

#function that reads file and completes gross data check. This function also returns the number os missing data values
def ReadData( fileName ):
    """This function takes a filename as input, and returns a dataframe with
    raw data read from that file in a Pandas DataFrame.  The DataFrame index
    should be the year, month and day of the observation.  DataFrame headers
    should be "agency_cd", "site_no", "Date", "Discharge", "Quality". The 
    "Date" column should be used as the DataFrame index. The pandas read_csv
    function will automatically replace missing values with np.NaN, but needs
    help identifying other flags used by the USGS to indicate no data is 
    availabiel.  Function returns the completed DataFrame, and a dictionary 
    designed to contain all missing value counts that is initialized with
    days missing between the first and last date of the file."""
    #global DataDF
    #global MissingValues
    # define column names
    colNames = ['agency_cd', 'site_no', 'Date', 'Discharge', 'Quality']
  
    # open and read the file
    DataDF = pd.read_csv(fileName, header=1, names=colNames,  
                         delimiter=r"\s+",parse_dates=[2], comment='#',
                         na_values=['Eqp'])
    DataDF = DataDF.set_index('Date')
    
    #gross data check
    for i in range (0,len(DataDF)-1):
          if 0 > DataDF['Discharge'].iloc[i]:
             DataDF['Discharge'].iloc[i]=np.nan
          
    
    # quantify the number of missing values
    MissingValues = DataDF["Discharge"].isna().sum()
    
    return( DataDF, MissingValues )

def ClipData( DataDF, startDate, endDate ):
    """This function clips the given time series dataframe to a given range 
    of dates. Function returns the clipped dataframe and and the number of 
    missing values."""
    DataDF=DataDF[startDate:endDate]
    MissingValues = DataDF["Discharge"].isna().sum()

    
    return( DataDF, MissingValues )

def CalcTqmean(Qvalues):
    """This function computes the Tqmean of a series of data, typically
       a 1 year time series of streamflow, after filtering out NoData
       values.  Tqmean is the fraction of time that daily streamflow
       exceeds mean streamflow for each year. Tqmean is based on the
       duration rather than the volume of streamflow. The routine returns
       the Tqmean value for the given data array."""
    a=Qvalues.dropna()
    Tqmean= ((a>a.mean()).sum()/len(a))
    return ( Tqmean )

def CalcRBindex(Qvalues):
    """This function computes the Richards-Baker Flashiness Index
       (R-B Index) of an array of values, typically a 1 year time
       series of streamflow, after filtering out the NoData values.
       The index is calculated by dividing the sum of the absolute
       values of day-to-day changes in daily discharge volumes
       (pathlength) by total discharge volumes for each year. The
       routine returns the RBindex value for the given data array."""
    
    a=Qvalues.dropna()
    RBindex=((abs(a.diff().dropna())).sum())/(a.sum())
    return ( RBindex )

def Calc7Q(Qvalues):
    """This function computes the seven day low flow of an array of 
       values, typically a 1 year time series of streamflow, after 
       filtering out the NoData values. The index is calculated by 
       computing a 7-day moving average for the annual dataset, and 
       picking the lowest average flow in any 7-day period during
       that year.  The routine returns the 7Q (7-day low flow) value
       for the given data array."""
   
    Qvalues=Qvalues.dropna()
    val7Q=(Qvalues.rolling(7).mean()).min()
        
        
    return ( val7Q )

def CalcExceed3TimesMedian(Qvalues):
    """This function computes the number of days with flows greater 
       than 3 times the annual median flow. The index is calculated by 
       computing the median flow from the given dataset (or using the value
       provided) and then counting the number of days with flow greater than 
       3 times that value.   The routine returns the count of events greater 
       than 3 times the median annual flow value for the given data array."""
    Qvalues=Qvalues.dropna()
    median_year=Qvalues.median() #calculate median for each year
    median3x=(Qvalues[Qvalues>(3*median_year)].shape[0]) #determine number of days where flow was 3 times the median
    return ( median3x )

def GetAnnualStatistics(DataDF):
    """This function calculates annual descriptive statistcs and metrics for 
    the given streamflow time series.  Values are retuned as a dataframe of
    annual values for each water year.  Water year, as defined by the USGS,
    starts on October 1."""
    
    global WYDataDF
    #create empty dataframe to fill with values
    colnames=['site_no', 'Mean Flow', 'Peak Flow', 'Median Flow', 'Coeff Var', 'Skew', 'Tqmean', 'R-B Index', '7Q', '3xMedian']
    year_data=DataDF.resample('AS-OCT').mean()
    WYDataDF = pd.DataFrame(0,index=year_data.index, columns=colnames)
    
    #fill data frame with annual values for the water year, water yer starts on october 1
    Mean_Flow= DataDF.resample("AS-OCT").mean()
    Peak_Flow=DataDF.resample("AS-OCT").max()
    Median_Flow=DataDF.resample("AS-OCT").median()
    Coeff_Var=(DataDF.resample("AS-OCT").std()/DataDF.resample("AS-OCT").mean())*100.0
    Skew_R=DataDF.resample("AS-OCT").apply(stats.skew)
    WYDataDF['Tqmean']=DataDF.resample("AS-OCT").apply({'Discharge':lambda x:CalcTqmean(x)}) #use .apply and lambda when using a custom function
    WYDataDF['3xMedian']=DataDF.resample("AS-OCT").apply({'Discharge':lambda x:CalcExceed3TimesMedian(x)})
    WYDataDF['7Q']=DataDF['Discharge'].resample("AS-OCT").apply({lambda x:Calc7Q(x)})
    WYDataDF['R-B Index']=DataDF['Discharge'].resample("AS-OCT").apply({lambda x:CalcRBindex(x)})
   
    WYDataDF['site_no']=DataDF.resample('AS-OCT')['site_no'].mean()   
    WYDataDF['Mean Flow']= Mean_Flow['Discharge']
    WYDataDF['Peak Flow']= Peak_Flow['Discharge']
    WYDataDF['Median Flow']= Median_Flow['Discharge']
    WYDataDF['Coeff Var']= Coeff_Var['Discharge']
    WYDataDF['Skew']= Skew_R['Discharge']
    return ( WYDataDF )

def GetMonthlyStatistics(DataDF):
    """This function calculates monthly descriptive statistics and metrics 
    for the given streamflow time series.  Values are returned as a dataframe
    of monthly values for each year."""
    global MoDataDF
    #create dataframe to fill with values
    colnames=['site_no', 'Mean Flow', 'Coeff Var', 'Tqmean', 'R-B Index']
    month_data=DataDF.resample('MS').mean()
    MoDataDF = pd.DataFrame(0,index=month_data.index, columns=colnames)
    #fill data frame with monthly statistics
    Mean_Flow= DataDF.resample("MS").mean()
    Coeff_Var=(DataDF.resample("MS").std()/DataDF.resample("MS").mean())*100.0
    Tqmean=DataDF.resample("MS").apply(CalcTqmean)
    RBindex=DataDF.resample("MS").apply(CalcRBindex)

    MoDataDF['site_no']=DataDF.resample('MS')['site_no'].mean()   

    MoDataDF['Mean Flow']= Mean_Flow['Discharge']
    MoDataDF['Coeff Var']= Coeff_Var['Discharge']
    MoDataDF['Tqmean']= Tqmean['Discharge']
    MoDataDF['R-B Index']= RBindex['Discharge']

    return ( MoDataDF )

def GetAnnualAverages(WYDataDF):
    """This function calculates annual average values for all statistics and
    metrics.  The routine returns an array of mean values for each metric
    in the original dataframe."""
    #create series of data that holds the average annual value for each metric in WYDataDF
    AnnualAverages=WYDataDF.mean(axis=0)
    
    return( AnnualAverages )

def GetMonthlyAverages(MoDataDF):
    """This function calculates annual average monthly values for all 
    statistics and metrics.  The routine returns an array of mean values 
    for each metric in the original dataframe."""
    #create empty dataframe to fill with values
    colnames=['site_no', 'Mean Flow', 'Coeff Var', 'Tqmean', 'R-B Index']
    MonthlyAverages=pd.DataFrame(0,index=range(1,13), columns=colnames)
    
    #fill dataframe with average avlues for each metric by month
    for n in range(0,12):
        MonthlyAverages.iloc[n,0]=MoDataDF['site_no'][::12].mean()
        
    index=[(0,3),(1,4),(2,5),(3,6),(4,7),(5,8),(6,9),(7,10),(8,11),(9,0),(10,1),(11,2)]
    for(n,m) in index:
        MonthlyAverages.iloc[n,1]=MoDataDF['Mean Flow'][m::12].mean()
    for(n,m) in index:
        MonthlyAverages.iloc[n,2]=MoDataDF['Coeff Var'][m::12].mean()
    for(n,m) in index:
        MonthlyAverages.iloc[n,3]=MoDataDF['Tqmean'][m::12].mean()
    for(n,m) in index:
        MonthlyAverages.iloc[n,4]=MoDataDF['R-B Index'][m::12].mean()
        
        
    return( MonthlyAverages )

# the following condition checks whether we are running as a script, in which 
# case run the test code, otherwise functions are being imported so do not.
# put the main routines from your code after this conditional check.

if __name__ == '__main__':

    # define filenames as a dictionary
    # NOTE - you could include more than jsut the filename in a dictionary, 
    #  such as full name of the river or gaging site, units, etc. that would
    #  be used later in the program, like when plotting the data.
    fileName = { "Wildcat": "WildcatCreek_Discharge_03335000_19540601-20200315.txt",
                 "Tippe": "TippecanoeRiver_Discharge_03331500_19431001-20200315.txt" }
    
    # define blank dictionaries (these will use the same keys as fileName)
    DataDF = {}
    MissingValues = {}
    WYDataDF = {}
    MoDataDF = {}
    AnnualAverages = {}
    MonthlyAverages = {}
    
    # process input datasets
    for file in fileName.keys():
        
        print( "\n", "="*50, "\n  Working on {} \n".format(file), "="*50, "\n" )
        
        DataDF[file], MissingValues[file] = ReadData(fileName[file])
        print( "-"*50, "\n\nRaw data for {}...\n\n".format(file), DataDF[file].describe(), "\n\nMissing values: {}\n\n".format(MissingValues[file]))
        
        # clip to consistent period
        DataDF[file], MissingValues[file] = ClipData( DataDF[file], '1969-10-01', '2019-09-30' )
        print( "-"*50, "\n\nSelected period data for {}...\n\n".format(file), DataDF[file].describe(), "\n\nMissing values: {}\n\n".format(MissingValues[file]))
        
        # calculate descriptive statistics for each water year
        WYDataDF[file] = GetAnnualStatistics(DataDF[file])
        
        # calcualte the annual average for each stistic or metric
        AnnualAverages[file] = GetAnnualAverages(WYDataDF[file])
        
        print("-"*50, "\n\nSummary of water year metrics...\n\n", WYDataDF[file].describe(), "\n\nAnnual water year averages...\n\n", AnnualAverages[file])

        # calculate descriptive statistics for each month
        MoDataDF[file] = GetMonthlyStatistics(DataDF[file])

        # calculate the annual averages for each statistics on a monthly basis
        MonthlyAverages[file] = GetMonthlyAverages(MoDataDF[file])
        
        print("-"*50, "\n\nSummary of monthly metrics...\n\n", MoDataDF[file].describe(), "\n\nAnnual Monthly Averages...\n\n", MonthlyAverages[file])
        
        
#run definitions for two rivers, Tippecanoe River and Wildcat River
ReadData('TippecanoeRiver_Discharge_03331500_19431001-20200315.txt')
DataDF, MissingValues=ClipData( DataDF,'1969-10-01', '2019-09-30')
Tippe_WY=GetAnnualStatistics(DataDF)
Tippe_Mo=GetMonthlyStatistics(DataDF)
Tippe_annualave=GetAnnualAverages(Tippe_WY).to_frame()#turn series into dataframe
Tippe_monthlyave= GetMonthlyAverages(Tippe_Mo)

ReadData('WildcatCreek_Discharge_03335000_19540601-20200315.txt')
DataDF, MissingValues= ClipData( DataDF,'1969-10-01', '2019-09-30')
Wildcat_WY=GetAnnualStatistics(DataDF)
Wildcat_Mo=GetMonthlyStatistics(DataDF)
Wildcat_annualave=GetAnnualAverages(Wildcat_WY).to_frame()
Wildcat_monthlyave= GetMonthlyAverages(Tippe_Mo)

#######Output results#######

#annual metrics
Annual_stats=Tippe_WY.append(Wildcat_WY) #put together data from both rivers
Annual_stats['Station']=['Tippe' if x==3331500 else 'Wildcat' for x in Annual_stats['site_no']]#add column with river name
Annual_stats.to_csv('Annual_Metrics.csv', sep=',', index=True)#write csv

#monthly metrics
Month_stats=Tippe_Mo.append(Wildcat_Mo)
Month_stats['Station']=['Tippe' if x==3331500 else 'Wildcat' for x in Month_stats['site_no']]
Month_stats.to_csv('Monthly_Metrics.csv', sep=',',index=True)


#annual average metrics
Annual_Ave=pd.merge(Tippe_annualave, Wildcat_annualave,left_index=True, right_index=True,  how='right')
Annual_Ave1=Annual_Ave.rename(columns={'0_x':"Tippe", '0_y':"Wildcat"})#rename columns of dataframe to indicate river name
Annual_Ave1.to_csv('Average_Annual_Metrics.txt', sep='\t', index=True)

#monthly average metrics
Month_Ave=Tippe_monthlyave.append(Wildcat_monthlyave)
Month_Ave['Station']=['Tippe' if x==3331500 else 'Wildcat' for x in Month_Ave['site_no']]
Month_Ave.to_csv('Average_Monthly_Metrics.txt', sep='\t',index=True)
