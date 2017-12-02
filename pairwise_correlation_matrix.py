'''
Author: Csaba Veraszto
Date: 18/05/2016
Content: This is how to read a csv file with Pandas from Python and compute correlations for scRNA sequencing data. 
This script should be run line by line. It uses pandas, but it works without the pandas package, example see below. 
pearson 	Standard correlation coefficient (default) 
kendall 	Kendall Tau correlation coefficient
spearman 	Spearman rank correlation coefficient
'''

cd '/where/your/files/are/' #Find your dataset (.csv format)

import pandas as pd #We will work with the pandas package, see pandas.pydata.org/

df1 = pd.read_csv('input_.csv', header=0, index_col=0) #import csv files with pandas as a dataframe with header and index
#df1.index
#df1.columns
'''
See the shape of your dataframe
In [49]: df1.shape
Out[49]: (45911, 168)
'''
pwcorr = df1.corr() # calculate pairwise correlation of DataFrame columns
'''
See the shape of your correlation matrix for columns
In [53]: pwcorr.shape
Out[53]: (168, 168)
'''
pwcorr.to_csv('correlations_sample_x_sample.csv') #Save your correlation matrix into a csv file.

genes = df1.transpose() #Transpose your dataframe to turn your rows into columns and save to another variable
pwcorr2 = genes.corr() #Calculate Pairwise correlation of your transposed DataFrame columns
'''
See the shape of your new correlation matrix
In [57]: pwcorr2.shape
Out[57]: (45911, 45911)
'''
pwcorr2.to_csv('correlations_gene_x_gene.csv') #Save your (second) correlation matrix into a new csv file.

#Let us draw a heatmap image of our pairwise correlation matrices. Alternative: seaborn https://seaborn.pydata.org/
import numpy as np 
import matplotlib.pyplot as plt

plt.pcolor(pwcorr)
plt.yticks(np.arange(0.5, len(pwcorr.index), 1), pwcorr.index)
plt.xticks(np.arange(0.5, len(pwcorr.columns), 1), pwcorr.columns)
plt.ion()
plt.show()

plt.pcolor(pwcorr2)
plt.yticks(np.arange(0.5, len(pwcorr2.index), 1), pwcorr2.index)
plt.xticks(np.arange(0.5, len(pwcorr2.columns), 1), pwcorr2.columns)
plt.ion()
plt.show()

#Parse your giant matrix 
df2 = pwcorr2.iloc[::100, ::100] #Take only every 100th element
df2 = pwcorr2.iloc[0:100,0:100] #Take only the first 100 element
'''
In [133]: df2.shape
Out[133]: (460, 460)
'''
plt.pcolor(df2)
plt.yticks(np.arange(0.5, len(df2.index), 1), df2.index)
plt.xticks(np.arange(0.5, len(df2.columns), 1), df2.columns)
plt.ion()
plt.show()


'''
Parsing arrays in numpy
import numpy as np
arr = np.random.randint(10, size=(100, 100))
arr [1:-1]
arr [1:-1].shape
arr [1:50:2].shape
arr [::10].shape
arr [::10]
arr [::10, ::10]
'''

#Once you have the BIG FILE, below is a (reliable) way to extract your genes' data with python without pandas on the cluster. The example below is for the gene X's transcripts(4). Remember that you can easily parse your database with pandas if you can read the file (takes roughly a day to read it into your v.memory), or do it on the cluster in the matter of minutes with a 900G q-node. If you can't use pandas there is an alternative: 

'''
Info on gene x

comp416373	 	4 transcripts of gene X row/column
comp416373_c0_seq2	17778 data1
comp416373_c0_seq3	17779 data2
comp416373_c0_seq4	17780 data3
comp416373_c0_seq8	17781 data4
header			0     data0	
'''

'''
ipython
import csv
csv.field_size_limit(1000000000) #there does seem to be a way of altering the limit, going by the module source code
cd '/where/your/files/are/'
r = csv.reader(open('correlations_gene_x_gene_names.csv', 'rb'))

for row in r:
     a = row[17780, 17781, 17782, 17783]


#Script below was taken from stackoverflow and modified to suit this example. To learn more about itertools, see https://docs.python.org/2/library/itertools.html
from itertools import islice
def consume(iterator, n): #Advance the iterator n-steps ahead. If n is none, consume entirely. # Use functions that consume iterators at C speed.
    if n is None:
        collections.deque(iterator, maxlen=0) # feed the entire iterator into a zero-length deque
    else:                                    
        next(islice(iterator, n, n), None) # advance to the empty slice starting at position n

import csv
csv.field_size_limit(1000000000) #there does seem to be a way of altering the limit, going by the module source code
cd '/where/your/files/are/'
r = csv.reader(open('correlations_gene_x_gene_names.csv', 'rb'))

with open("correlations_gene_x_gene_names.csv") as f:
    consume(f,17781) #get the row with your gene of interest
    line = next(f)

line = line.split(',') #convert a string to an array, here you run into a few(~30) bugged lines thanks to badly named genes
f = open('data4.csv', 'wb') #write your gene of interest's corr. values into another csv file. 
w = csv.writer(f, delimiter = ',')
w.writerows([x.split(',') for x in line]) 
f.close()

#Don't forget to add the blast data so you know what is what.
'''
