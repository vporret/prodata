#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Hierarchical clustering is a method for grouping similar data points together.

The agglomerative hierarchical clustering technique is one part of Hierarchical clustering:
where each point is considered as an individual cluster. 
After each iteration, the similar cluster merge with others until one cluster is formed.

The method for Agglomerative hierarchical clustering:
1) Compute proximity matrix
2) Define each data point as a cluster
3) Iterate: merge the two closest cluster together and update matrix
4) Stops when one single cluster remains


"""


# In[5]:


# We will use the data given in exercice 1. 

# Import the libraries
import numpy as np

X = np.array([[0,0],
            [0.4,0],
            [0,0.4],
            [1,1],
            [0.8,1],
            [1,0.8],
            [0.6,0.6],])

# Preview the data
print(X)



# In[6]:


# Let's plot the data points

import matplotlib.pyplot as plt

labels = range(1, 11)  
plt.figure(figsize=(10, 7))  
plt.subplots_adjust(bottom=0.1)  
plt.scatter(X[:,0],X[:,1], label='True Position')

for label, x, y in zip(labels, X[:, 0], X[:, 1]):  
    plt.annotate(
        label,
        xy=(x, y), xytext=(-3, 3),
        textcoords='offset points', ha='right', va='bottom')
plt.show()  


# In[ ]:


"""
From this plot, it becomes evident that the points forms two clusters.
One onsisting of points 1-3, and another from points 4-7.

"""


# In[7]:


# Use dendograms in hierarchical clustering

# Import the libraries

from scipy.cluster.hierarchy import dendrogram, linkage  
from matplotlib import pyplot as plt


linked = linkage(X, 'single')

labelList = range(1, 11)

plt.figure(figsize=(10, 7))  
dendrogram(linked,  
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()  


# In[ ]:


"""
The code starts by finding the two closest point in term of Euclidea distance,
then it joint the cluster formed by joining two points to the nearest cluster, etc.
The algorithm stops when all the points are joined together to form one cluster.

"""


# In[ ]:


## Agglomerative Hierarchical Clusterig 


# In[3]:


# Import the libraries

import pandas as pd
get_ipython().magic(u'matplotlib inline')


# In[8]:


# import the class for clustering 

from sklearn.cluster import AgglomerativeClustering

# Define the cluster using sklearn.cluster library
cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')  

# Returns the names of the clusters that each point belongs to
cluster.fit_predict(X)

print(cluster.labels_) 


# In[ ]:


"""
The first three points have been clustered together 
and the last 4 together.
"""


# In[9]:


# PLot the cluster 

plt.scatter(X[:,0],X[:,1], c=cluster.labels_, cmap='rainbow')


# In[10]:


# Agglomerative clustering with complete linkage

complete = AgglomerativeClustering(n_clusters=2,linkage='complete')

# Fit the dataset & predict the cluster labels
complete_pred = complete.fit_predict(X)


# In[12]:


# Dendogram for complete

linked = linkage(X, 'complete')

labelList = range(1, 11)

plt.figure(figsize=(10, 7))  
dendrogram(linked,  
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()  


# In[11]:


# Agglomerative clustering with average linkage

avg = AgglomerativeClustering(n_clusters=2,linkage='average')

# Fit the dataset & predict the cluster labels
avg_pred = avg.fit_predict(X)


# In[16]:


# Dendogram for average

linked = linkage(X, 'average')

labelList = range(1, 11)

plt.figure(figsize=(10, 7))  
dendrogram(linked,  
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show() 


# In[ ]:




