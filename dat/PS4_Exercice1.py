#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Advanced Data Analytics - Problem Set 3


# In[ ]:


#### Team Members: Anais Berney, Louis Droz, Jorgen Mellem


# In[ ]:


## Exercice 1


# In[ ]:


"""
K-Means is part of the partition-based clustering methods.

The K-mean method identifies k number of centroids, 
and allocates every point to the nearest cluster, while keeping
the centroids as small as possible.

The 'means' in the K-means refers to averaging of the data; 
thus finding the centroid.

"""


# In[22]:


# We will use the data given in the exercice set. 

# Import the libraries
import pandas as pd
import numpy as np


data = pd.DataFrame(
    data={
        'X_value':[0,0.4,0,1,0.8,1,0.6],
        'Y_value':[0,0,0.4,1,1,0.8,0.6]
    }
)

# Create a pandas DataFrame with the names and (x, y) coordinates
X_value = data.iloc[:, 0]
Y_value = data.iloc[:, 1]


# Preview the data
print(data)



# In[25]:


"""
Using K-Means

First of all, we need to specify the number of K (number of clusters)
we want out of the data.

We use K = 2.

""" 


# In[19]:


# Initialise centroids

# Centroids are centers

c1 = (1.0, 0.8)
c2 = (0.6, 0.6)



# In[23]:


# Write a function to calculate the Euclidean distances between points en centroids

def calculate_distance(centroid, X, Y):
    distances = []
        
    # Unpack the x and y coordinates of the centroid
    c_x, c_y = centroid
        
    # Iterate over the data points and calculate the distance using the           # given formula
    for x, y in list(zip(X, Y)):
        root_diff_x = (x - c_x) ** 2
        root_diff_y = (y - c_y) ** 2
        distance = np.sqrt(root_diff_x + root_diff_y)
        distances.append(distance)
        
    return distances


# In[24]:


# Apply this function to the data points and assign results in DataFrame

data['C1_Distance'] = calculate_distance(c1, data.X_value, data.Y_value)
data['C2_Distance'] = calculate_distance(c2, data.X_value, data.Y_value)

# Preview the data
print(data.head())


# In[25]:


# Compare the distances and take the smallest ones
# The centroid with the smallest distance get assigns as the cluster for that data point

# Get the min distance centroids
data['Cluster'] = data[['C1_Distance', 'C2_Distance']].apply(np.argmin, axis =1)
    
# Map the centroids accordingly and rename them
data['Cluster'] = data['Cluster'].map({'C1_Distance': 'C1', 'C2_Distance': 'C2'})
    
# Get a preview of the data
print(data.head(10))


# In[26]:


# Update the centroids by determining the mean (K-means)


# Calculate the coordinates of the new centroids from cluster 1

x_new_centroid1 = data[data['Cluster']=='C1']['X_value'].mean()
y_new_centroid1 = data[data['Cluster']=='C1']['Y_value'].mean()

# Calculate the coordinates of the new centroid from cluster 2

x_new_centroid2 = data[data['Cluster']=='C2']['X_value'].mean()
y_new_centroid2 = data[data['Cluster']=='C2']['Y_value'].mean()


# Print the coordinates of the new centroids

print('Centroid 1 ({}, {})'.format(x_new_centroid1, y_new_centroid1))
print('Centroid 2 ({}, {})'.format(x_new_centroid2, y_new_centroid2))


# In[ ]:


""" 
The K-means is a partition-based method of clustering.
The main objective of this method is to partition n observations into k clusters, in which each obersvation
belongs to the cluster with the nearest mean.

We observed that when comparing the distance between each x(i) points and the centroids, 
both C1 and C2 appears the same number of times. 

"""


# In[27]:


# Generate plot
import matplotlib.pyplot as plt

plt.scatter(X_value, Y_value, s =50, c='b')
plt.scatter(0.933333333333, 0.933333333333, s=200, c='g', marker='s') # C1
plt.scatter(0.25, 0.25, s=200, c='r', marker='s')# C2

plt.title('Data points and cluster centroids')
plt.show()


# In[ ]:


"""
In this method, we first used x6 and x7 as centroids:
x6 = c1 = (1.0, 0.8)
x7 = c2 = (0.6, 0.6)

which have been used as beginning point for every cluster.
Then the program performed iterative calculation to optimize the position of the centroids.

The k-means algorithm stops when:
- The centroids have been stabilized: no change in their values
- The defined number of iteration has been achieved (which was 10 in this example)

We observe that the centroids based on K-means algorithm: 
Centroid 1 (0.933333333333, 0.933333333333)
Centroid 2 (0.25, 0.25)
are different from the randomly selected one. 

"""

