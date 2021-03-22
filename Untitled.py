#!/usr/bin/env python
# coding: utf-8

# We have been given the optimization problem
# Maximize Z = 6x_1 + 8x_2  <br>
# subject to the boundary conditions : <br>
#  1x_1 + 1x_2 <= 10 <br>
#  2x_1 + 3x_2 <= 25 <br>
#  1x_1 + 5x_2 <= 35 <br>

# In[1]:


n = 2
m = 3


# In[2]:


import numpy as np


# In[3]:


C_j = np.array([6,8])


# In[6]:


a_ij = np.array([[1,1,1,0,0],[2,3,0,1,0],[1,5,0,0,1]])


# In[8]:


X_B = np.array([[10],[25],[35]])


# In[9]:


B = np.eye(3)


# In[11]:


np.linalg.inv(B)


# In[ ]:


y = 

