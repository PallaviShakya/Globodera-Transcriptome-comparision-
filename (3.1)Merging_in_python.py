
# coding: utf-8

# In[13]:


import os
os.getcwd()


# In[5]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.core.display import display, HTML, Latex, Markdown

# Set the visualization settings
matplotlib.rcParams['axes.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)


# ## Load the data

# In[6]:


pal_on_rosto = pd.read_csv("pal_on_rosto2.csv")
rosto_on_pal = pd.read_csv("rosto_on_pa_rol2.csv")

display(pal_on_rosto.head())
display(rosto_on_pal.head())

display(pal_on_rosto.columns)
display(rosto_on_pal.columns)


# ## Data Analysis for determining thresholds

# In[8]:





sns.distplot(rosto_on_pal[rosto_on_pal.maxhit<0.001].maxhit)
print(rosto_on_pal.maxhit.median())
print(rosto_on_pal.maxhit.mean())


# In[10]:



sns.distplot(pal_on_rosto[pal_on_rosto.maxhit<0.001].maxhit)
print(pal_on_rosto.maxhit.median())
print(pal_on_rosto.maxhit.mean())
#sns.distplot(pal_on_rosto.maxhit)


# ## Filtering the data to thresholds

# In[12]:



rosto_on_pal_filtered = rosto_on_pal[(rosto_on_pal.maxhit<0.001)]
pal_on_rosto_filtered = pal_on_rosto[(pal_on_rosto.maxhit<0.001)]

print(len(rosto_on_pal_filtered))
print(len(pal_on_rosto_filtered))


# ## merging the filtered data

# In[ ]:


# compuationally intensive

pal_palid_merged = pd.merge(pal_on_rosto_filtered, rosto_on_pal_filtered, left_on='palid', right_on='pal_hit', suffixes=('_pal_on_rosto', '_rosto_n_pal'))
pal_palid_merged.to_csv('pal_palid_merged.csv', index=False)


# # Final step of Filtering in Python (Just a script that I did not run) 

# In[ ]:


pal_id = []

for line in pal_palid_merged:
    # filter(firsttry, palid==pal_hit) %>%
    if palid != pal_hit: # so if not equal to each other we skip the whole line
        continue
 

    else: # we continue
        palid.split(palid("."))
        # separate(palid, into = c("palid", "bla"), sep = "\\.") %>%
        # dplyr::select(-bla) %>%
        #remove the transcript number  (e.g. Gpal_D383_g0001.t1 becomes Gpal_D383_g0001) 
        pal_id.append(pal_id) #store gene_id in a list 
        sorted(pal_id) # Alphabetical sort lines by palId


# pal_palid_merged['gene_id'] = pal_palid_merged['palid']+"test"
pal_palid_merged = pal_palid_merged[pal_palid_merged['palid']!=pal_palid_merged['pal_hit']]  # filter(firsttry, palid==pal_hit) %>%
pal_palid_merged['gene_id'] = pal_palid_merged.apply(lambda row: row['palid'].split('.')[0], axis=1)

pal_palid_merged['maxhit'].max()
pal_palid_merged['maxhit'].min()

list_of_gene_ids = pd.unique(pal_palid_merged['gene_id'])

#Gpal_D383_g0001.t1. > [Gpal_D383_g0001, t1, ""]
    
#mutate(palbest= min(maxhit)) %>%
#identify the best rosto blast hit
        for p_id in pal_id:
            palbest = min(maxhit_pal_on_rosto)
            rostobest = min(maxhit_rosto_n_pal)
            palid = list(dict.fromkeys(palid)) #Remove any duplicate
        

display(pal_palid_merged.head)

 



    

 

