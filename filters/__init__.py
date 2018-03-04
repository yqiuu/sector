import os

packageDir = os.path.dirname(os.path.abspath(__file__))
filterDict = {"B435":os.path.join(packageDir, "HST_ACS_F435W.npy"), 
              "V606":os.path.join(packageDir, "HST_ACS_F606W.npy"), 
              "i775":os.path.join(packageDir, "HST_ACS_F775W.npy"), 
              "I814":os.path.join(packageDir, "HST_ACS_F814W.npy"), 
              "z850":os.path.join(packageDir, "HST_ACS_F850LP.npy"), 
              "Y098":os.path.join(packageDir, "HST_IR_F098M.npy"), 
              "Y105":os.path.join(packageDir, "HST_IR_F105W.npy"), 
              "J125":os.path.join(packageDir, "HST_IR_F125W.npy"), 
              "H160":os.path.join(packageDir, "HST_IR_F160W.npy"), 
              "3.6":os.path.join(packageDir, "HST_IRAC_3.6.npy")}


