import os
import sys
import time

fileList = []
rootdir = raw_input("Enter root directory: ")
#save=open(raw_input("Enter file name to save & press enter: "),'w')

for each in open("filenames_hpm.txt"):
	f=each.strip().split("\t")
	oldfilename,newdir,newfilename,olddir=f[0],f[1],f[2],f[3] #change according to filecolumns
	for root, subFolders, files in os.walk(rootdir):
		for f in subFolders:
			print "current ",os.path.basename(f) 
			print "old", olddir.split("/")[-1]
			print "olddir", os.path.join(root,os.path.basename(f))
			print "newdir",os.path.join(root,newdir)
			if  os.path.basename(f)==olddir.split("/")[-1]:
				
				os.rename(os.path.join(root,os.path.basename(f)),os.path.join(root,newdir))
			
		
for each in open("filenames_hpm.txt"):
	f=each.strip().split("\t")
	oldfilename,newdir,newfilename,olddir=f[0],f[1],f[2],f[3] #change according to filecolumns
	for root, subFolders, files in os.walk(rootdir):
		for f in files:
			#print os.path.join(root,f),os.path.dirname(olddir)
			if  oldfilename.split("/")[-1]==f:
				os.rename(os.path.join(root,f),os.path.join(root,newfilename))				

		
				
		
				
			


	

