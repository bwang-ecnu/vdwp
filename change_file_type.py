# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 10:20:51 2020

@author: bwang
"""
import os,re
import numpy as np
"""
change pdb file type of amber and rosetta
amber to rosetta : AtoR 
rosetta to amber : RtoA
"""
def AtoR(sourcepath,rosettapath,fileAR,typelines):
        """
        sourcepath : amber pdb file path
        rosettapath : changed rosetta pdb file path
        fileAR : atom type file
        typelines : read pdb file in a list
        """
        lines=open(sourcepath+'/'+str(fileAR),'r').readlines()
        Lines=[]
        Changepdbname=os.path.join(rosettapath,fileAR)
        mutfile=open(Changepdbname,'w')
        for line in lines:
            if line.strip().split()[0]=='ATOM':
                if line.strip().split()[3]=='CYX':
                    Line=line.replace('CYX','CYS')
                    Lines.append(Line)
                elif line.strip().split()[3]=='HID':
                    Line=line.replace('HID','HIS')
                    Lines.append(Line)
                elif line.strip().split()[3]=='HIP':
                    Line=line.replace('HIP','HIS')
                    Lines.append(Line)
                elif line.strip().split()[3]=='HIE':
                    Line=line.replace('HIE','HIS')
                    Lines.append(Line)
                else:
                    Lines.append(line)
        for line in Lines:
            if str(line[:4])=='ATOM':
                resname=line[17:20]
                if 'H' in line[12:16]:
                    for ll in typelines:
                        ll=ll.strip().split()
                        if str(line[12:16].strip())==str(ll[1]) and resname==str(ll[3]):
                            name1=line[12:16].strip()
                            name2=ll[2]
                            if str(name1)=='H':
                                mutfile.write(line[:56] +'\n')
                            elif len(name1)==2:
                                if name1[1] in ['1','2','3']:
                                    line=line[:12]+str(name2)+'  '+line[16:56]
                                    mutfile.write(line[:56] +'\n')
                                elif name1[1] not in ['1','2','3']:
                                    mutfile.write(line[:56] +'\n')
                            elif len(name1)==3:
                                if 'A' in str(name1) or 'B' in str(name1):
                                    line=line[:12]+str(name2)+' '+line[16:56]
                                    mutfile.write(line[:56] +'\n')
                                else:
                                    mutfile.write(line[:56] +'\n')
                            elif len(name1)==4:
                                line=line[:12]+str(name2)+line[16:56]
                                mutfile.write(line[:56] +'\n')
                else:
                    mutfile.write(line[:56] +'\n')
        mutfile.write('TER' +'\n')
        
def RtoA(rosettapath,amberpath,fileAR,typelines):
        """
        rosettapath : rosetta pdb file path
        amberpath : changed amber pdb file path
        fileAR : atom type file
        typelines : read pdb file in a list
        """
        rosettapdb=os.path.join(rosettapath,fileAR)
        print(rosettapdb)
        LL=[];Lines1=[];lines=[]
        RESIDname=[];RESIDno=[]
        Lines1=open(rosettapdb+'_0001.pdb','r').readlines()
        for line in Lines1:
            if 'ATOM' in line[:4]:
                lines.append(line)
        amberpdb=os.path.join(amberpath,fileAR)
        mutfile=open(amberpdb,'w')
        for line in lines:
            if str(line[:4])=='ATOM':
                resname=line[17:20]
                if 'H' in line[12:16].strip():
                    for typeline in typelines:
                        tll=typeline.strip().split()
                        if str(line[12:16].strip())==str(tll[2]) and str(tll[3])==str(resname):
                            name1=tll[2].strip()
                            name2=tll[1]
                            if str(name1)=='H':
                                LL.append(line)
                            elif len(name1)==2:
                                line=line[:12]+' '+str(name2)+'  '+line[17:]
                                LL.append(line)
                            elif len(name1)==3:
                                line=line[:12]+' '+str(name2)+' '+line[17:]
                                LL.append(line)
                            elif len(name1)==4:
                                line=line[:12]+str(name2)+' '+line[17:]
                                LL.append(line)
                else:
                    LL.append(line) 
        for line in LL:
            RESIDname.append(line[17:21].strip())
            lll=re.sub("\D","",line[23:28].strip())
            RESIDno.append(int(lll))
        his=[]
        for i in range(len(RESIDname)):
            if str(RESIDname[i])=='HIS':
                his.append(RESIDno[i])
        his=list(np.unique(his))
        his_name=[]
        for i in range(len(his)):
            flag1=0;flag2=0;flag3=0
            for line in LL:
                if str(line[23:27].strip())==str(his[i]) and str(line[12:17].strip())=='HE1':
                    flag1=1
                if str(line[23:27].strip())==str(his[i]) and str(line[12:17].strip())=='HE2':
                    flag2=1
                if str(line[23:27].strip())==str(his[i]) and str(line[12:17].strip())=='HD1':
                    flag3=1
            if flag1==1 and flag2==1 and flag3==1:
                his_name.append('HIP')
            elif flag1==1 and flag2==0 and flag3==1:
                his_name.append('HID')
            elif flag1==1 and flag2==1 and flag3==0:
                his_name.append('HIE')
        for line in LL:
            Ll=line.strip().split()
            Nno=re.sub("\D","",str(Ll))
            if int(Nno) in his:
                for i in range(len(his)):
                    if str(Ll[5])==str(his[i]):
                        Line=line.replace('HIS',str(his_name[i]))
                        mutfile.write(Line)
            else:
                mutfile.write(line)