# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 10:19:54 2020

@author: bwang
"""
import math
import numpy as np
def Readparm(compparm):
    lines=[pl.rstrip() for pl in open(compparm,'r').readlines()]
    n=lines.index("%FLAG POINTERS")
    Natom=int(lines[n+2].split()[0])
    Ntype=int(lines[n+2].split()[1])
    Nco=Ntype*Ntype;Nvdwtype=(Ntype*(Ntype+1))/2
    Nres=int(lines[n+3].split()[1])
    m=lines.index('%FLAG ATOM_NAME');n=lines.index('%FLAG CHARGE')
    line=lines[m+2:n];pl='';atomname=[]
    for pline in line:
        pl=pl+pline;
        for iii in range(80-len(pline)):pl=pl+' '
    for i in range(Natom): atomname.append(pl[i*4:i*4+4].strip())
    m=lines.index('%FLAG CHARGE');n=lines.index('%FLAG ATOMIC_NUMBER')
    line=lines[m+2:n];pl='';charge=[]
    for pline in line:pl=pl+pline
    for i in range(Natom): charge.append(float(pl[i*16:i*16+16]))
    m=lines.index('%FLAG RESIDUE_LABEL');n=lines.index('%FLAG RESIDUE_POINTER')
    line=lines[m+2:n];pl='';resname=[]
    for pline in line: pl=pl+pline+' '
    for i in range(Nres): resname.append(pl[i*4:i*4+4].strip())
    m=lines.index('%FLAG RESIDUE_POINTER');n=lines.index('%FLAG BOND_FORCE_CONSTANT')
    line=lines[m+2:n];pl='';tmp=[];resid=[0 for i in range(Natom)]
    for pline in line: pl=pl+pline
    for i in range(Nres): tmp.append(int(pl[i*8:i*8+8].strip()))
    for i in range(Nres-1):
        for j in range(tmp[i]-1,tmp[i+1]-1): resid[j]=i
    for j in range(tmp[Nres-1]-1,Natom): resid[j]=Nres-1
    m=lines.index('%FLAG ATOM_TYPE_INDEX');n=lines.index('%FLAG NUMBER_EXCLUDED_ATOMS')
    line=lines[m+2:n];pl='';iac=[]
    for pline in line: pl=pl+pline
    for i in range(Natom): iac.append(int(pl[i*8:i*8+8].strip()))
    m=lines.index('%FLAG NUMBER_EXCLUDED_ATOMS');n=lines.index('%FLAG NONBONDED_PARM_INDEX')
    line=lines[m+2:n];pl='';enum=[]
    for pline in line: pl=pl+pline
    enum_total=0
    for i in range(Natom): enum.append(pl[i*8:i*8+8].strip());enum_total=enum_total+int(enum[i])
    enum=map(int,enum)
    m=lines.index('%FLAG DIHEDRALS_INC_HYDROGEN');n=lines.index('%FLAG DIHEDRALS_WITHOUT_HYDROGEN')
    lline=lines[m+2:n];endihe=[]
    if len(lline)>1:
        for line in lline:
            line=line.strip().split()
            part1=line[:4]
            if float(part1[2])<0 or float(part1[1])<0:
                pass
            else:
                pair11_4=[]
                pair11_4.append(int(abs(int(part1[0]))/3+1));pair11_4.append(int(abs(int(part1[3]))/3+1))
                endihe.append(pair11_4)
            if len(line)>5:
                part2=line[5:9]
                if float(part2[2])<0 or float(part2[1])<0:
                    pass
                else:
                    pair21_4=[]
                    pair21_4.append(int(abs(int(part2[0]))/3+1));pair21_4.append(int(abs(int(part2[3]))/3+1))
                    endihe.append(pair21_4)
                   
    m=lines.index('%FLAG NONBONDED_PARM_INDEX');n=lines.index('%FLAG RESIDUE_LABEL')
    line=lines[m+2:n];pl='';ico=[]
    for pline in line: pl=pl+pline
    for i in range(Nco): ico.append(int(pl[i*8:i*8+8].strip()))
    m=lines.index('%FLAG LENNARD_JONES_ACOEF');n=lines.index('%FLAG LENNARD_JONES_BCOEF')
    line=lines[m+2:n];pl='';cna=[]
    for pline in line: pl=pl+pline
    for i in range(int(Nvdwtype)): cna.append(float(pl[i*16:i*16+16].strip()))
    m=lines.index('%FLAG LENNARD_JONES_BCOEF');n=lines.index('%FLAG BONDS_INC_HYDROGEN')
    line=lines[m+2:n];pl='';cnb=[]
    for pline in line: pl=pl+pline
    for i in range(int(Nvdwtype)): cnb.append(float(pl[i*16:i*16+16].strip()))
    return iac,ico,cna,cnb,Ntype,charge,Natom,Nres,resname,resid,atomname,endihe

def Readtrj(Natom,filename):
    lines=open(filename,'r').readlines()
    lines=lines[1:]
    num=int(math.ceil(Natom*3/10.)+1)
    #print('num:',num)
    Nframe=len(lines)/num
    #print('Nframe:',Nframe)
    xcrd=[];ycrd=[];zcrd=[]
    crd1=[]
    rows=math.ceil(Natom*float(3)/10.0)
    a=int(rows)+1
    for fr in range(int(Nframe)):
        num0=int(fr*a)
        num1=int((fr+1)*a-1)
        xyz=lines[num0:num1]    
        for line in xyz:
            crd1+=line.strip().split()
    coord=np.array(crd1).reshape(-1,3)
    x=coord[:,0];y=coord[:,1];z=coord[:,2]
    xcrd.extend(x)
    ycrd.extend(y)
    zcrd.extend(z)        
    return xcrd,ycrd,zcrd
