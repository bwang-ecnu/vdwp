import os,re
import numpy as np
from sys import argv
from multiprocessing import Pool
from multiprocessing import Manager
from multiprocessing import Process
from change_file_type import AtoR,RtoA
from tools import Readtrj,Readparm

class vdwid(object):
    """This is the main class for vdwP prediction
    top,crd,expttxt,mutaframes is the Parmtop file, MD crd file, mutation file and the frames you want to use 
    """
    def __init__(self,top,crd,mutainfo,frames,prod):
        self.top = top #'./example/1stn.parm7'  
        self.crd = crd #'./example/md.crd'
        self.expttxt = mutainfo #'expt.txt'
        self.mainpath = os.getcwd()
        self.mutaframes = int(frames) #10
        self.parall = True
        self.prod = int(prod)
        
    def _findfile(self,path,label):
        """find files that have the same suffix
        """
        files=[];filenames=os.listdir(path)
        for name in filenames:
            if os.path.splitext(name)[0]==str(label):
                files.append(name)
        return files
    
    def main(self):
        """Whether to use parallel computing
        """
        lines=open(self.expttxt,'r').readlines()
        if self.parall:
            pool=Pool(self.prod)
            pool.map(self.mutatraj,lines)
            pool.close()
            pool.join()
        else:
            for line in lines:
                self.mutatraj(line)

    def minimize(self,amberpdbpath,pdbfile,amberminpath):
        """Minimize the mutation frame
        """
        minfile=os.path.join(self.mainpath,'min.in')
        res='';pdbname=os.path.join(amberpdbpath,pdbfile);leap=os.path.join(amberpdbpath,'leap.in')
        top=os.path.join(amberpdbpath,pdbfile+'.prmtop');crd=os.path.join(amberpdbpath,pdbfile+'.inpcrd')
        rst=os.path.join(amberpdbpath,pdbfile+'min.rst');minout=os.path.join(amberpdbpath,'min.out')
        res+='source leaprc.protein.ff14SB' +'\n'
        res+='source leaprc.gaff' +'\n'
        res+='source leaprc.water.tip3p' +'\n'
        res+='com=loadpdb '+pdbname +'\n'
        res+='setbox com "vdw"' +'\n'
        res+='addIons com Na+ 0' +'\n'
        res+='addIons com Cl- 0'+'\n'
        res+='saveamberparm com '+top+' '+crd +'\n'
        res+='quit'
        f=open(leap,'w');f.write(res);f.close()
        os.system('tleap -s -f '+leap)
        os.system('pmemd -O -i '+minfile+' -p '+top+' -c '+crd+' -r '+rst+' -ref '+pdbname+'.inpcrd -o '+minout)
        res='';res+='parm '+top+'\n';res+='trajin '+rst+'\n';res+='strip :WAT,Na+,Cl- \n';res+='trajout '+rst+'stri \n'
        f=open(leap,'w');f.write(res);f.close()
        os.system('cpptraj < '+leap) 
        res='';res+='parm '+top+'\n';res+='parmstrip :WAT,Na+,Cl- \n';res+='parmwrite out '+top+'stri \n'
        f=open(leap,'w');f.write(res);f.close()
        os.system('cpptraj < '+leap)        
        os.system('ambpdb -p '+top+'stri -c '+rst+'stri > '+amberminpath+'/'+pdbfile+'.pdb')
        
    def mutatraj(self,line):
        """use rosetta program fixbb to make mutation frame
        """
        line=line.strip().split()
        no=line[1];mutares=line[2]
        typelines=open(self.mainpath+'/A-R.dat','r').readlines()  
        mutapath=os.path.join(self.mainpath,no+mutares)
        if not os.path.exists(mutapath):os.makedirs(mutapath)
        rosettapdb=os.path.join(mutapath,'rospdbs')
        rosettamuta=os.path.join(mutapath,'rosmutapdbs')
        amberpdb=os.path.join(mutapath,'ambpdbs')
        minpdb=os.path.join(mutapath,'minpdbs')
        if not os.path.exists(rosettapdb):os.makedirs(rosettapdb)
        if not os.path.exists(rosettamuta):os.makedirs(rosettamuta)
        if not os.path.exists(amberpdb):os.makedirs(amberpdb)
        if not os.path.exists(minpdb):os.makedirs(minpdb)
        res=''
        res+='parm '+self.mainpath+'/'+str(self.top)+'\n';res+='trajin '+self.mainpath+'/'+self.crd+' 1 '+str(self.mutaframes)+'\n'
        res+='trajout snap.pdb multi chainid A\n'
        f=open(mutapath+'/traj.in','w')
        f.write(res)
        f.close()
        os.system('cpptraj < '+no+mutares+'/traj.in')
        res=''
        res+='NATRO'+'\n';res+='start'+'\n'
        res+='    '+str(no)+' A '+' PIKAA '+str(mutares)+'\n'
        resfile=mutapath+'/resfile.in'
        f=open(resfile,'w')
        f.write(res)
        f.close()
        inputpdbs=self._findfile('./','snap.pdb')
        print(inputpdbs)
        for i in range(len(inputpdbs)):
            AtoR(self.mainpath,rosettapdb,inputpdbs[i],typelines)
            inputf=os.path.join(rosettapdb,inputpdbs[i])
            fileno=re.sub("\D","",inputpdbs[i])
            cmd='fixbb.default.linuxgccrelease -in:file:s '+inputf+' -out:path:all '+rosettamuta+' -resfile '+resfile+' -out:suffix .'+str(fileno)+' -overwrite -ex1:level 2 -ex2:level 2 -ex3:level 0 -ex4:level 0' 
            os.system(cmd)
            RtoA(rosettamuta,amberpdb,inputpdbs[i],typelines)
            self.minimize(amberpdb,inputpdbs[i],minpdb)
        self.calint(no,mutares,'wild')
        self.calint(no,mutares,'muta')
        
    def calint(self,no,mutares,choice):
        """calculate Van der Waals interaction energy of the wild and mutation frame 
        """
        if choice=='wild':
            iac,ico,cna,cnb,Ntype,charge,Natom,Nres,resname,resid,atomname,endihe=Readparm(self.top)
            xcrd,ycrd,zcrd=Readtrj(Natom,self.crd)
            f1=open(self.mainpath+'/'+str(no)+'vdw.txt','w')
        elif choice=='muta':
            f2=open(self.mainpath+'/'+no+mutares+'/'+str(no)+'vdw.txt','w')
            mutatop=os.path.join(self.mainpath,no+mutares,'ambpdbs','snap.pdb.1.prmtopstri')
            iac,ico,cna,cnb,Ntype,charge,Natom,Nres,resname,resid,atomname,endihe=Readparm(mutatop)
            pdbfiles=os.listdir(os.path.join(self.mainpath,no+mutares,'minpdbs'))
            xcrd=[];ycrd=[];zcrd=[]
            for pdb in pdbfiles:
                pdb=os.path.join(self.mainpath,no+mutares,'minpdbs',pdb)
                lines=open(pdb,'r').readlines()
                for line in lines:
                    if 'ATOM' in line[:4]:
                        xcrd.append(float(line[29:38].strip()))
                        ycrd.append(float(line[38:46].strip()))
                        zcrd.append(float(line[46:54].strip()))        
        vdwlist=[]
        for frame in range(self.mutaframes):
            Ele=0;Vdw=0
            for i in range(Natom):
                if resid[i]==int(no)-1: 
                    xi=float(xcrd[Natom*frame+i]);yi=float(ycrd[Natom*frame+i]);zi=float(zcrd[Natom*frame+i])
                    for j in range(Natom):
                        if resid[j] not in [int(no)-2,int(no)-1,int(no)]:
                            pairi_j=[]
                            pairi_j.append(i+1);pairi_j.append(j+1)
                            xj=float(xcrd[Natom*frame+j]);yj=float(ycrd[Natom*frame+j]);zj=float(zcrd[Natom*frame+j])
                            dis=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                            ele=charge[i]*charge[j]/np.sqrt(dis)
                            r6=dis*dis*dis;r12=r6*r6
                            n1=iac[i];n2=iac[j];n3=Ntype*(n1-1)+n2
                            a=cna[ico[n3-1]-1];b=cnb[ico[n3-1]-1]
                            vdw=(a/r12)-(b/r6)
                            #if pairi_j not in endihe:
                            #    pass
                            #else:
                            #    ele=ele/1.2
                            #    vdw=vdw/2.0
                            Ele+=ele;Vdw+=vdw
            print(Vdw)
            vdwlist.append(Vdw) 
        if choice=='wild':
            for vdwv in vdwlist:
                f1.write(str(vdwv)+'\n')
            f1.close()
        elif choice=='muta':
            for vdwv in vdwlist:
                f2.write(str(vdwv)+'\n')
            f2.close()                           

if __name__=='__main__':
    print('###########Welcome To vdwP################')
    top=argv[1];crd=argv[2];mutainfo=argv[3];frames=argv[4];prod=argv[5]
    VDWID=vdwid(top,crd,mutainfo,frames,prod)
    VDWID.main()
