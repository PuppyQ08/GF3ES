#Parse the ERI overlap core hamiltonian from g16 files generally.
import os
import pandas as pd
import sys


#XXX path of the data file
abs_cwd = "/home/xiuyiqin/GF3ES/Test_SCF/src/" #absolute current directory
abs_data = "/home/xiuyiqin/GF3ES/Test_SCF/data/" #absolute data directory
logfile = "inputfile.log"
comfile = "inputfile.com"
#read the object name
with open(abs_data+comfile) as comf:
    comparse = comf.readlines() 
for ii in range(len(comparse)):
    eachline = comparse[ii].split()
    if eachline and eachline[-1] == "calculation": #it could be empty
        objname = eachline[1] + "_"
        print("Done fetch object name")
comf.close()
#read log file
with open(abs_data+logfile) as f:
    parse  = f.readlines()
if len(parse) < 1:
    sys.exit("Error: Fail to open file")

#For the 2 dimensional matrix parsing in gaussian log file.
def parsematrx(iptname, Nbasis,matrxidx,search):
    rowdone = 0
    colndone = 0
    #total rows number:
    TotalRow = (Nbasis//5)*(Nbasis + 1) + Nbasis%5 +1
    ff = open(abs_data+objname+iptname+'.dat','w')
    count = 0
    for ii in range(TotalRow):#parse all the data area (including data and index row and column
        matrxidx+=1 # move to current line
        eachline = search[matrxidx].split()
        if 'D' in eachline[-1]:
            colndone+=1  #found a new line can be parsed!
            for xx in range(1,len(eachline)): #we parse line by line(5 elements)
                output = str(eachline[0])+' '+str(rowdone+xx)+' ' +eachline[xx].replace('D','E')+'\n'
                ff.write(output)
                count+=1
            if colndone%Nbasis == 0: #means we move to next 5-element column area
                rowdone+=5
    if count!= Nbasis*(Nbasis+1)/2:
        sys.exit("Error:Fail to form correct"+iptname+"matrix")
    ff.close()

def ERIparse(iptname,idx,search):
    align1 = search[idx+6].split()[0]
    if (align1 != "IntCnt="):
        sys.exit("Error:Fail to align ERI data_1")
    Parsing = True
    data_idx = idx + 7
    ff = open(abs_data+objname+iptname+'.dat','w')
    while (Parsing):
        eachline = search[data_idx].split()
        if (eachline[0] == "I="):
            output = str(eachline[1])+" "+str(eachline[3])+" "+str(eachline[5])+" "+str(eachline[7])+" "+str(eachline[9].replace('D','E'))+'\n'
            ff.write(output)
            data_idx+=1
        else:
            Parsing = False
    if (eachline[0] == "Leave"):
        print("Done parsing ERI")
    ff.close()

def writeenuc_Nbasis(iptname,Nb,enu):
    ff = open(abs_data+objname+iptname+".dat","w")
    ff.write(str(enu)+'\n')
    ff.write(str(Nb)+'\n')
    ff.close()

#XXX:Here, I assumed Nbasis will apear before ERI core hamiltonian data 
# So I just need to parse the array once
#Otherwise you need to parse Nbasis first 
for idx in range(len(parse)):
    eachline =  parse[idx].split()
    if all(x in eachline for x in ["nuclear","repulsion","energy"]):
        enuc = eachline[3] 
    if all(x in eachline for x in ["NBasis", "="]):
        Nbasis = int(eachline[2])
        writeenuc_Nbasis("enuc_Nb",Nbasis,enuc)
    if all(x in eachline for x in ["***","Overlap","***"]):
        #each row only have 5 items each column will have Nbasis
        parsematrx("overlap",Nbasis,idx,parse) 
    if all(x in eachline for x in ["Core","Hamiltonian"]):
        parsematrx("coreH",Nbasis,idx,parse)
    if all(x in eachline for x in ["Dumping","Two-Electron","integrals"]):
        ERIparse("eri",idx,parse)

