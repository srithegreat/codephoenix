import sys

def getpsms(fname):
    mytemp=open("PSMs_SetNB.tsv","w")
    counter=0
    noise=500
    
    for each in open(fname):
        
        header=each.strip().replace("\"","").split("\t")

        if counter==0:
            sequenceindex=header.index("Annotated Sequence")
            rank=header.index("Rank")
            chan1=header.index("126")
            chan2=header.index("127_N")
            chan3=header.index("127_C")
            chan4=header.index("128_N")
            chan5=header.index("128_C")
            chan6=header.index("129_N")
            chan7=header.index("129_C")
            chan8=header.index("130_N")
            chan9=header.index("130_C")
            chan10=header.index("131")
            modification=header.index("Modifications")
            area=header.index("Area")
            spectrumfile=header.index("Spectrum File")
            scanid=header.index("First Scan")
            counter+=1
        else:
            if header[rank]=="1":
                if header[chan1]=="":
                    header[chan1]="0"
                if header[chan2]=="":
                    header[chan2]="0"
                if header[chan3]=="" :
                    header[chan3]="0"
                if header[chan4]=="" :
                    header[chan4]="0"
                if header[chan5]=="" :
                    header[chan5]="0"
                if header[chan6]=="" :
                    header[chan6]="0"
                if header[chan7]=="" :
                    header[chan7]="0"
                if header[chan8]=="" :
                    header[chan8]="0"
                if header[chan9]=="" :
                    header[chan9]="0"
                if header[chan10]=="" :
                    header[chan10]="0"  
                """
                #SetA                                                              
                new_ch1=float(header[chan1])-0.002*float(header[chan3])
                new_ch2=float(header[chan2])-0.009*float(header[chan4])
                new_ch3=float(header[chan3])-0.469*float(header[chan1])-float(header[chan5])*0.0053
                new_ch4=float(header[chan4])-0.054*float(header[chan2])-float(header[chan6])*0.007
                new_ch5=float(header[chan5])-0.046*float(header[chan3])-float(header[chan7])*0.013
                new_ch6=float(header[chan6])-0.047*float(header[chan4])-float(header[chan8])*0.012
                new_ch7=float(header[chan7])-0.0259*float(header[chan5])-float(header[chan9])*0.029
                new_ch8=float(header[chan8])-0.03*float(header[chan6])-float(header[chan10])*0.0236
                new_ch9=float(header[chan9])-0.025*float(header[chan7])
                new_ch10=float(header[chan10])-0.028*float(header[chan8])
                
                #setB_C
                                                                              
                new_ch1=float(header[chan1])-0.002*float(header[chan3])
                new_ch2=float(header[chan2])-0.009*float(header[chan4])
                new_ch3=float(header[chan3])-0.05*float(header[chan1])-float(header[chan5])*0.005
                new_ch4=float(header[chan4])-0.05*float(header[chan2])-float(header[chan6])*0.007
                new_ch5=float(header[chan5])-0.046*float(header[chan3])-float(header[chan7])*0.013
                new_ch6=float(header[chan6])-0.047*float(header[chan4])-float(header[chan8])*0.012
                new_ch7=float(header[chan7])-0.032*float(header[chan5])-float(header[chan9])*0.015
                new_ch8=float(header[chan8])-0.033*float(header[chan6])-float(header[chan10])*0.021
                new_ch9=float(header[chan9])-0.025*float(header[chan7])
                new_ch10=float(header[chan10])-0.028*float(header[chan8])
                #setD
                                                                         
                new_ch1=float(header[chan1])-0.002*float(header[chan3])
                new_ch2=float(header[chan2])-0.009*float(header[chan4])
                new_ch3=float(header[chan3])-0.05*float(header[chan1])-float(header[chan5])*0.005
                new_ch4=float(header[chan4])-0.05*float(header[chan2])-float(header[chan6])*0.007
                new_ch5=float(header[chan5])-0.046*float(header[chan3])-float(header[chan7])*0.013
                new_ch6=float(header[chan6])-0.047*float(header[chan4])-float(header[chan8])*0.012
                new_ch7=float(header[chan7])-0.032*float(header[chan5])-float(header[chan9])*0.015
                new_ch8=float(header[chan8])-0.033*float(header[chan6])-float(header[chan10])*0.021
                new_ch9=float(header[chan9])-0.025*float(header[chan7])
                new_ch10=float(header[chan10])-0.028*float(header[chan8])
                #setNA_NB
                
                """                                                               
                new_ch1=float(header[chan1])-0.002*float(header[chan3])
                new_ch2=float(header[chan2])-0.009*float(header[chan4])
                new_ch3=float(header[chan3])-0.0469*float(header[chan1])-float(header[chan5])*0.0053
                new_ch4=float(header[chan4])-0.054*float(header[chan2])-float(header[chan6])*0.007
                new_ch5=float(header[chan5])-0.046*float(header[chan3])-float(header[chan7])*0.013
                new_ch6=float(header[chan6])-0.047*float(header[chan4])-float(header[chan8])*0.012
                new_ch7=float(header[chan7])-0.0259*float(header[chan5])-float(header[chan9])*0.029
                new_ch8=float(header[chan8])-0.03*float(header[chan6])-float(header[chan10])*0.0236
                new_ch9=float(header[chan9])-0.025*float(header[chan7])
                new_ch10=float(header[chan10])-0.028*float(header[chan8])
                
                
                if new_ch1=="" or float(new_ch1)-noise<0:
                    new_ch1="0"
                else: 
                    new_ch1=str(new_ch1-noise)
                if new_ch2=="" or new_ch2-noise<0:
                    new_ch2="0"
                else: 
                    new_ch2=str(new_ch2-noise)
                if new_ch3=="" or float(new_ch3)-noise<0:
                    new_ch3="0"
                else: 
                    new_ch3=str(float(new_ch3)-noise)
                if new_ch4=="" or float(new_ch4)-noise<0:
                    new_ch4="0"
                else: 
                    new_ch4=str(float(new_ch4)-noise  )  
                    
                if new_ch5=="" or float(new_ch5)-noise<0:
                    new_ch5="0"
                else: 
                    new_ch5=str(float(new_ch5)-noise)
                if new_ch6=="" or float(new_ch6)-noise<0:
                    new_ch6="0"
                else: 
                    new_ch6=str(float(new_ch6)-noise)
                if new_ch7=="" or float(new_ch7)-noise<0:
                    new_ch7="0"
                else: 
                    new_ch7=str(float(new_ch7)-noise)
                if new_ch8=="" or float(new_ch8)-noise<0:
                    new_ch8="0"
                else: 
                    new_ch8=str(float(new_ch8)-noise)
                if new_ch9=="" or float(new_ch9)-noise<0:
                    new_ch9="0"
                else: 
                    new_ch9=str(float(new_ch9)-noise)
                if new_ch10=="" or float(new_ch10)-noise<0:
                    new_ch10="0"                                                            
                else: 
                    new_ch10=str(float(new_ch10)-noise)
                    
                mytemp.write( "\t".join([header[sequenceindex].upper(),str(new_ch1),str(new_ch2),str(new_ch3),str(new_ch4),str(new_ch5),str(new_ch6),str(new_ch7),str(new_ch8),str(new_ch9),str(new_ch10),header[modification],header[area],header[scanid],header[spectrumfile]])+"\n")
    mytemp.close()
    #group_psms("psms_temp.txt")
    


    

def group_psms(f):
    peptide_gene={}
    for each in open("/home/srikanth/Downloads/CellMap/PSMs/Final_MS3_SEQUESTHT/Only_geneSpecific_peptides_all.txt"):
        if len(each.strip().split("\t"))>1:
            peptide_gene[each.strip().split("\t")[0]]=each.strip().split("\t")[1] 
            
    
    
    
    pepdict={}
    for each in open(f):
        if each.strip().split("\t")[0] in pepdict:
            pepdict[each.strip().split("\t")[0]].append([float(f) for f in each.strip().split("\t")[1:]])
        else:
            pepdict[each.strip().split("\t")[0]]=[]
            pepdict[each.strip().split("\t")[0]].append([float(f) for f in each.strip().split("\t")[1:]])
            
            
    for peptide,values in pepdict.iteritems():
        summed=[sum(l) for l in zip(*values)]
        if peptide in peptide_gene:
            print peptide_gene[peptide]+"\t"+peptide+"\t"+"\t".join([str(h) for h in summed])

getpsms(sys.argv[1])
#group_psms(sys.argv[1])
    