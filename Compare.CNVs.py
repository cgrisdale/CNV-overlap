#!/usr/bin/python
#Cameron Grisdale
#Nov 30, 2017

import sys
import re
#from collections import Counter

#Compare output of CNV software control-freec. Look for recurrent CNVs etc.

def Read_in_files(File):
  '''Read in CNV file'''

  cnvl,cnvd,chrm,start,end,copyn,typea,tml=[],{},'','','','','',[]

  with open(File, 'r') as f:

    fname=str(File)

    for line in f:
      line=line.strip()
      chrm,start,end,copyn,typea=line.split('\t')

      #Define ranges #ranges end at -1 ie. (5,8) is 5,6,7
      cnvrange=xrange(int(start),int(end))

      tml=[chrm,cnvrange,typea,fname]

      #Load dictionary using chrm number/letter as key
      cnvd=Load_dict(cnvd,tml)

      #Add to query list
      cnvl.append(tml)

    #print len(cnvd[fname])

  return cnvd,cnvl


def Load_dict(d,l):
  '''Load dictionary, adding lists to value for each chrm key'''

  if l[0] in d:
    d[l[0]].append(l)

  else:
    d[l[0]]=[l]

  return d


def Compare_ranges(cnvlist,qlist):
  '''Compare genomic coordinate ranges looking for overlap of CNV regions and frequency of recurrence'''
  cnvl,tmp=[],[]

  #Go through all CNV lines once (query list)
  for x in qlist:

    chrm,cnvr,etype,filen=x[0:]

    #Go through list of CNV dictionaries (key=chrm)
    for y in cnvlist:

      #Search against all items in dictionary of relevant chromosome (ie y[chrm])
      for i in y[chrm]:

        if range(max(cnvr[0],i[1][0]),min(cnvr[-1],i[1][-1])+1): #if they overlap

          if filen!=i[3] and etype==i[2]:#if not from the same file and are both gain or loss

            overlap=range(max(cnvr[0],i[1][0]),min(cnvr[-1],i[1][-1])+1)

            tmp=[chrm,xrange(min(cnvr[0],i[1][0]),max(cnvr[-1],i[1][-1])+1),etype,filen] #take min and max to make largest range

            print "Found match and not from same file and both gain or loss",i,x

            #Need to check if it's in final list before adding it

            if cnvl:

              for z in range(0,len(cnvl)):

                if range(max(cnvl[z][1][0],tmp[1][0]),min(cnvl[z][1][-1],tmp[1][-1]+1)): #overlap with range in final list, add to it
                  cnvl[z].append(tmp[-1]) #add file name for counting later

                else:

                  cnvl.append(tmp)

            else: #for first round when cnvl is still empty

              cnvl.append(tmp)

        else:

          pass
#########NOT MAKING OUTPUT########
  return cnvl


if __name__ == "__main__":

  if len(sys.argv)>2:
    Files=sys.argv[1:] #

  else:
    sys.exit("Script works for multiple files only; exiting")

  CV,QL=[],[]

  for i in range(len(Files)):
    cv,cl=Read_in_files(Files[i])
    QL.append(cl)
    CV.append(cv)
  print CV

  #Flatten list of lists of lists down one level to list of lists
  QL=[item for sublist in QL for item in sublist]
  print '\n'
  print QL

  c=Compare_ranges(CV,QL)
  print '\n'
  print c
#  outf=open('CV.outfile.tsv', 'w
