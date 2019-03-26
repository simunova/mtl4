#!/usr/bin/env python
import os, os.path, shutil, sys, re, fileinput

if len(sys.argv) != 2:
  print "syntax:", sys.argv[0], "tutorialfile"
  sys.exit(1)



f= open(sys.argv[1])
lines= f.readlines()
f.close()
pages, navis= [], []
for i in xrange(len(lines)):
  if r"\if Navigation \endif" in lines[i]:
    navis.append(i)
    for j in xrange(i-1, -1, -1):
      if r"\page" in lines [j]:
        tmatch= re.search(r'\\page *(\S*)', lines[j])
        if not tmatch:
          print "No preceeding page entry found"
          sys.exit(1)
        #print tmatch.group(1)
        pages.append(tmatch.group(1))
        break
#print "Pages ", pages
#print "Navis ", navis

refname= os.path.splitext(os.path.basename(sys.argv[1]))[0]
for n in xrange(len(navis)):
  s= r" "
  if n > 0: s= s + r" Return to \ref " + pages[n-1] + " "
  s= s + 30 * "&nbsp;" + r' \ref ' + refname + r' "Table of Content" ' + 30 * "&nbsp;"
  if n < len(navis)-1: s= s + r" Proceed to \ref " + pages[n+1]
  lines[navis[n] + 1]= s + r" " + "\n"

f= open(sys.argv[1], "w")
f.writelines(lines)
f.close()
