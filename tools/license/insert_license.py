#!/usr/bin/env python
import os, os.path, shutil, sys, re, fileinput

if len(sys.argv) < 2:
  print "syntax:", sys.argv[0], "licensefile [filepattern] .."
  sys.exit(1)

# replace some patterns or wildcard for all files
if len(sys.argv) == 2:
    filepatterns= ["^.*$"]
else:
    filepatterns= []
    for i in range(2, len(sys.argv)):
        if sys.argv[i] == "cpppattern":
            filepatterns.append("^\w*\..pp$|^\w*\.cc")
        elif sys.argv[i] == "cpattern":
            filepatterns.append("^\w*\.h$|^\w*\.c$")
        elif sys.argv[i] == "scriptpattern":
            filepatterns.append("^\w*\.csh$|^\w*\.tcsh$|^\w*\.sh$|^\w*\.bash$")
        else:
            filepatterns.append(sys.argv[i])


# read files of dir without stopping on empty dirs
def myFiles():
    files= []
    try:
        files= os.listdir(".")
    except OSError:
        pass   # ignore No such file Error
    return files

# filter list of files w.r.t. list of patterns
# if filepattern isn't a list it sort a weird error message
def filterFiles(files, filepatterns):
    #print files, filepatterns
    return filter(lambda s : reduce (lambda x,y: x or y,
                                     map (lambda i : re.search(i, s), filepatterns)), files)

def dirs(files):
    return filter(lambda f : os.path.isdir(f), files)

def regfiles(files):
    return filter(lambda f : os.path.isfile(f), files)

def recFilterFiles( dir, path, filepatterns ) :
    olddir= os.getcwd()
    os.chdir(dir)
    files= myFiles()
    myAbsFiles= map (lambda f : os.path.join(path, f), filterFiles(regfiles(files), filepatterns))
    # print "in", path, "matches", myAbsFiles
    for d in dirs(files) :
        if not re.search( "^\.svn$|^CVS$", d ) :     # filter out svn and CVS dirs
            myAbsFiles.extend( recFilterFiles( d, os.path.join(path, d), filepatterns ) )
    os.chdir(olddir)
    return myAbsFiles

def allFiles(filepatterns = ["^.*$"]):
    return recFilterFiles(".", "", filepatterns)
    
allMyFiles= allFiles(filepatterns)


def insertCopyright(licName, files) :
    for fname in files :        
        fcopyname= fname + ".new"
        fcopy= open(fcopyname, "w")
        f= open(fname)
        for l in f :
            match= re.search("^(.*)\$COPYRIGHT\$", l)
            if match :
                comment= match.group(1)
                fLic= open(licName)
                # now copy license as comment
                for lLic in fLic :
                    fcopy.write(comment + lLic)
                fLic.close()
                # copy rest of file 1 to 1
                for ll in f :
                    fcopy.write(ll)
            else:
                # before match copy 1 to 1
                fcopy.write(l)
        f.close()
        fcopy.close()
        os.remove(fname)
        os.rename(fcopyname, fname)

insertCopyright(sys.argv[1], allMyFiles)
