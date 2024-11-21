#!/usr/bin/python


import sys
import os



def CmdLineFindIndex( tag ):
	for i in range(len(sys.argv)):
		if sys.argv[i] == tag:
			return i
	return -1

def CmdLineFind( tag, defaultvalue ):
	i = CmdLineFindIndex(tag)
	if i > 0:
		if i < len(sys.argv)-1:
			return sys.argv[i+1]
	return defaultvalue


doit = CmdLineFindIndex("-h")
if doit>0:
	print "makepbsscripts:"
	print "-script        Name of script to run                                          [volumerender.py]"
	print "-nbnodes       Number of nodes to ask pbs for                                 [1]"
	print "-job           Name of pbs job                                                [MyJob]"
	print "-cores         Number of cores per node to ask pbs for                        [8]"
	print "-mem           Memory in GB to ask pbs for                                    [16]"
	print "-base          Base directory for all files                                   [scratch]"
	print "-first         First pbs file to generate                                     [1]"
	print "-last          Last pbs file to generate                                      [first]"
	print "-walltime      Walltime in hours to ask for                                   [20]"
	print "-extras        Everything after this option is passed verbatim to the script  []"
	print ""
	sys.exit(0)


scriptToRun = CmdLineFind( "-script", "nothing.py" )
nbNodes =  int(CmdLineFind( "-nbnodes", 1 ))
jobLabel = str( CmdLineFind( "-job", "MyJob" ) )
ppn = int( CmdLineFind("-cores", 8 ) )
memory = int( CmdLineFind("-mem", 16 ) )
base = str( CmdLineFind("-base", "scratch" ) )


first = int( CmdLineFind("-first", 1) )
last = int( CmdLineFind("-last", first))



walltime = int( CmdLineFind("-walltime", 20) )
wally = str(walltime).zfill(2)

user = str(os.getenv("USER"))

jobextras = ""
extras = int( CmdLineFindIndex("-extras") )
if extras > 0:
	jobextras = " -extras "
	for task in sys.argv[extras+1:]:
		jobextras = jobextras + str(task) + " "

imageDirectory = "/" + base + "/" + user + "/projects/images/" + jobLabel
if os.path.exists(imageDirectory) == False:
        os.mkdir(imageDirectory)


for frame in range(first,last+1):
	padframe = str(frame).zfill(4)
	pbsname = jobLabel + padframe
	boilerplate = "#PBS -N " + pbsname + "\n"
	boilerplate = boilerplate + "#PBS -l select=" + str(nbNodes) + ":ncpus=" + str(ppn)
	boilerplate = boilerplate + ":mem=" + str(memory) + "gb"
	boilerplate = boilerplate + ",walltime=" + wally + ":00:00\n"
	boilerplate = boilerplate + "#PBS -j oe\n"
	job = boilerplate + "\n"
	job = job + str(scriptToRun) 
	job = job + " -turntableframe " + str(frame)
	job = job + " -case " + jobLabel + " "
	job = job + " -nbthreads " + str(ppn) + " "
	job = job + " -location " + str(imageDirectory) + "/ "
	job = job + jobextras + "\n"
	pbsfilename = "/" + base + "/" + user + "/projects/submitfiles/" + jobLabel + "." + padframe + ".pbs"
	print jobLabel + ": " + padframe + " ------> " + pbsfilename
	#print "Job:\n" + job + "\n\n"
	jobfile = open( pbsfilename, 'w' )
	jobfile.write( job )
	jobfile.close()






