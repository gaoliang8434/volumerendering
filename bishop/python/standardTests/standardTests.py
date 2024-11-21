#!/usr/bin/python

import os
import sys


commonOptions = "-nbthreads 10 "

testList = [
           #[ testScriptName.py, testScriptOptions, testOutputName ]
	   [ "csg.py",              commonOptions + "-csgoperation union -case union",             "csg.exr"  ],
	   [ "csg.py",              commonOptions + "-csgoperation blend -case blend",             "csg.exr"  ],
	   [ "csg.py",              commonOptions + "-csgoperation cutout -case cutout",           "csg.exr"  ],
	   [ "csg.py",              commonOptions + "-csgoperation othercutout -case othercutout", "csg.exr"  ],
	   [ "csg.py",              commonOptions + "-csgoperation intersect -case intersect",     "csg.exr"  ],
	   [ "renderpyro.py",       commonOptions + "-case implicitOctaves2 -octaves 2",           "pyro.exr"  ],
	   [ "renderpyro.py",       commonOptions + "-case implicitOctaves6 -octaves 6",           "pyro.exr"  ],
	   [ "renderjack.py",       commonOptions + "-case implicit",                              "jack.exr"  ],
	   [ "renderbunny.py",      commonOptions + "-case levelset",                              "bunny.exr"  ],
	   [ "sphere.py",           commonOptions + "-case simple",                                "sphere.exr"  ],
	   [ "box.py",              commonOptions + "-case rounded",                               "box.exr"  ],
	   [ "travelingpyro.py",    commonOptions + "-case implicit -octaves 6",                   "travelingpyro.exr"  ]
        ]


def makeTest( t ):
	testcmd = "./" + t[0] + " " + t[1] + " -name " + t[2] 
	return testcmd




if len(sys.argv) <= 1:
	for test in testList:
		testcmd = makeTest( test )
		print testcmd
		os.system(testcmd)

if len(sys.argv) > 1:
	for case in sys.argv[1:]:
		test = testList[int(case)]
		testcmd = makeTest( test )
		print testcmd
		os.system(testcmd)

