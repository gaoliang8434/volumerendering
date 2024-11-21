#!/usr/bin/python

import os
import sys


tests = [
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0001", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0001/", "-source":"0.1", "-sourceradius":"0.2" },
           {"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0001", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0001/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0003", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0003/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0002", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0002/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0003", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0003/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0004", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0004/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0005", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0005/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0006", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0006/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0007", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0007/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0008", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0008/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0009", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0009/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0010", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0010/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0011", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0011/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0012", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0012/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0013", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0013/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0002", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0002/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0004", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0004/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0005", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0005/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0014", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0014/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0015", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0015/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0016", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0016/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0017", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0017/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0018", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0018/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0019", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0019/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0020", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0020/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0021", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett2", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0021/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0022", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0022/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0023", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett6", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0023/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0024", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0024/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0", "-st_sigma":"0.1" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0025", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0025/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0006", "-dx":"0.01171875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0006/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0026", "-dx":"0.01171875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0026/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0", "-st_sigma":"0.1" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0027", "-dx":"0.01171875", "-igthreshold":"1.0", "-solver":"garrett4", "-grad_width":"4", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0027/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./planeTest8.py", "-name":"smallPlane.vdb", "-case":"1", "-version":"0028", "-dx":"0.01171875", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"8", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallPlane/0028/", "-source":"0.1", "-sourceradius":"0.2", "-erode":"1.0" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0008", "-dx":"0.046875", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0008/", "-source":"0.1", "-sourceradius":"0.2" },
           #{"executable":"./sphereTest4.py", "-name":"smallSphere.vdb", "-case":"1", "-version":"0009", "-dx":"0.0234375", "-igthreshold":"1.0", "-solver":"garrett8", "-grad_width":"6", "-nbthreads":"4", "-substeps":"1", "-location":"./tests/smallSphere/0009/", "-source":"0.1", "-sourceradius":"0.2" },
        ]






for test in tests:
	cmd = "time " + test["executable"] + " "
	for label in test.keys():
		if label != "executable":
			cmd = cmd + label + " " + test[label] + " "
	cmd = cmd + " >> testmatrix.out "
	os.system( "mkdir -p " + test["-location"] )
	print "Running " + cmd
	os.system(cmd)

