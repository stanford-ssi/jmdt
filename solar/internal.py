import sys
import time
import numpy as np

with open(sys.argv[1]) as f:
	t = f.readlines()

mode = sys.argv[2]

# Ayyy one liners FTW
vertices = dict(map(lambda (x,y): (x, map(float,y.split(','))), map(lambda x: x.split("$"), t[0].split(" "))))
parsed = map(lambda x: [x.split(": ")[0], map(lambda v: vertices[v], x.split(": ")[1].strip().split(','))], t[1:])

vdict = vertices

import os
os.environ["PYOPENCL_CTX"] = "0"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "0"
import numpy as np
import pyopencl as cl

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)
mf = cl.mem_flags

types = ', '.join(map(lambda (i, x): '1' if x[0] == "solar" else '0', enumerate(parsed)))
vv = sum(map(lambda x: x[1], parsed), [])

norm = lambda x: x/np.linalg.norm(x)
vertices = ', '.join(['(float4)(%ff,%ff,%ff,0)'%tuple(x) for x in vv])
normals = ', '.join(['(float4)(%ff,%ff,%ff,0)'%tuple(list(norm(np.cross(np.array(x[1])-np.array(x[0]),np.array(x[2])-np.array(x[0]))))) for _,x in parsed])

src = """

__constant float TWOPI = 6.28318530718f;

bool intriangle(float4 p, float4 p0, float4 p1, float4 p2) {
    /*float A = length(cross(p1-p0, p2-p0))/2.0f;
	float alpha = length(cross(p1-p,p2-p))/(2.0f*A);
	float beta = length(cross(p2-p,p0-p))/(2.0f*A);*/

	float4 v0 = p1-p0;
	float4 v1 = p2-p0;
	float4 v2 = p-p0;
    float d00 = dot(v0, v0);
    float d01 = dot(v0, v1);
    float d11 = dot(v1, v1);
    float d20 = dot(v2, v0);
    float d21 = dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float alpha = (d11 * d20 - d01 * d21) / denom;
    float beta = (d00 * d21 - d01 * d20) / denom;
	float gamma = 1 - alpha - beta;
    
    return alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1;
}

__kernel void solar(__global const float4* initials, __global float* output, float4 origin, float4 rayvec) {
	int N = %d;
	uchar types[] = {%s};
	float4 vertices[] = {%s};
	float4 normals[] = {%s};

	int gid = get_global_id(0);

	float d;
	float tclose = 1000000.0;
	float valclose = 0.0f;
	float t;

	float m1, m2, anglesum, costheta;
	float4 p1, p2, q;
	int i;
	for (i=0; i<N; i++) {
		d = dot(normals[i], rayvec); // cos(theta), if they're unit vectors
		if (fabs(d) > 0.0001f) {
			t = dot(vertices[4*i]-initials[gid], normals[i])/d;
			if (t >= 0) {
				if (t < tclose) {
					q = initials[gid]+t*rayvec;
					if (intriangle(q, vertices[4*i+0], vertices[4*i+1], vertices[4*i+3]) || intriangle(q, vertices[4*i+1], vertices[4*i+2], vertices[4*i+3]) ) {
						tclose = t;
						valclose = %s;
					}
				}
			}
		}
	}
	output[gid] = valclose; //(initials[gid]+t*rayvec).z;
}

""" % (len(parsed), types, vertices, normals, 'types[i]*fabs(d)' if mode == 'solar' else '1.0')
#print src

prg = cl.Program(ctx, src).build()

d = 0.7;
NN=400
output_np = np.zeros(NN*NN, dtype=np.float32)
output_g = cl.Buffer(ctx, mf.WRITE_ONLY, output_np.nbytes)

inp = sys.stdin.read().split(" ")[:-1]
inp = map(lambda x: map(float, x.split(',')), inp)
#inp = [('a','b')]
t0 = time.time()
for a,b in inp:
	#
	#vector = np.array([1,0,0])
	vector = np.array([np.sin(a)*np.cos(b),np.sin(a)*np.sin(b),np.cos(a)])
	vector = vector/np.linalg.norm(vector)
	p0 = -2*vector
	D = np.dot(vector, p0)

	if vector[0] == 0:
		v2 = np.cross(vector, np.array([1,0,0]))
	else:
		v2 = np.cross(vector, np.array([0,0,1]))

	v2 = v2/np.linalg.norm(v2)
	v3 = np.cross(vector, v2)
	v3 = v3/np.linalg.norm(v3)

	a = np.linspace(-d/2,d/2,NN,dtype=np.float32)
	b = np.linspace(-d/2,d/2,NN,dtype=np.float32)
	z = np.zeros(shape=(2,NN*NN),dtype=np.float32)
	z[0] = np.repeat(a, NN)
	z[1] = np.tile(a, NN)
	mat = np.array([[v2[0], v3[0]],[v2[1], v3[1]],[v2[2], v3[2]],[0,0]], dtype=np.float32)

	vertices_np = np.dot(mat,z).T+np.array([p0[0],p0[1],p0[2],0], dtype=np.float32)
	vertices_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=vertices_np.flatten())

	origin = np.zeros(shape=4, dtype=np.float32)
	origin[0:3] = p0
	rayvec = np.zeros(shape=4, dtype=np.float32)
	rayvec[0:3] = vector

	prg.solar(queue, [len(vertices_np)], None, vertices_g, output_g, origin, rayvec)
	cl.enqueue_copy(queue, output_np, output_g)
	print np.sum(output_np)*d*d/NN/NN
#print time.time()-t0

exit()
output_np.shape = (NN, NN)

#print np.min(output_np%(2*3.1415926))
#exit()

#exit()
import matplotlib.pyplot as plt


from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
vertices=vv
#colors = colors = "bgrcmykw"
for i, poly in enumerate(parsed):
    ax.add_collection3d(Poly3DCollection([poly[1]],alpha=0.5))
#ax.scatter([x[0] for i,x in enumerate(vertices)],[x[1] for x in vertices],[x[2] for x in vertices])
ax.scatter([x[0] for x in vv],[x[1] for x in vv],[x[2] for x in vv])
ax.scatter(vertices_np[:,0], vertices_np[:,1], vertices_np[:,2])
ax.auto_scale_xyz([-a,a], [-a,a],[-a,a])
for a, p in vdict.iteritems():
    tr = proj3d.proj_transform(p[0],p[1],p[2], ax.get_proj())
    ax.annotate(str(a), (tr[0], tr[1]))
a=0.6
plt.show()
exit()


print output_np
fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(output_np, interpolation='none')
numrows, numcols = output_np.shape
def format_coord(x, y):
    col = int(x+0.5)
    row = int(y+0.5)
    if col>=0 and col<numcols and row>=0 and row<numrows:
        z = output_np[row,col]
        return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
    else:
        return 'x=%1.4f, y=%1.4f'%(x, y)

ax.format_coord = format_coord
plt.show()
