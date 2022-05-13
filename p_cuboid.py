##########################################################################################
##                                                                                      ##
##   Cuboid generation script for Blender 2.83    Copyright (C) 2020  Pan Thistle       ##
##                                                                                      ##
##   This program is free software: you can redistribute it and/or modify               ##
##   it under the terms of the GNU General Public License as published by               ##
##   the Free Software Foundation, either version 3 of the License, or                  ##
##   (at your option) any later version.                                                ##
##                                                                                      ##
##   This program is distributed in the hope that it will be useful,                    ##
##   but WITHOUT ANY WARRANTY; without even the implied warranty of                     ##
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                      ##
##   GNU General Public License for more details.                                       ##
##                                                                                      ##
##   You should have received a copy of the GNU General Public License                  ##
##   along with this program.  If not, see <https://www.gnu.org/licenses/>.             ##
##                                                                                      ##
##########################################################################################


import bpy
import bmesh
from mathutils import Vector


# ------------------------------------------------------------------------
#    GENERAL PURPOSE UTILITY FUNCTIONS
# ------------------------------------------------------------------------

def get_seam_ids(apts, bpts, nvs, sid = 0):
    ''' Return list of boundary-vertex indices for subdivided plane
        INPUT:
            - apts = integer, side-a points 
            - bpts = integer, side-b points
            - nvs = integer, total points
            - sid = integer, start index
        OUTPUT:
            - ids = list of boundary-vertex indices
        NOTES: 
            - general purpose helper function (quad-based, rectangular planes only)
            - does NOT check the validity of input data
    '''
    aseg = apts - 1
    ids = [i for i in range(sid, sid + apts)]
    for i in range(1, bpts):
        loop = sid + aseg + i * apts
        ids += [loop]
    ids += [i for i in range(loop - aseg, nvs - 1)][::-1]
    red = []
    for i in range(1, bpts - 1):
        red += [sid + i * apts]
    ids += red[::-1]
    return ids

def get_bridge_faces(aid, bid, pts, flip = False, close = True):
    ''' Return faces for bridge loop
        INPUT:
            - aid = list of integers, indices of plane-a boundary vertices 
            - bid = list of integers, indices of plane-b boundary vertices 
            - pts = integer, bridge points (pts = len(aid) = len(bid))
            - flip = boolean, flip order of verts (face normal direction)
            - close = boolean, include end face
        OUTPUT:
            - faces = list of lists of face vertex indices
        NOTES: 
            - general purpose helper function
            - input vertex indices must be in the same order
            - does NOT check the validity of input data
    '''
    faces = []
    for i in range(pts - 1):
        face = [aid[i], bid[i], bid[i + 1], aid[i + 1]]
        face = face[::-1] if flip else face
        faces.append(face)
    if close:
        face = [aid[-1], bid[-1], bid[0], aid[0]]
        face = face[::-1] if flip else face
        faces.append(face)
    return faces

def get_ring_faces(iter, rpts, sid = 0, flip = False, close = True):
    ''' Return faces for ring loop
        INPUT:
            - iter = integer, number of ring iterations 
            - rpts = integer, ring points
            - sid = integer, start index
            - flip = boolean, flip order of verts (face normal direction)
            - close = boolean, include end face
        OUTPUT:
            - faces = list of lists of face vertex indices
        NOTES: 
            - general purpose helper function
            - does NOT check the validity of input data
    '''
    faces = []    
    for i in range(iter):
        j_loop = rpts * i
        for j in range(sid, sid + rpts - 1):
            face = [j + j_loop, j + j_loop + 1, j + j_loop + rpts + 1, j + j_loop + rpts]
            face = face[::-1] if flip else face
            faces.append(face)
        if close:
            face = [sid + j_loop, sid + j_loop + rpts, sid + j_loop + 2 * rpts - 1, sid + j_loop + rpts - 1]
            face = face[::-1] if flip else face
            faces.append(face)
    return faces

def get_plane_faces(iter, lpts, sid = 0, flip = False):
    ''' Return faces for subdivided plane
        INPUT:
            - iter = integer, number of line iterations 
            - lpts = integer, line points
            - sid = integer, start index
            - flip = boolean, flip order of verts (face normal direction)
        OUTPUT:
            - faces = list of lists of face vertex indices
        NOTES: 
            - general purpose helper function (quad-based, rectangular planes only)
            - does NOT check the validity of input data
    '''
    faces = []
    for i in range(iter):
        j_loop = lpts * i
        for j in range(sid, sid + lpts - 1):
            face = [j + j_loop, j + j_loop + 1, j + j_loop + lpts + 1, j + j_loop + lpts]
            face = face[::-1] if flip else face
            faces.append(face)
    return faces


# ------------------------------------------------------------------------
#    CUBOID DATA [MAIN FUNCTION]
# ------------------------------------------------------------------------

def cuboid_mesh_data(rad, seg):
    ''' Return verts, faces for cuboid of custom dimensions and subdivisions
        INPUT:
            - rad = list of 3 floats (x-y-z radius)
            - seg = list of 3 integers (x-y-z segments) 
        OUTPUT:
            - verts = list of all vertex coordinates (mathutils Vector) 
            - faces = list of lists of face vertex indices
        NOTES: 
            - start/end faces are built on Y-Z axes, middle face-rings run from left to right on X axis
            - uses 'get_seam_ids', 'get_bridge_faces', 'get_plane_faces', 'get_ring_faces' functions
            - does NOT check the validity of input data
    '''
    # number of points (per axis)
    pts = [seg[i] + 1 for i in range(3)]
    # distance between points (per axis) 
    fac = [2 * rad[i] / seg[i] for i in range(3)]
    # left side in Y-Z (verts and faces)
    verts = [Vector((-rad[0], -rad[1] + j * fac[1], -rad[2] + i * fac[2])) for i in range(pts[2]) for j in range(pts[1])]
    faces = get_plane_faces(seg[2], pts[1], 0, True)
    nvs = len(verts)
    # seam indices
    sid = get_seam_ids(pts[1], pts[2], nvs)
    # ring points
    rpts = len(sid)
    # store total side points and ring-excluded points (will use later)
    yz_pts = nvs
    cpts = nvs - rpts
    if seg[0] > 1:
        # path along X, ring in Y-Z
        path = [Vector((i * fac[0], 0, 0)) for i in range(pts[0])]
        ring = [verts[sid[j]] for j in range(rpts)]
        rid = [nvs + i for i in range(rpts)]
        # bridge left side with first ring loop (faces only)
        faces += get_bridge_faces(rid, sid, rpts)
        # middle ring (verts and faces)
        verts += [ring[j] + path[i] for i in range(1, seg[0]) for j in range(rpts)]
        faces += get_ring_faces(seg[0]-2, rpts, nvs)
        nvs = len(verts)
    # right side in Y-Z (verts and faces)
    verts += [Vector((rad[0], -rad[1] + j * fac[1], -rad[2] + i * fac[2])) for i in range(pts[2]) for j in range(pts[1])]
    faces += get_plane_faces(seg[2], pts[1], nvs)
    nvs = len(verts)
    # adjust start index for last seam 
    tmp_id = nvs - yz_pts
    # seam indices
    eid = get_seam_ids(pts[1], pts[2], nvs, tmp_id)    
    if seg[0] > 1:
        # adjust start index for last ring loop
        tmp_id = nvs + cpts - 2 * yz_pts
        rid = [tmp_id + i for i in range(rpts)]
        # bridge last ring loop with right side (faces only)
        faces += get_bridge_faces(eid, rid, rpts)
    else:
        # bridge left side with right side (faces only)
        faces += get_bridge_faces(eid, sid, rpts)
    return verts, faces


# ------------------------------------------------------------------------
#    MESH GENERATION FUNCTIONS
# ------------------------------------------------------------------------

def cuboid_mesh_bnops(name, rad, seg):
    ''' Using bmesh module - safer/slower
    '''
    verts, faces = cuboid_mesh_data(rad, seg)
    me = bpy.data.meshes.new(name)
    bm = bmesh.new(use_operators = False)
    bmv = [bm.verts.new(v) for v in verts]
    for f in faces:        
        bm.faces.new((bmv[f[0]], bmv[f[1]], bmv[f[2]], bmv[f[3]]))
    bm.to_mesh(me)
    me.update()
    bm.free()
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)    
    
def cuboid_mesh_pydat(name, rad, seg):
    ''' Using from_pydata - risky/faster
    '''
    verts, faces = cuboid_mesh_data(rad, seg)
    me = bpy.data.meshes.new(name)
    me.from_pydata(verts, [], faces)
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)


# ------------------------------------------------------------------------
#    TESTING
# ------------------------------------------------------------------------

seg = [10, 7, 4]           # sequence of 3 integers (resolution)
rad = [1.2, 1.0, 0.76]     # sequence of 3 floats   (dimensions)
name = 'cuboid_tmp'        # string

# 1. Play it safe: using bmesh module is generally safer in that it 
#    will not freeze/crash Blender if your script has errors. I tend
#    to use it during the development process

cuboid_mesh_bnops(name, rad, seg)

# 2. Go for speed: using from_pydata is less forgiving to script errors 
#    but it gets much faster as you up the segment values. I would use 
#    it only after I have thoroughly tested a script

#cuboid_mesh_pydat(name, rad, seg)
