# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

import sympy as sm


def check_simplepolygon(VL):
    '''
    Parameters
    ----------
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    True/False depending on whether VL is a simple polygon or not.

    '''

    edge = []
    for i in range(1,VL.shape[0]):
        edge.append(sm.Segment(sm.Point(VL[i-1][0],VL[i-1][1]),sm.Point(VL[i][0],VL[i][1])))
    edge.append(sm.Segment(sm.Point(VL[VL.shape[0]-1][0],VL[VL.shape[0]-1][1]),sm.Point(VL[0][0],VL[0][1])))
    
    if VL.shape[0]<4:
        if VL.shape[0]==3:
            return not sm.Point.is_collinear((VL[0][0],VL[0][1]),(VL[1][0],VL[1][1]),(VL[2][0],VL[2][1]))
        return False
    
    check=True
    k=1
    for i in range(len(edge)-2):
        if i!=0:
            k=0
        for j in range(i+2,len(edge)-k):
            if edge[i].intersection(edge[j]) !=[]:
                check=False
                break
        if check==False:
            break
            
    return check


def check_convexity(VL):
    '''
    Parameters
    ----------
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    True/False depending on whether VL forms a boundary of a convexy polygon.

    '''
    points = []
    for i in range(VL.shape[0]):
        points.append(sm.Point(VL[i][0],VL[i][1]))
        
    p = sm.Polygon(*points)
    vertex = p.vertices
    orientation = p._isright(vertex[-2:][0],vertex[-2:][1],vertex[0])
    check = True
    for i in range(1,len(vertex)):
        if orientation != p._isright(vertex[i-2],vertex[i-1],vertex[i]):
            check = False
            break
    if check_simplepolygon(VL):
        return check
    return False

def point_membership(P,VL):
    '''
    Parameters
    ----------
    P : a 2D point example, P = np.array([1,2])
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    Should an integer type 
    1 if the P is inside the boundaries defined by VL
    0 if the P is outside the boundaries defined by VL
    -1 if the P is on the boundary defined by VL

    '''
    points = []
    for i in range(VL.shape[0]):
        points.append(sm.Point(VL[i][0],VL[i][1]))
    
    p = sm.Polygon(*points)
    
    edges = p.sides
    
    for i in range(len(edges)):
        if P in edges[i]:
            return -1
    cnt=0
    k=1
    while k<100:
        repeat = False
        r = sm.Ray((P[0],P[1]),angle = (((sm.pi/8)*(k+2))))
        intersect = r.intersection(p)
        #print(intersect)
        if any(isinstance(a,sm.Segment) for a in intersect):
            continue
        else:
            for i in intersect:
                if i in points :
                    repeat = True
                    break
            if repeat==False:
                cnt = len(intersect)
                break
        k=k+1
        
    if cnt%2!=0 :
        return 1
    
    else :
        return 0

def find_intersection(VL1,VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the intersection 
     of the two solids 1 and 2.

    '''
    newpoints1 = findIntersection(VL1,VL2)
    newpoints2 = findIntersection(VL2, VL1)
    
    
    vec1 = compute_ds(newpoints1, VL2 , 0)
    vec2 = compute_ds(newpoints2, VL1 , 1)
  

    inside = [x for x in vec1 if x[2]==0]
    
    for i in vec2:
        if i[2]==0:
            inside.append(i)
    
    if len(inside)==0:
        return np.array([])
    
    VL_int = joinAllPoints(inside)
    
    return VL_int


def find_union(VL1,VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the union 
     of the two solids 1 and 2.

    '''
    newpoints1 = findIntersection(VL1,VL2)
    newpoints2 = findIntersection(VL2, VL1)
    
    
    
    vec1 = compute_ds(newpoints1, VL2, 0)
    vec2 = compute_ds(newpoints2, VL1, 1)
    
    
    inside = [x for x in vec1 if x[2]==1]
    
    for i in vec2:
        if i[2]==1:
            inside.append(i)
    
    #return inside
    if len(inside)==0:
        return np.array([])
    
    VL_int = joinAllPoints(inside)
    
    return VL_int

def find_difference(VL1, VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the difference  
     of the two solids 1 and 2.
     S1-S2.
    '''
    newpoints1 = findIntersection(VL1,VL2)
    newpoints2 = findIntersection(VL2, VL1)
    
    vec1 = compute_ds(newpoints1, VL2, 0)
    vec2 = compute_ds(newpoints2, VL1, 1)
    
    inside = [x for x in vec1 if x[2]==1]
    
    for i in vec2:
        if i[2]==0:
            inside.append(i)
    
    #return inside
    if len(inside)==0:
        return np.array([])
    
    VL_int = joinAllPoints(inside)
    
    return VL_int


def findIntersection(VL1,VL2):
    points1 = []
    points2 = []
    for i in range(VL1.shape[0]):
        points1.append((VL1[i][0],VL1[i][1]))
    
    for i in range(VL2.shape[0]):
        points2.append((VL2[i][0],VL2[i][1]))
    
    polygon1 = sm.Polygon(*points1)
    polygon2 = sm.Polygon(*points2)
    newpoints1 = []
    for i in polygon1.sides:
        intersec = [i.intersection(j) for j in polygon2.sides]
        intersec = [intersectionPoint for subPoint in intersec for intersectionPoint in subPoint]
        intersec = [currPoint for currPoint in intersec if currPoint != []]
        intersec =sorted(intersec,key=lambda e: i.p1.distance(e))
        if isinstance(i,sm.Segment):
            newpoints1.append(i.p1)
            #print(newpoints1)
        for k in intersec:
            if k!=[] and not isinstance(k,sm.Segment) and  not k in newpoints1:
                newpoints1.append(k)
    return newpoints1

# 0 is edge inside of polygon
# 1 is edge outside the polygon


def compute_ds(newpoints1,VL2,inout):
    vec1 = []
    for i in range(1,len(newpoints1)):
        a = newpoints1[i-1]
        b = newpoints1[i]
        aLocation = point_membership(np.array([(a.x + b.x)/2,(a.y + b.y)/2]), VL2)
        
        if aLocation == 1:
            vec1.append([a,b,0])
        elif aLocation == 0:
            vec1.append([a,b,1])
        elif aLocation == -1:
            if inout == 0:
                vec1.append([a,b,0])
            elif inout == 1:
                vec1.append([a,b,1])
    
    aLocation = point_membership(np.array([(newpoints1[len(newpoints1)-1].x + newpoints1[0].x) /2 , (newpoints1[len(newpoints1)-1].y + newpoints1[0].y)/2]) ,VL2)
    
    if aLocation == 1 :
        vec1.append([newpoints1[len(newpoints1)-1],newpoints1[0],0])
    elif aLocation == 0 :
        vec1.append([newpoints1[len(newpoints1)-1],newpoints1[0],1])
    #elif aLocation == -1 :
        #vec1.append([newpoints1[len(newpoints1)-1],newpoints1[0],1])
        #vec1.append([newpoints1[len(newpoints1)-1],newpoints1[0],0])
    
    return vec1


def joinAllPoints(inside):
    combined = [inside[0][0]]
    current = inside[0][1]
    del inside[0]
    
    while len(inside)!=0:
        for i in range(len(inside)):
            if inside[i][0] == current:
                combined.append(current)
                current = inside[i][1]
                del inside[i]
                break
        for i in range(len(inside)):
            if inside[i][1] == current:
                combined.append(current)
                current = inside[i][0]
                del inside[i]
                break
        for i in combined:
            if i==current:
                if len(inside) !=0:
                    #combined.append(current)
                    combined.append(inside[0][0])
                    current = inside[0][1]
                    del inside[0]
    
    if combined[0] == combined[len(combined)-1]:
        combined.pop()
        
    VL_int = np.zeros((len(combined),2))
    
    for i in range(len(combined)):
        VL_int[i][0] = combined[i].x
        VL_int[i][1] = combined[i].y
    
    return VL_int