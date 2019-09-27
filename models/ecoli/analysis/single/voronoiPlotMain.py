"""
Plot Voronoi diagram

@author: Ray Chang
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019

Algorithm Reference: A. Nocaj and U. Brandes, "Computing Voronoi Treemaps: Faster, Simpler, and Resolution-independent, Computer Graphics Forum, vol. 31, no. 31, pp. 855-864, 2012. Available: 10.1111/j.1467-8659.2012.03078.x.
"""

import os
import numpy as np
import math as math
import pandas as pd
from scipy.spatial import distance_matrix
from scipy.spatial import ConvexHull
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib
import time

"""
Classes
"""
class POLYGON(object):
    #please always reorder the corners before initiating a new polygon object!
    def __init__(self, xy):
        self.xy = xy;
        self.n_corners = len(xy);
        self.edges = [];
        for i in range(self.n_corners):
            self.edges.append(self.xy[[i-1,i],:])
        self.area = PolygonArea(self.xy);
            
class VORONOI(object):
    def __init__(self, polygons, Sites, Weights, Values, canvas_obj):
        self.polygons = polygons;
        self.n_sites = len(polygons);
        self.Sites = Sites;
        self.Weights = Weights;
        self.Values = Values;
        self.canvas_obj = canvas_obj;

class RAY(object):
    def __init__(self, origin, tangent, main_sites, adjunct_site):
        #ray: origin + a*tangent, a>=0;
        self.origin = origin;
        self.tangent = tangent;
        #ray: a_r*x + b_r*y + c_r = 0; t_y*x - t_x*y + (t_x*y_int - t_y*x_int) = 0;
        self.a_r = tangent[1];
        self.b_r = -tangent[0];
        self.c_r = (tangent[0]*origin[1] - tangent[1]*origin[0]);
        self.main_sites = main_sites;
        self.adjunct_site = adjunct_site;
"""
1st tier Core Functions
"""
def Voronoi_Main_functions(Points, Values, canvas):
    tic()
    err_thres = 1E-6;
    error = float("inf");
    i_max = 75;
    canvas_obj = POLYGON(reorderPoints(canvas));
    Voronoi, error = Voronoi_treemap(canvas_obj, Points, Values, i_max, err_thres);
    print(error);
    while (error > 0.15)|(error == 0):
        Voronoi, error_new = Voronoi_treemap_recal(Voronoi, i_max, err_thres);
        print(error_new);
        if error_new <= 0.15:
            error = error_new;
            break;
        else:
            if (error_new >= error):
                Voronoi, error_new = Voronoi_treemap(canvas_obj, Points, Values, i_max, err_thres);
                print(error_new);
        error = error_new; 
    Voronoi_treemap_plot(Voronoi);
    toc()
    return Voronoi, error

def Voronoi_Main_functions_for_layer(Points, Values, canvas_obj):
    tic()
    err_thres = 1E-6;
    error = float("inf");
    i_max = 75;
    Voronoi, error = Voronoi_treemap(canvas_obj, Points, Values, i_max, err_thres);
    print(error);
    while (error > 0.15)|(error == 0):
        print('error > 0.15:', error > 0.15)
        Voronoi, error_new = Voronoi_treemap_recal(Voronoi, i_max, err_thres);
        print(error_new);
        if error_new <= 0.15:
            error = error_new;
            break;
        else:
            if (error_new >= error):
                print('error_new >= error');
                Voronoi, error_new = Voronoi_treemap(canvas_obj, Points, Values, i_max, err_thres);
                print(error_new);
        error = error_new; 
    toc()
    return Voronoi, error

def Layered_Voronoi(df, canvas, i_max, err_thres):
    df = df.sort_index();
    index = df.index;
    n_layers = len(index[0]);
    n_Points_total = len(index);
    canvas_obj = POLYGON(reorderPoints(canvas));
    for i_level in np.arange(n_layers+1):
        if i_level == 0:
            Values_0 = list(df.sum(level=0)['values']);
            #Points_0 = list(df.index.unique(level=0)); #for python 3
            Points_0 = list(index.get_level_values(0).unique()); #for python 2.7
            error_0 = np.inf;
            Voronoi_0, error_0 = Voronoi_Main_functions_for_layer(Points_0, Values_0, canvas_obj)
        elif i_level == 1:
            n_categories_0 = len(Values_0);
            Voronoi_1_all = [[] for _ in np.arange(n_categories_0)];
            Voronoi_2_all = [[] for _ in np.arange(n_categories_0)];
            n_categories_1_all = np.zeros(n_categories_0);
            Points_1_all = [[] for _ in np.arange(n_categories_0)];
            for i in np.arange(n_categories_0):
                category_0 = Points_0[i];
                Values_1 = list(df.loc[category_0].sum(level=0)['values']);
                n_categories_1 = len(Values_1);
                n_categories_1_all[i] = n_categories_1;
                Voronoi_2_all[i] = [[] for _ in np.arange(n_categories_1)];
                #Points_1 = list(df.loc[category_0].index.unique(level=0)); #for python 3
                Points_1 = list(df.loc[(category_0)].index.get_level_values(0).unique()); #for python 2.7
                Points_1_all[i] = Points_1;
                canvas_1 = Voronoi_0.polygons[i];
                Voronoi_1_all[i], error_1 = Voronoi_Main_functions_for_layer(Points_1, Values_1, canvas_1);
        elif i_level == 2:
            for i in np.arange(n_categories_0):
                category_0 = Points_0[i];
                n_categories_1 = int(n_categories_1_all[i]);
                for j in np.arange(n_categories_1):
                    canvas_2 = Voronoi_1_all[i].polygons[j];
                    category_1 = Points_1_all[i][j];
                    #Points_2 = list(df.loc[(category_0, category_1)].index.unique(level=0));#for python 3
                    Points_2 = list(df.loc[(category_0, category_1)].index.get_level_values(0).unique());#for python 2.7
                    Values_2 = list(df.loc[(category_0, category_1)].sum(level=0)['values']);
                    Voronoi_2_all[i][j], error_2 = Voronoi_Main_functions_for_layer(Points_2, Values_2, canvas_2);
    return Voronoi_0, Voronoi_1_all, Voronoi_2_all

def Layered_Voronoi_plot(Voronoi_0, Voronoi_1_all, Voronoi_2_all):
    error = 0;
    total_value = sum(Voronoi_0.Values);
    fig = plt.figure(figsize=(5,5), dpi=200)
    ax = fig.add_axes([0,0,1,1])
    n_categories_0 = len(Voronoi_2_all);
    patches = []
    COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
    [55,126,184],
    [77,175,74],
    [255,127,0],
    [152,78,163],
    [255,255,51],
    [166,86,40],
    [247,129,191],
    [228,26,28]
    ]
    COLORS = [[colorValue/255. for colorValue in color] for color in COLORS_256]
    colors_all = [];
    for k in np.arange(n_categories_0):
        category_0 = Voronoi_0.Points[k];
        n_categories_1 = len(Voronoi_2_all[k]);
        for m in np.arange(n_categories_1):
            category_1 = Voronoi_1_all[k].Points[m];
            #plot Voronoi Treemap
            canvas_obj = Voronoi_2_all[k][m].canvas_obj; 
            Sites = Voronoi_2_all[k][m].Sites; 
            Values = Voronoi_2_all[k][m].Values; 
            Points = Voronoi_2_all[k][m].Points;
            #plot canvas
            for edge_canvas in canvas_obj.edges:
                ax.plot(edge_canvas[:,0],edge_canvas[:,1],'black', lw=1.5, solid_capstyle='round', zorder=2);
            #plot polygons
            n_polygons = len(Voronoi_2_all[k][m].polygons)
            for i in range(n_polygons):
                polygon = Voronoi_2_all[k][m].polygons[i];
                polygon_plot_obj = Polygon(polygon.xy, True);
                ax.plot(polygon.xy[:,0], polygon.xy[:,1], color='black', alpha=1,linewidth=1, solid_capstyle='round', zorder=2);
                patches.append(polygon_plot_obj);
                colors = COLORS[k] + list(np.random.rand(1,3)/5)[0]
                if max(colors) >= 1:
                    colors = colors/max(colors)
                colors_all.append(colors)
                error += np.abs( polygon.area - Voronoi_0.canvas_obj.area*Voronoi_2_all[k][m].Values[i]/total_value );
            for i in range(Voronoi_2_all[k][m].n_sites):
                plt.text(Sites[i,0], Sites[i,1], Points[i], fontsize=8, horizontalalignment='center', verticalalignment='center');
    error = error/(2*Voronoi_0.canvas_obj.area);
    print('error of the whole Voronoi diagram:', error);
    p = PatchCollection(patches, facecolors=colors_all, alpha = 1)
    ax.add_collection(p)
    #ax.set_title('Biomass components');
    plt.title("Biomass components")
    ax.set_aspect('equal');
    plt.axis('off');
    return error

def Voronoi_treemap(canvas_obj, Points, Values, i_max, err_thres):
    #initialization
    n_Points = len(Points);
    if n_Points == 1:
        Voronoi = VORONOI([canvas_obj], findCOM(canvas_obj.xy).reshape((-1,2)), canvas_obj.area, Values, canvas_obj);
        error = 0.0000001;
        Voronoi.Points = Points;
        return Voronoi, error
    else:
        Sites = RandomPointsInCanvas(n_Points, canvas_obj.xy);
        Weights = np.ones(n_Points)*canvas_obj.area/n_Points;
        error = float("inf");
        #calculate initial Voronoi plot
        Entry = 0;
        Voronoi = ComputePowerDiagram(canvas_obj, Sites, Weights, Entry, Values);
        Entry = 1;
        #adjust position & weight until error becomes acceptable
        i = 0;
        while i < i_max:
            Voronoi = AdaptPositionsWeights(Voronoi)
            Voronoi = reComputePowerDiagram(Voronoi);
            Voronoi = AdaptWeights(Voronoi);
            Voronoi = reComputePowerDiagram(Voronoi);
            Voronoi, error = ComputeError(Voronoi);
            i += 1;
            if error < err_thres:
                break
        Voronoi.Points = Points;
    return Voronoi, error

def Voronoi_treemap_recal(Voronoi, i_max, err_thres):
    Points = Voronoi.Points;
    i = 0;
    while i < i_max:
        Voronoi = AdaptPositionsWeights(Voronoi)
        Voronoi = reComputePowerDiagram(Voronoi);
        Voronoi = AdaptWeights(Voronoi);
        Voronoi = reComputePowerDiagram(Voronoi);
        Voronoi, error = ComputeError(Voronoi);
        i += 1;
        if error < err_thres:
            break
    Voronoi.Points = Points;
    return Voronoi, error

def Voronoi_treemap_plot(Voronoi):
    canvas_obj = Voronoi.canvas_obj; Sites = Voronoi.Sites; Values = Voronoi.Values; Points = Voronoi.Points;
    fig = plt.figure(figsize=(5,5), dpi=200)
    ax = fig.add_axes([0,0,1,1])
    #plot canvas
    for edge_canvas in canvas_obj.edges:
        ax.plot(edge_canvas[:,0],edge_canvas[:,1],'black', lw=2, solid_capstyle='round', zorder=2);
    #plot polygons
    patches = []
    for polygon in Voronoi.polygons:
        polygon_plot_obj = Polygon(polygon.xy, True);
        ax.plot(polygon.xy[:,0], polygon.xy[:,1], color='black', alpha=1,linewidth=1, solid_capstyle='round', zorder=2);
        patches.append(polygon_plot_obj);
    colors = np.array(Values)
    p = PatchCollection(patches, cmap=matplotlib.cm.Blues, alpha=1)
    p.set_array(colors)
    ax.add_collection(p)
    #add label of polygons
    for i in range(Voronoi.n_sites):
        plt.text(Sites[i,0], Sites[i,1], Points[i], horizontalalignment='center', verticalalignment='center');
    #show axis
    ax.set_title('Biomass components');
    ax.set_aspect('equal');
    plt.axis('off');

"""
2nd tier Core funnction
"""
def ComputePowerDiagram(canvas_obj, Sites, Weights, Entry, Values):
    n_polygon = len(Weights);
    if n_polygon < 4:
        #directly compute by intersecting the bisectors
        if n_polygon == 2:
            x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
            [[x1], [x2]] = x_col; [[y1], [y2]] = y_col; [w1, w2] = Weights;
            #formula = 2(x1-x2)x + 2(y1-y2)y = (x1^2+y1^2-w1) - (x2^2+y2^2-w2)
            bisector = np.array([2*(x1-x2), 2*(y1-y2), (x1**2+y1**2-w1)-(x2**2+y2**2-w2)]);
            #find intersect point with canvas
            edge_new = findIntersectWithCanvas(bisector, canvas_obj);
            while len(edge_new) == 0: 
                Sites = RandomPointsInCanvas(n_polygon, canvas_obj.xy);
                x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
                [[x1], [x2]] = x_col; [[y1], [y2]] = y_col; [w1, w2] = Weights;
                bisector = np.array([2*(x1-x2), 2*(y1-y2), (x1**2+y1**2-w1)-(x2**2+y2**2-w2)]);
                edge_new = findIntersectWithCanvas(bisector, canvas_obj);    
            Polygons_all = Divide2Polygons(Sites, edge_new, canvas_obj);
            Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
            Voronoi = findSitesCausingDisplacement(Voronoi);
        else:
            #n_polygon = 3
            site_label = [0,1,2];
            Sites, ray_obj_all = findBisectorPositiveRay(Sites, Weights, site_label, canvas_obj);
            Polygons_all = Divide3Polygons(Sites, ray_obj_all, canvas_obj);
            Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
            Voronoi = findSitesCausingDisplacement(Voronoi);
    else:
        x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
        #1. map sites into dual space
        Dual_Sites = np.hstack((x_col,y_col,x_col**2+y_col**2-Weights.reshape(-1,1)));
        #2. find the convex hull of the points in dual space
        Faces = ConvexHull(Dual_Sites);
        simplices = Faces.simplices;
        #3. select only the lower convex hull
        #and convert the triangle in dual space into the intersection points(iPoints) in normal space
        iPoints, simplices_prune = PruneAndConvertToIPoints(Dual_Sites, simplices, canvas_obj);
        #4-1. compute the ray objects of each iPoints
        ray_obj_all = findMultiplePositiveRay(x_col, y_col, Sites, Weights, iPoints, simplices_prune, canvas_obj);
        #4-2. compute the coordinate of each polygon
        Polygons_all = DivideNPolygons(Sites, ray_obj_all, canvas_obj, simplices_prune, iPoints);     
        Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
        Voronoi = findSitesCausingDisplacement(Voronoi);
    return Voronoi

def reComputePowerDiagram(Voronoi):
    Sites = Voronoi.Sites; Weights = Voronoi.Weights; canvas_obj = Voronoi.canvas_obj; Values = Voronoi.Values;
    n_polygon = len(Weights);
    if n_polygon < 4:
        #directly compute by intersecting the bisectors
        if n_polygon == 2:
            x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
            [[x1], [x2]] = x_col; [[y1], [y2]] = y_col; [w1, w2] = Weights;
            #formula = 2(x1-x2)x + 2(y1-y2)y = (x1^2+y1^2-w1) - (x2^2+y2^2-w2)
            bisector = np.array([2*(x1-x2), 2*(y1-y2), (x1**2+y1**2-w1)-(x2**2+y2**2-w2)]);
            #find intersect point with canvas
            edge_new = findIntersectWithCanvas(bisector, canvas_obj);
            while len(edge_new) == 0: 
                Sites = RandomPointsInCanvas(n_polygon, canvas_obj.xy);
                x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
                [[x1], [x2]] = x_col; [[y1], [y2]] = y_col; [w1, w2] = Weights;
                bisector = np.array([2*(x1-x2), 2*(y1-y2), (x1**2+y1**2-w1)-(x2**2+y2**2-w2)]);
                edge_new = findIntersectWithCanvas(bisector, canvas_obj);    
            Polygons_all = Divide2Polygons(Sites, edge_new, canvas_obj);
            Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
            Voronoi = findSitesCausingDisplacement(Voronoi);
        else:
            #n_polygon = 3
            site_label = [0,1,2];
            Sites, ray_obj_all = findBisectorPositiveRay(Sites, Weights, site_label, canvas_obj);
            Polygons_all = Divide3Polygons(Sites, ray_obj_all, canvas_obj);
            Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
            Voronoi = findSitesCausingDisplacement(Voronoi);
    else:
        x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
        #1. map sites into dual space
        Dual_Sites = np.hstack((x_col,y_col,x_col**2+y_col**2-Weights.reshape(-1,1)));
        #2. find the convex hull of the points in dual space
        Faces = ConvexHull(Dual_Sites);
        simplices = Faces.simplices;
        #3. select only the lower convex hull
        #and convert the triangle in dual space into the intersection points(iPoints) in normal space
        iPoints, simplices_prune = PruneAndConvertToIPoints(Dual_Sites, simplices, canvas_obj);
        #4-1. compute the ray objects of each iPoints
        ray_obj_all = findMultiplePositiveRay(x_col, y_col, Sites, Weights, iPoints, simplices_prune, canvas_obj);
        #4-2. compute the coordinate of each polygon
        Polygons_all = DivideNPolygons(Sites, ray_obj_all, canvas_obj, simplices_prune, iPoints);
        Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);
        Voronoi = findSitesCausingDisplacement(Voronoi);
    return Voronoi

def AdaptPositionsWeights(Voronoi):
    K = Voronoi.siteCausingDisp;
    Sites = Voronoi.Sites;
    Weights = Voronoi.Weights;
    for i in range(Voronoi.n_sites):
        if not Voronoi.polygons[i]: #rescue the sites
            Sites[i,:] = RandomPointsInCanvas(1, Voronoi.canvas_obj.xy);
            Weights[i] = Voronoi.canvas_obj.area/Voronoi.n_sites;
        else:
            Sites[i,:] = findCOM(Voronoi.polygons[i].xy); #adapt the position to the centroid
    for i in range(Voronoi.n_sites):
        if Voronoi.polygons[i]:
            if Voronoi.moveMe[i] == True: #speed-up heuristic#
                dist_vec = Sites[i,:] - Sites[K,:];
                dist_abs = math.sqrt(dist_vec[0]**2 + dist_vec[1]**2);
                if dist_abs > 0:
                    ds_abs = (math.sqrt(Voronoi.polygons[i].area/math.pi))**2/(5*dist_abs);
                    Site_new = Sites[i,:] + dist_vec*ds_abs/dist_abs;
                    if PointIsWithinCanvas(Site_new, Voronoi.polygons[i]):
                        Sites[i,:] = Site_new;
            distanceBorder = findMinDistToBorder(Sites[i,:], Voronoi.polygons[i]);
            Weights[i] = (np.nanmin( [math.sqrt(Weights[i]) , distanceBorder] ))**2;
    Voronoi.Sites = Sites;
    Voronoi.Weights = Weights;
    return Voronoi

def AdaptWeights(Voronoi):
    Sites = Voronoi.Sites;
    Weights = Voronoi.Weights;    
    NN_dist = NearestNeighbor(Sites);
    epsilon = Voronoi.canvas_obj.area/1000; #minimum weight of each site
    for i in range(Voronoi.n_sites):
        if not Voronoi.polygons[i]: #rescue the sites
            Sites[i,:] = RandomPointsInCanvas(1, Voronoi.canvas_obj.xy);
            Weights[i] = Voronoi.canvas_obj.area/Voronoi.n_sites;
        else:
            A_current = Voronoi.polygons[i].area;
            A_target = Voronoi.canvas_obj.area*Voronoi.Values[i]/sum(Voronoi.Values);
            f_adapt = A_target/A_current;
            sqrt_w_new = math.sqrt(Weights[i])*f_adapt;
            sqrt_w_max = NN_dist[i];
            Weights[i] = (np.nanmin([sqrt_w_new, sqrt_w_max]))**2;
            Weights[i] = np.nanmax([Weights[i], epsilon]);
    Voronoi.Sites = Sites;
    Voronoi.Weights = Weights;
    return Voronoi

def ComputeError(Voronoi):
    error = 0;
    total_value = sum(Voronoi.Values);
    for i in range(Voronoi.n_sites):
        if not Voronoi.polygons[i]: #rescue the sites if the polygon is empty
            Voronoi.Sites[i,:] = RandomPointsInCanvas(1, Voronoi.canvas_obj.xy);
            Voronoi.Weights[i] = Voronoi.canvas_obj.area/Voronoi.n_sites;
            error += 100;
        else:
            error += np.abs( Voronoi.polygons[i].area - Voronoi.canvas_obj.area*Voronoi.Values[i]/total_value );
    error = error/(2*Voronoi.canvas_obj.area)
    return Voronoi, error

def findSitesCausingDisplacement(Voronoi):
    Voronoi.sitesAreaRatio = [[] for _ in range(Voronoi.n_sites)];
    Voronoi.moveMe = [True for _ in range(Voronoi.n_sites)];    
    for i in range(Voronoi.n_sites):
        if not Voronoi.polygons[i]: #rescue the sites
            Voronoi.Sites[i,:] = RandomPointsInCanvas(1, Voronoi.canvas_obj.xy);
            Voronoi.Weights[i] = Voronoi.canvas_obj.area/Voronoi.n_sites;
        else:
            A_current = Voronoi.polygons[i].area;
            A_target = Voronoi.canvas_obj.area*Voronoi.Values[i]/sum(Voronoi.Values);
            Voronoi.sitesAreaRatio[i] = A_target/A_current;
            if A_current/Voronoi.canvas_obj.area < 0.05:
                Voronoi.moveMe[i] = False;
    try:
        Voronoi.siteCausingDisp = np.argmax(Voronoi.sitesAreaRatio);
        Voronoi.moveMe[Voronoi.siteCausingDisp] = False;
    except ValueError:
        Voronoi.siteCausingDisp = 0; #rescue;
    return Voronoi

def NearestNeighbor(Sites):
    n_sites = len(Sites);
    NN_dist = np.nanmin(distance_matrix(Sites, Sites)+np.diag([np.nan]*n_sites), axis = 0);
    #use np.nanargmin if you want to get the location
    return NN_dist

def Divide2Polygons(Sites, edge, canvas_obj):
    #only for n_site == 2;
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);
    #convert edge to Ax+By = C
    [[x1, y1],[x2, y2]] = edge;
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    #calculate Ax+By-C >0 or <0
    AboveBelow_sites = (np.dot(Sites, np.array([A,B])) - C) >= 0; #each site is above or below the edge
    AboveBelow_corners = (np.dot(canvas_obj.xy, np.array([A,B])) - C) >= 0; #each corner of the canvas is above or below the edge
    Polygons_all = [];
    if AboveBelow_sites[0] != AboveBelow_sites[1]:
        for i in range(n_sites):
            corner_polygon = canvas_obj.xy[AboveBelow_corners == AboveBelow_sites[i]];
            corner_polygon = np.vstack([corner_polygon, edge]);
            corner_polygon = (corner_polygon*(10**9)).astype(int)/(10.**9);
            corner_polygon = np.unique(corner_polygon, axis=0);
            corner_polygon = reorderPoints(corner_polygon);
            polygon = POLYGON(corner_polygon);
            Polygons_all.append(polygon);
    return Polygons_all

def Divide3Polygons(Sites, ray_obj_all, canvas_obj):
    #only for n_sites = 3
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);
    #convert each ray to a_r1*x + b_r1*y + c_r1 = 0
    a_r1, a_r2, a_r3 = ray_obj_all[0].a_r, ray_obj_all[1].a_r, ray_obj_all[2].a_r;
    b_r1, b_r2, b_r3 = ray_obj_all[0].b_r, ray_obj_all[1].b_r, ray_obj_all[2].b_r;
    c_r1, c_r2, c_r3 = ray_obj_all[0].c_r, ray_obj_all[1].c_r, ray_obj_all[2].c_r;
    Coefficient_matrix = np.array([[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]]);
    C_r_matrix = np.array([c_r1, c_r2, c_r3]);
    #calculate the edge created by each ray
    for i in range(3):
        ray_obj_all[i] = FindRayIntersectWithCanvas(ray_obj_all[i], canvas_obj);
    #calculate Ax+By+C >0 or <0 (AB table = above/below table)
    AboveBelow_corners = np.dot(canvas_obj.xy,Coefficient_matrix) + np.tile(C_r_matrix, (n_corners,1)) >= 0; #each corner of the canvas is above or below the edge
    AboveBelow_sites = np.dot(Sites,Coefficient_matrix) + np.tile(C_r_matrix, (n_sites,1)) >= 0; #each site is above or below the edge
    Polygons_all = [[] for _ in range(3)];
    indices = np.arange(3);
    for i in range(3):
        AB_corners_temp = AboveBelow_corners[:,indices!=i];
        AB_site_temp = AboveBelow_sites[i, indices!=i];
        target_corner = (np.equal(AB_corners_temp, AB_site_temp).sum(axis = 1) == 2);
        corner_polygon = canvas_obj.xy[target_corner,:]; #append the included corners by 2 rays
        j_list = np.arange(3); j_list = j_list[j_list!=i];
        for j in j_list:
            if not len(ray_obj_all[j].edge) == 0:
                if len(corner_polygon) == 0:
                    corner_polygon = ray_obj_all[j].edge;
                else:
                    corner_polygon = np.vstack([corner_polygon, ray_obj_all[j].edge]); #append the edges
        if len(corner_polygon) >=3 :
            corner_polygon = (corner_polygon*(10**9)).astype(int)/(10.**9);
            corner_polygon = np.unique(corner_polygon, axis=0);
            corner_polygon = reorderPoints(corner_polygon);
            polygon = POLYGON(corner_polygon);
            Polygons_all[i] = polygon;    
    return Polygons_all

def DivideNPolygons(Sites, ray_obj_all, canvas_obj, simplices_prune, iPoints):
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);
    n_iPoints = len(simplices_prune);
    AB_corners_table = [[] for _ in range(n_iPoints)]; #each corner is above or below the bisector
    AB_sites_table = [[] for _ in range(n_iPoints)]; #each site is above or below the bisector
    for i in range(n_iPoints):
        if not len(ray_obj_all[i]) == 0:
            #convert each ray to a_r1*x + b_r1*y + c_r1 = 0
            a_r1, a_r2, a_r3 = ray_obj_all[i][0].a_r, ray_obj_all[i][1].a_r, ray_obj_all[i][2].a_r;
            b_r1, b_r2, b_r3 = ray_obj_all[i][0].b_r, ray_obj_all[i][1].b_r, ray_obj_all[i][2].b_r;
            c_r1, c_r2, c_r3 = ray_obj_all[i][0].c_r, ray_obj_all[i][1].c_r, ray_obj_all[i][2].c_r;
            Coefficient_matrix = np.array([[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]]);
            C_r_matrix = np.array([c_r1, c_r2, c_r3]);
            #calculate the edge created by each ray
            for j in range(3):
                ray_obj_all[i][j] = FindRayIntersectWithCanvasANDiPoints(ray_obj_all[i][j], canvas_obj, iPoints, simplices_prune);
            #calculate Ax+By+C >0 or <0 (AB table = above/below table)
            AB_corners_table[i] = np.dot(canvas_obj.xy,Coefficient_matrix) + np.tile(C_r_matrix, (n_corners,1)) >= 0; #each corner of the canvas is above or below the edge
            AB_sites_table[i] = np.dot(Sites,Coefficient_matrix) + np.tile(C_r_matrix, (n_sites,1)) >= 0; #each site is above or below the edge
    Polygons_all = [[] for _ in range(n_sites)];
    corner_polygon_all = [[] for _ in range(n_sites)];
    indices = np.arange(3);
    for k in range(n_sites):
        iP_index_required = np.where(np.sum((simplices_prune == k),axis = 1) == 1)[0];
        target_corner_test = [];
        for m in iP_index_required:
            if not len(ray_obj_all[m]) == 0:
                simplex = simplices_prune[m];
                [loc_012] = np.where(simplex == k)[0]
                #finding potential corners:
                AB_corners_temp = AB_corners_table[m][:,indices!=loc_012];
                AB_site_temp = AB_sites_table[m][k, indices!=loc_012];
                target_corner_test.append(set(np.where(np.equal(AB_corners_temp, AB_site_temp).sum(axis = 1) == 2)[0]));
                #append the edges
                n_list = np.arange(3); n_list = n_list[n_list!=loc_012];
                for n in n_list:
                    if not len(ray_obj_all[m][n].edge) == 0:
                        if len(corner_polygon_all[k]) == 0:
                            corner_polygon_all[k] = ray_obj_all[m][n].edge;
                        else:
                            corner_polygon_all[k] = np.vstack([corner_polygon_all[k], ray_obj_all[m][n].edge]);
        #reconcile the potential corners across different iPoints, find target_corner
        if not len(target_corner_test) == 0:
            for p in range(1,len(target_corner_test)):
                target_corner_test[0] = target_corner_test[0].intersection(target_corner_test[p])        
        #append the included corners
        if not len(target_corner_test) == 0:
            if len(corner_polygon_all[k]) == 0:
                corner_polygon_all[k] = canvas_obj.xy[list(target_corner_test[0]),:];
            else:
                corner_polygon_all[k] = np.vstack([corner_polygon_all[k], canvas_obj.xy[list(target_corner_test[0]),:]]);
        if not len(corner_polygon_all[k]) == 0:
            #remove redundant corners
            corner_polygon_all[k] = (corner_polygon_all[k]*(10**9)).astype(int)/(10.**9);
            corner_polygon_all[k] = np.unique(corner_polygon_all[k], axis=0);
            #create and store the polygon object for each site
            if len(corner_polygon_all[k]) >= 3:
                corner_polygon_all[k] = reorderPoints(corner_polygon_all[k]);
                polygon = POLYGON(corner_polygon_all[k]);
                Polygons_all[k]=polygon;     
    return Polygons_all

"""
Auxiliary function
"""
def findMinDistToBorder(site, polygon):
    minDist_list = [];
    for edge in polygon.edges:
        minDist_list.append(findMinDistToEdge(site, edge));
    distanceBorder = min(minDist_list);
    return distanceBorder

def findMinDistToEdge(site, edge):
    [[x1, y1],[x2, y2]] = edge;
    [xs, ys] = site;
    #convert to Ax+By+C=0
    A = (y2-y1); B = -(x2-x1); C = y1*(x2-x1)-x1*(y2-y1);
    dist_perp = np.abs(A*xs + B*ys + C)/math.sqrt(A**2+B**2);
    xp = (B*(B*xs - A*ys) - A*C)/(A**2+B**2);
    yp = (A*(-B*xs + A*ys) - B*C)/(A**2+B**2);
    if isOnSegment([xp, yp], edge):
        minDist = dist_perp;
    else:
        minDist = min(math.sqrt((x1-xs)**2 + (y1-ys)**2), math.sqrt((x2-xs)**2 + (y2-ys)**2));
    return minDist

def findBisectorPositiveRay(Sites, Weights, site_label, canvas_obj):    
    x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
    [[x0], [x1], [x2]] = x_col; [[y0], [y1], [y2]] = y_col; [w0, w1, w2] = Weights;
    A1 = 2*x0; A2 = 2*x1; A3 = 2*x2;
    B1 = 2*y0; B2 = 2*y1; B3 = 2*y2; #C1 = -1; C2 = -1; C3 = -1;
    D1 = -(x0**2+y0**2-w0); D2 = -(x1**2+y1**2-w1); D3 = -(x2**2+y2**2-w2);
    Det = -A1*B2 - A2*B3 - A3*B1 + A3*B2 + A1*B3 + A2*B1;
    while Det == 0:
        Sites = RandomPointsInCanvas(3, canvas_obj.xy);
        x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
        [[x0], [x1], [x2]] = x_col; [[y0], [y1], [y2]] = y_col; [w0, w1, w2] = Weights;
        A1 = 2*x0; A2 = 2*x1; A3 = 2*x2;
        B1 = 2*y0; B2 = 2*y1; B3 = 2*y2; #C1 = -1; C2 = -1; C3 = -1;
        D1 = -(x0**2+y0**2-w0); D2 = -(x1**2+y1**2-w1); D3 = -(x2**2+y2**2-w2);
        Det = -A1*B2 - A2*B3 - A3*B1 + A3*B2 + A1*B3 + A2*B1;
    x_int = (-D1*B2 - D2*B3 - D3*B1 + D3*B2 + D1*B3 + D2*B1)/(-Det);
    y_int = (-A1*D2 - A2*D3 - A3*D1 + A3*D2 + A1*D3 + A2*D1)/(-Det);
    iPoint = np.array([x_int, y_int]);
    while PointIsWithinCanvas(iPoint, canvas_obj) == False:
        Det = 0;
        while Det == 0:
            Sites = RandomPointsInCanvas(3, canvas_obj.xy);
            x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
            [[x0], [x1], [x2]] = x_col; [[y0], [y1], [y2]] = y_col; [w0, w1, w2] = Weights;
            A1 = 2*x0; A2 = 2*x1; A3 = 2*x2;
            B1 = 2*y0; B2 = 2*y1; B3 = 2*y2; #C1 = -1; C2 = -1; C3 = -1;
            D1 = -(x0**2+y0**2-w0); D2 = -(x1**2+y1**2-w1); D3 = -(x2**2+y2**2-w2);
            Det = -A1*B2 - A2*B3 - A3*B1 + A3*B2 + A1*B3 + A2*B1;
        x_int = (-D1*B2 - D2*B3 - D3*B1 + D3*B2 + D1*B3 + D2*B1)/(-Det);
        y_int = (-A1*D2 - A2*D3 - A3*D1 + A3*D2 + A1*D3 + A2*D1)/(-Det);
        iPoint = np.array([x_int, y_int]);
    #tangent vectors
    ray_obj_all = [];
    if PointIsWithinTriangle(iPoint, Sites):
        for k in range(3):
            site_remain = np.arange(3); site_remain = site_remain[site_remain!=k];
            tag1 = site_remain[0]; tag2 = site_remain[1];
            t_bisector = np.array([2*(y_col[tag2][0]-y_col[tag1][0]),2*(x_col[tag1][0]-x_col[tag2][0])]);
            if np.dot(t_bisector,(Sites[k,:] - iPoint)) > 0:
                t_bisector = - t_bisector;
            ray_obj_all.append(RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]));
    else:
        COM_line = (Sites[[1,2,0],:] + Sites[[2,0,1],:])/2;
        point_opposite = np.argmin((COM_line - iPoint)[:,0]**2 + (COM_line - iPoint)[:,1]**2);
        for k in range(3):
            site_remain = np.arange(3); site_remain = site_remain[site_remain!=k];
            tag1 = site_remain[0]; tag2 = site_remain[1];
            t_bisector = np.array([2*(y_col[tag2][0]-y_col[tag1][0]),2*(x_col[tag1][0]-x_col[tag2][0])]);
            if k == point_opposite:
                if np.dot(t_bisector,(Sites[k,:] - iPoint)) > 0:
                    t_bisector = - t_bisector;
                ray_obj_all.append(RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]));     
            else:
                if np.dot(t_bisector,(COM_line[k,:] - iPoint)) < 0:
                    t_bisector = - t_bisector;
                ray_obj_all.append(RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]));    
    return Sites, ray_obj_all

def findMultiplePositiveRay(x_col, y_col, Sites, Weights, iPoints, simplices_prune, canvas_obj):
    n_iPoints = len(simplices_prune);
    ray_obj_all = [[] for _ in range(n_iPoints)];
    for i in range(n_iPoints):
        iPoint = iPoints[i,0:2]; 
        site_label = simplices_prune[i];
        Sites_temp = Sites[simplices_prune[i], :];
        x_col_temp = Sites[simplices_prune[i],0].reshape((-1,1)); 
        y_col_temp = Sites[simplices_prune[i],1].reshape((-1,1));
        ray_obj_iPoint = [];
        if PointIsWithinTriangle(iPoint, Sites_temp):
            for k in range(3):
                site_remain = np.arange(3); site_remain = site_remain[site_remain!=k];
                tag1 = site_remain[0]; tag2 = site_remain[1];
                t_bisector = np.array([2*(y_col_temp[tag2][0]-y_col_temp[tag1][0]),2*(x_col_temp[tag1][0]-x_col_temp[tag2][0])]);
                if np.dot(t_bisector,(Sites_temp[k,:] - iPoint)) > 0:
                    t_bisector = - t_bisector;
                ray_obj = RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]);
                ray_obj.index_iP = i; ray_obj.index_ray = k;
                ray_obj_iPoint.append(ray_obj);
        else:
            COM_line = (Sites_temp[[1,2,0],:] + Sites_temp[[2,0,1],:])/2;
            point_opposite = np.argmin((COM_line - iPoint)[:,0]**2 + (COM_line - iPoint)[:,1]**2);
            for k in range(3):
                site_remain = np.arange(3); site_remain = site_remain[site_remain!=k];
                tag1 = site_remain[0]; tag2 = site_remain[1];
                t_bisector = np.array([2*(y_col_temp[tag2][0]-y_col_temp[tag1][0]),2*(x_col_temp[tag1][0]-x_col_temp[tag2][0])]);
                if k == point_opposite:
                    if np.dot(t_bisector,(Sites_temp[k,:] - iPoint)) > 0:
                        t_bisector = - t_bisector;
                    ray_obj = RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]);
                    ray_obj.index_iP = i; ray_obj.index_ray = k;
                    ray_obj_iPoint.append(ray_obj);     
                else:
                    if np.dot(t_bisector,(COM_line[k,:] - iPoint)) < 0:
                        t_bisector = - t_bisector;
                    ray_obj = RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]);
                    ray_obj.index_iP = i; ray_obj.index_ray = k;
                    ray_obj_iPoint.append(ray_obj);
        ray_obj_all[i] = ray_obj_iPoint;
    return ray_obj_all

def FindRayIntersectWithCanvasANDiPoints(ray_obj, canvas_obj, iPoints, simplices_prune):
    #get the label of each ray object
    index_iP = ray_obj.index_iP;
    index_ray = ray_obj.index_ray;     
    intersect_coordinates = [];
    simplex = simplices_prune[index_iP];
    a = (simplices_prune == simplex[0]) + (simplices_prune == simplex[1]) + (simplices_prune == simplex[2]);
    arr = np.nonzero((a.sum(axis = 1) == 2))[0]
    for i in arr:
        ipoint = iPoints[i, 0:2];
        if isOnRay(ipoint, ray_obj):
            if len(intersect_coordinates) == 0:
                intersect_coordinates = ipoint;
            else:
                intersect_coordinates = np.vstack([intersect_coordinates, ipoint]); 
    #find intersection with canvas
    if len(intersect_coordinates) == 0:
        intersect_TF_list = [];
        for edge_canvas in canvas_obj.edges:
            result, intersectPoint = rayIntersectWithEdge(ray_obj, edge_canvas);
            intersect_TF_list.append(result);
            if (result == True) :
                if len(intersectPoint) == 2 :
                    if len(intersect_coordinates) == 0:
                        intersect_coordinates = intersectPoint;
                    else:
                        intersect_coordinates = np.vstack([intersect_coordinates, intersectPoint]);
    #only keep the nearest one to the origin
    if len(intersect_coordinates) > 0:
        intersect_coordinates = intersect_coordinates.reshape((-1,2));
    if len(intersect_coordinates) > 1:
        intersect_coordinates = keepNearestPointsOnRay(intersect_coordinates, ray_obj);
    #append origin and create an edge
    if len(intersect_coordinates) == 0:
        intersect_coordinates = ray_obj.origin;
    else:
        intersect_coordinates = np.vstack([intersect_coordinates, ray_obj.origin]);
    intersect_coordinates = intersect_coordinates.reshape((-1,2));
    if len(intersect_coordinates) > 0:
        intersect_coordinates = (intersect_coordinates*(10**9)).astype(int)/(10.**9)
        unique_intersect_coordinate = np.unique(intersect_coordinates, axis=0);
        if len(unique_intersect_coordinate) == 2:
            edge = unique_intersect_coordinate;
        else:
            edge = [];
    else:
        edge = [];
    ray_obj.edge = edge;
    return ray_obj

def keepNearestPointsOnRay(intersect_coordinates, ray_obj):
    n_candidates = len(intersect_coordinates);
    #find candidates with shortest distance from the origin;
    dist = np.sum((intersect_coordinates - np.tile(ray_obj.origin,(n_candidates,1)))**2, axis = 1);
    dist = (dist*(10**9)).astype(int)/(10.**9);
    dist[dist == 0] = np.inf; #avoid the intersection point that is the origin itself
    intersect_coordinates = intersect_coordinates[np.argmin(dist), :];
    return intersect_coordinates

def lineIntersectWithEdge(bisector, edge_canvas):
    [P1, P2] = edge_canvas;
    [[x1, y1],[x2, y2]] = edge_canvas;
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    [D, E, F] = bisector;
    Det = A*E - B*D; #np.linalg.det(np.array([[A,B],[D,E]]));
    DetX = C*E - F*B; #np.linalg.det(np.array([[C,B],[F,E]]));
    DetY = A*F - C*D; #np.linalg.det(np.array([[A,C],[D,F]]));
    #if their slopes are the same but they don't intersect, they are parallel
    if Det == 0:
        if (DetX == 0) & (DetY == 0):
            result = True;
            intersectPoint = [np.inf];
        else:
            result = False;
            intersectPoint = [];
    else:
        x_int = DetX/Det;
        y_int = DetY/Det;
        if isOnSegment([x_int, y_int], edge_canvas):
            result = True;
            intersectPoint = np.array([x_int, y_int]);
        else:
            result = False;
            intersectPoint = [];
    return result, intersectPoint

def findIntersectWithCanvas(bisector, canvas_obj):
    n_corners = canvas_obj.n_corners;
    intersect_TF_list = [];
    intersect_coordinates = [];
    for i in range(n_corners):
        edge_canvas = canvas_obj.edges[i];
        result, intersectPoint = lineIntersectWithEdge(bisector, edge_canvas);
        intersect_TF_list.append(result);
        if (result == True):
            if len(intersectPoint) == 2:
                if len(intersect_coordinates) == 0:
                    intersect_coordinates = intersectPoint;
                else:
                    intersect_coordinates = np.vstack([intersect_coordinates, intersectPoint]);
    if len(intersect_coordinates) > 0:
        intersect_coordinates = intersect_coordinates.reshape((-1,2));
        intersect_coordinates = (intersect_coordinates*(10**9)).astype(int)/(10.**9);
        unique_intersect_coordinate = np.unique(intersect_coordinates, axis=0);
        if len(unique_intersect_coordinate) == 2:
            edge = unique_intersect_coordinate;
        else:
            edge = [];
    else:
        edge = [];
    return edge

def PruneAndConvertToIPoints(Dual_Sites, simplices, canvas_obj):
    COM_DualSites = findCOM(Dual_Sites);
    simplices_prune = simplices;
    iPoints = np.array([]);
    n_faces = len(simplices);
    prune_list = [];
    for i in range(n_faces):
        simplex = simplices[i];
        x1, y1, z1 = Dual_Sites[simplex[0]];
        x2, y2, z2 = Dual_Sites[simplex[1]];
        x3, y3, z3 = Dual_Sites[simplex[2]];
        alpha = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
        beta = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
        gamma = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
        delta = x1*(y2*z3-y3*z2) + x2*(y3*z1-y1*z3) + x3*(y1*z2-y2*z1);
        a = - alpha/gamma;
        b = - beta/gamma;
        c = delta/gamma; #z = ax+by+c => n = (a, b, -1)
        #determine if the face belongs to the lower convex Hull
        COM_face = findCOM(Dual_Sites[simplex])
        if np.dot([a, b, -1], (COM_face - COM_DualSites))>0:
            #nf = np.array([a,b,-1]);
            #determine if the iPoint is within canvas
            if PointIsWithinCanvas(np.array([a/2, b/2]), canvas_obj):
                iPoints = np.append(iPoints, np.array([a/2, b/2, -c]), axis=0)
            else:
                prune_list.append(i);
        else:
            #nf = np.array([-a,-b,1]);
            prune_list.append(i);
    iPoints = iPoints.reshape(-1,3);
    simplices_prune = np.delete(simplices_prune, prune_list, axis = 0);
    return iPoints, simplices_prune

def isIntersect(edge1, edge2):
    #both edges are 2x2 np array
    [P1, P2] = edge1; [P3, P4] = edge2;
    [[x1, y1],[x2, y2]] = edge1
    [[x3, y3],[x4, y4]] = edge2
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    D = (y4-y3); E = -(x4-x3); F = x3*(y4-y3)-y3*(x4-x3);
    Det = A*E - B*D;#np.linalg.det(np.array([[A,B],[D,E]]));
    DetX = C*E - B*F;#np.linalg.det(np.array([[C,B],[F,E]]));
    DetY = A*F - C*D;#np.linalg.det(np.array([[A,C],[D,F]]));
    #if their slopes are the same but they don't intersect, they are parallel
    if Det == 0:
        if (DetX == 0) & (DetY == 0):
            if isOnSegment(P1, edge2) | isOnSegment(P2, edge2):
                if len(np.unique(np.vstack([edge1, edge2]), axis=0)) == 3:
                    result = True;
                    intersectPoint = P1*np.array_equal(P1, P3) + P1*np.array_equal(P1, P4) + P2*np.array_equal(P2, P3) + P2*np.array_equal(P2, P4)
                else:
                    result = True;
                    intersectPoint = [np.inf];
            else:
                result = False;
                intersectPoint = [];        
        else: #(Det == 0) & ((DetX != 0) | (DetY != 0))
            result = False;
            intersectPoint = []; 
    else:
        x_int = DetX/Det;
        y_int = DetY/Det;
        if isOnSegment([x_int, y_int], edge1):
            if isOnSegment([x_int, y_int], edge2):
                result = True;
                intersectPoint = np.array([x_int, y_int]);
            else:
                result = False;
                intersectPoint = [];
        else:
            result = False;
            intersectPoint = [];
    return result, intersectPoint

def reorderPoints(corners):
    #input: numpy array of nx2
    orderedPoints = corners;
    COM = orderedPoints.mean(axis = 0);#find center of mass
    orderedPoints=orderedPoints[np.argsort(np.angle((orderedPoints - COM)[:,0] + (orderedPoints - COM)[:,1]*1j))]
    return orderedPoints

def reorderPointsAndSimplices(corners, simplices):
    #input: numpy array of nx2
    orderedPoints = corners;
    COM = orderedPoints.mean(axis = 0);#find center of mass
    order = np.argsort(np.angle((orderedPoints - COM)[:,0] + (orderedPoints - COM)[:,1]*1j));
    orderedPoints=orderedPoints[order];
    orderedSimplices = simplices[order];
    return orderedPoints, orderedSimplices

def RandomPointsInCanvas(n_Points, canvas):
    #divide the canvas into multiple triangels
    #weight each triangle according to the area
    #randomly placed points within the triangle
    #Please make sure that the points have been reordered properly;
    vec_all = canvas[1:,:] - canvas[0,:];
    n_triangle = len(canvas) - 2;
    area_triangle = np.zeros(n_triangle);
    for i in range(n_triangle):
        area_triangle[i] = PolygonArea(np.vstack([[0,0],vec_all[i],vec_all[i+1]]));
    randScale = np.sum(np.tril(area_triangle), axis = 1)/sum(area_triangle);
    randNum = np.hstack([np.random.rand(n_Points, 3), np.zeros([n_Points,1])]);
    Sites = np.zeros([n_Points, 2]);
    for i in range(n_triangle):
        mask = (randNum[:,0] <= randScale[i]) & (randNum[:,3] == 0);
        vec1_tile_xy = np.tile(vec_all[i,:],(sum(mask),1));
        vec2_tile_xy = np.tile(vec_all[i+1,:],(sum(mask),1));    
        randNum_masked = randNum[mask,1].reshape((-1, 1));
        randLen_masked = np.sqrt(randNum[mask,2].reshape((-1, 1)));
        Sites[mask,:]=randNum_masked*vec1_tile_xy*randLen_masked+(1-randNum_masked)*vec2_tile_xy*randLen_masked;
        randNum[mask,3]=1
    Sites = Sites + np.tile(canvas[0,:], (n_Points,1));
    return Sites

def PolygonArea(corners):
    #calculate polygon area using shoelace formula
    #please make sure that the corners are reordered before calling PolygonArea function!
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def findCOM(corners):
    COM = corners.mean(axis = 0);
    return COM

def isOnSegment(P, edge):
    [P1, P2] = edge
    [[x1, y1],[x2, y2]] = edge
    [x, y] = P
    #convert to Ax+By=C
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1) - y1*(x2-x1);
    if math.sqrt(A**2+B**2) == 0:
        if (x == x1) & (y == y1):
            result = True;
        else:
            result = False
    else:
        tester = (A*x+B*y-C);
        if int(tester*(10**9))/(10.**9) == 0:
            x = int(x*(10**9)+0.5)/(10.**9); x1 = int(x1*(10**9)+0.5)/(10.**9); x2 = int(x2*(10**9)+0.5)/(10.**9);
            y = int(y*(10**9)+0.5)/(10.**9); y1 = int(y1*(10**9)+0.5)/(10.**9); y2 = int(y2*(10**9)+0.5)/(10.**9);
            if (x >= min(x1, x2)) & (x <= max(x1, x2)) & (y >= min(y1, y2)) & (y <= max(y1, y2)):
                result = True;
            else:
                result = False;
        else:
            result = False;        
    return result

def isOnRay(P, ray_obj):
    #P: 1x2 np array;  edge: 2x2 np array
    #convert to Ax+By=C
    A = ray_obj.a_r; B = ray_obj.b_r; C = -ray_obj.c_r;
    t_x = ray_obj.tangent[0]; t_y = ray_obj.tangent[1];
    test = (A*P[0]+B*P[1]-C);
    if math.isnan(test):
        result = False;
    else:
        if (int(test*(10**9))/(10.**9) == 0):
            if t_x != 0:
                if int((P[0]-ray_obj.origin[0])/t_x*(10**9))/(10.**9) >=0:
                    result = True;
                else:
                    result = False;
            else:
                if int((P[1]-ray_obj.origin[1])/t_y*(10**9))/(10.**9) >=0:  
                    result = True;
                else:
                    result = False;
        else:
            result = False;
    return result

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def angle(v1, v2):
    cosine = dotproduct(v1, v2) / (math.hypot(v1[0],v1[1]) * math.hypot(v2[0],v2[1]));
    if (cosine > 1)|(cosine < -1):
        return 100
    else:
        return math.acos(cosine)

def PointIsOnCornerOfCanvas(P, canvas_obj):
    a = (canvas_obj.xy == P)
    if len(np.nonzero(a[:,0] & a[:,1])[0]) == 1:
        return True
    else:
        return False

def PointIsWithinCanvas(P, canvas_obj):
    if PointIsOnCornerOfCanvas(P, canvas_obj):
        return True;
    else:
        sum_angle = 0;
        for i in range(canvas_obj.n_corners):
            edge_canvas = canvas_obj.edges[i];
            a = edge_canvas - P;
            sum_angle += angle(a[0,:], a[1,:]);
        if not math.isnan(sum_angle):
            if int((sum_angle - 2*math.pi)*(10**9))/(10.**9) == 0:
                return True
            else:
                return False
        else:
            return False

def PointIsWithinTriangle(P, triangle):
    a = (triangle == P)
    if len(np.nonzero(a[:,0] & a[:,1])[0]) == 1:
        return true
    else:
        sum_angle = 0;
        for i in range(3):
            edge = triangle[[i-1,i],:];
            a = edge - P;
            sum_angle += angle(a[0,:], a[1,:]);
        if not math.isnan(sum_angle):
            if int((sum_angle - 2*math.pi)*(10**9))/(10.**9) == 0:
                return True;
            else:
                return False;
        else:
            return False;

def EdgeCutWithinCanvas(edge, canvas_obj):
    [P1, P2] = edge
    P1_In = PointIsWithinCanvas(P1, canvas_obj);
    P2_In = PointIsWithinCanvas(P2, canvas_obj);
    if (P1_In == True) & (P2_In == True):
        edge_cut = edge;
    elif (P1_In == False) & (P2_In == False):
        edge_cut = [];
    else:
        for i in range(canvas_obj.n_corners):
            edge_canvas = canvas_obj.edges[i];
            result, intersectPoint = isIntersect(edge, edge_canvas);
            if result == True:
                P_In = P1*P1_In + P2*P2_In
                edge_cut = np.vstack([P_In, intersectPoint]);
    return edge_cut

def FindRayIntersectWithCanvas(ray_obj, canvas_obj):
    intersect_TF_list = [];
    intersect_coordinates = [];
    for edge_canvas in canvas_obj.edges:
        result, intersectPoint = rayIntersectWithEdge(ray_obj, edge_canvas);
        intersect_TF_list.append(result);
        if result == True: 
            if len(intersectPoint) == 2:
                if len(intersect_coordinates) == 0:
                    intersect_coordinates = intersectPoint;
                else:
                    intersect_coordinates = np.vstack([intersect_coordinates, intersectPoint]);
    if len(intersect_coordinates) == 0:
        intersect_coordinates = ray_obj.origin;
    else:
        intersect_coordinates = np.vstack([intersect_coordinates, ray_obj.origin]);
    if len(intersect_coordinates) > 0:
        intersect_coordinates = intersect_coordinates.reshape((-1,2));
        intersect_coordinates = (intersect_coordinates*(10**9)).astype(int)/(10.**9);
        unique_intersect_coordinate = np.unique(intersect_coordinates, axis=0);
        if len(unique_intersect_coordinate) == 2:
            edge = unique_intersect_coordinate;
        else:
            edge = [];
    else:
        edge = [];
    ray_obj.edge = edge;
    return ray_obj

def rayIntersectWithEdge(ray_obj, edge_canvas):
    [P1, P2] = edge_canvas;
    [[x1, y1],[x2, y2]] = edge_canvas;
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    D, E, F = ray_obj.a_r, ray_obj.b_r, -ray_obj.c_r;
    Det = A*E - B*D; #np.linalg.det(np.array([[A,B],[D,E]]));
    DetX = C*E - B*F; #np.linalg.det(np.array([[C,B],[F,E]]));
    DetY = A*F - C*D;#np.linalg.det(np.array([[A,C],[D,F]]));
    if Det == 0:
        if ((DetX != 0) | (DetY != 0)):
            result = False;
            intersectPoint = [];
        else: #(Det == 0) & (DetX == 0) & (DetY == 0)
            P1isOnRay = isOnRay(P1, ray_obj); P2isOnRay = isOnRay(P2, ray_obj);
            OriginIsOnSeg = isOnSegment(ray_obj.origin, edge_canvas);
            if P1isOnRay & P2isOnRay:
                result = True;
                intersectPoint = [np.inf];
            elif (P1isOnRay != P2isOnRay):
                if np.array_equal(P1, ray_obj.origin) | np.array_equal(P2, ray_obj.origin):
                    result = True;
                    intersectPoint = P1*np.array_equal(P1, ray_obj.origin) + P2*np.array_equal(P2, ray_obj.origin);
                else:
                    result = True;
                    intersectPoint = [np.inf];
            else:
                result = False;
                intersectPoint = [];               
    else:
        x_int = DetX/Det;
        y_int = DetY/Det;
        if isOnSegment([x_int, y_int], edge_canvas):
            if isOnRay(np.array([x_int, y_int]), ray_obj):
                result = True;
                intersectPoint = np.array([x_int, y_int]);
            else:
                result = False; intersectPoint = [];
        else:
            result = False;
            intersectPoint = [];    
    return result, intersectPoint
"""
Timer
"""
def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference
TicToc = TicTocGenerator() # create an instance of the TicTocGen generator
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

"""
main program
"""
if __name__ == "__main__":
    import string
    import random
    def randomString(stringLength):
        """Generate a random string of fixed length """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))
    canvas = np.array([[0,0],[4,0],[4,4],[0,4]]);
    i_max = 75;
    err_thres = 1E-6;
    #create dataframe
    IDs = np.array( [randomString(3) for _ in range(160)] )
    index1 = np.array(('P1 '*32+'P2 '*32+'P3 '*32+'M1 '*32+'M2 '*32).split( ));
    index0 = np.array(('protein '*96+'mRNA '*64).split( ));
    values = np.random.rand(160, 1);
    df = pd.DataFrame(values.reshape((-1,1)), index = [index0, index1, IDs], columns = ['values'])
    df = df.sort_index();
    Voronoi_0, Voronoi_1_all, Voronoi_2_all = Layered_Voronoi(df, canvas, i_max, err_thres);
    error_all = Layered_Voronoi_plot(Voronoi_0, Voronoi_1_all, Voronoi_2_all);