"""
Plot Voronoi diagram

@author: Ray Chang
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019

Algorithm Reference: A. Nocaj and U. Brandes, "Computing Voronoi Treemaps: Faster, Simpler, and 
Resolution-independent, Computer Graphics Forum, vol. 31, no. 31, pp. 855-864, 2012. 
Available: 10.1111/j.1467-8659.2012.03078.x.
"""

"""
This program generate a Voronoi diagram, in which the canvas of a plot is divided in to n polygons.
Each data point (or "Points" in this program) corresponds to one polygon, and the centroid of each
polygon is called a "site". The area of each polygon is proportional to the value each point wants 
to represent. For example, it can be the mass compositions of a cell, the relative abundance of 
different proteins, or the energy cost of different cellular processes.

Dividing a canvas (in this program we use np.array([[0,0],[4,0],[4,4],[0,4]])) into multiple polygons
is a complicated computational geometry problem. The main working horse of this task is the function
"Voronoi_treemap". In "Voronoi_treemap", we did the following things:
(1) initialize random sites within the canvas (using "RandomPointsInCanvas") and assign an initial
weight to each site. The initial weights of each site are set equal so that it will be easier to 
compute a draft.
(2) compute the draft of our Voronoi diagram, which is called "Power Diagram" following the notation
of our algorithm reference. (using "ComputePowerDiagram")
(3) adjust the positions of each site and the weights of each site until the program reaches the
maximum iterations (i_max = 75) or until error < err_thres.
(4) return the final error of the whole Voronoi diagram.

In Voronoi diagram, certain amount of error in area representation is expected. According to the 
analysis in our algorithm reference, an error between 5~23% is expected. However, an error up to 20%
can severely distort the diagram. We therefore set a stricter limit (15%) on our error in the function
"Voronoi_Main_functions". The program will be trapped in a while loop if the final error of the whole
diagram is greater than 15%. If the newly computed Voronoi diagram has a lower error compared to the
previous one, it will continuously increase the number of iterations until the error falls below 15%.
However, if the error of a newly computed diagram is greater than the previous one, the program will
plot the diagram all over again. Generally, a final error <=10% will give you a nice representation.

In this program, we use 3 types of objects to facilitate the calculation. A POLYGON stores the coordinates 
of the corners, the edges, and the area of the polygon. A VORONOI stores all the polygons of a voronoi
diagram, the coordinates/weights/values of each site, and the canvas of the whole plot. A RAY object 
stores the origin and the tangent vector of a ray, and the site it represents.

Using this algorithm, the program is capable of dividing 32 polygons in 10 seconds, and the program 
scales with O(n*log(n)). As a result, it would be preferable to introduce layered structure if the number
of polygons are more than 50. To fulfill this need, we use function "Voronoi_Main_functions_for_layer"
and "Layered_Voronoi" to accomplish this need.

If you directly run this program, it will give you a demo Voronoi plot with 32x5 polygons. If you want to 
implement the functions that are created here and make a new Voronoi plot with your own dataset, please 
consider the following:
(1)Place this in the heading, or other ways that can import every function in this file.
import sys
sys.path.append('wholecell/utils')
from voronoiPlotMain import *

(2)Determine whether your data is layered/hiearchical or not.
- If your data is not layered/hiearchical, please appropriately modify the following code to meet your 
needs:

canvas = np.array([[0,0],[4,0],[4,4],[0,4]]); #you can use other shape as well.
Points = [THE IDS OF YOUR DATA POINTS];
Values = np.array([THE VALUES OF EACH DATA POINTS, IT COULD BE ABUNDANCE, ENERGY COST, FLUX,...]);
error = Voronoi_Main_functions(Points, Values, canvas);

- If your data is layered/hiearchical, please store your data into a pandas dataframe, and appropriately
modify the following code to meet your needs:

canvas = np.array([[0,0],[4,0],[4,4],[0,4]]); #you can use other shape as well.
i_max = 75;
err_thres = 1E-6;
IDs = np.array([THE IDS OF YOUR DATA POINTS])
index1 = np.array([THE SECOND LAYER HIEARCHY OF EACH POINT: E.G. HOUSEKEEPING, METABOLISM,...]);
index0 = np.array([THE FIRST LAYER HIEARCHY OF EACH POINT: E.G. PROTEINS, mRNAs, ...]);
values = np.array([THE VALUES OF EACH DATA POINTS, IT COULD BE ABUNDANCE, ENERGY COST, FLUX,...]);
df = pd.DataFrame(values.reshape((-1,1)), index = [index0, index1, IDs], columns = ['values'])
df = df.sort_index();
Voronoi_0, Voronoi_1_all, Voronoi_2_all = Layered_Voronoi(df, canvas, i_max, err_thres);
error_all = Layered_Voronoi_plot(Voronoi_0, Voronoi_1_all, Voronoi_2_all);
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
    """ 
    Call Voronoi_treemap to compute the Voronoi diagram, with baseline number of iterations = 75. 
    If the computed Voronoi diagram has error >= 0.15, increase the iterations to 150 and so on
    until the error < 0.15. If the error increases after increasing the number of iterations, 
    recompute the whole Voronoi diagram all over again.

    Please be careful when interpreting the error output. An error between 0.5~0.15 is expected.
    If the error = 0, that means the whole Voronoi diagram is not properly generated.
    """
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
    """
    The same Voronoi main function as Voronoi_Main_functions but designed for layered structures.
    The only differences is that this function will not generate plot.
    """
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
    """
    Main function used for computing layered Voronoi diagram. This function takes a pandas dataframe
    and canvas as input and generate a layered Voronoi diagram.
    """
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
    """
    Generate the layered Voronoi diagram from the output of Layered_Voronoi and compute the overall
    error of the final plot.
    """
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

            # plot Voronoi Treemap
            canvas_obj = Voronoi_2_all[k][m].canvas_obj; 
            Sites = Voronoi_2_all[k][m].Sites; 
            Values = Voronoi_2_all[k][m].Values; 
            Points = Voronoi_2_all[k][m].Points;

            # plot canvas
            for edge_canvas in canvas_obj.edges:
                ax.plot(edge_canvas[:,0],edge_canvas[:,1],'black', lw=1.5, solid_capstyle='round', zorder=2);

            # plot polygons
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

            # add text label to each polygon
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
    """
    The main working horse of the whole function. This function compute the initial Voronoi diagram,
    adapt the location of each site and the weight of each site, re-compute the Voronoi diagram, and
    compute the error of the Voronoi diagram. If the error is too high, repeat this process until the
    error falls below error threshold or the number of iteration exceeds i_max(= 75).

    Certain amount of error is intrinsic to this plot. Therefore, blindly increasing the number of maximum
    iteration does not guarantee a decrease in error.
    """
    n_Points = len(Points);

    if n_Points == 1:
        # if there is only 1 point, we don't need to compute the diagram: the point occupies the whole canvas
        Voronoi = VORONOI([canvas_obj], findCOM(canvas_obj.xy).reshape((-1,2)), canvas_obj.area, Values, canvas_obj);
        error = 0.0000001;
        Voronoi.Points = Points;
        return Voronoi, error

    else:
        # initialization
        Sites = RandomPointsInCanvas(n_Points, canvas_obj.xy);
        Weights = np.ones(n_Points)*canvas_obj.area/n_Points;
        error = float("inf");

        # calculate initial Voronoi plot
        Voronoi = ComputePowerDiagram(canvas_obj, Sites, Weights, Values);

        # adjust position & weight until error becomes acceptable
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
    """
    Basically the same function as Voronoi_treemap but skip the initialization part. This ensures
    the new plot is built on previous results and the previous iterations are not wasted.
    """
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
    """
    Plot the Voronoi diagram for non-layered structure.
    """
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
def ComputePowerDiagram(canvas_obj, Sites, Weights, Values):
    """
    Compute the draft of the Voronoi diagram, or the power diagram following the nomeclature of our
    code reference. This function seperate the cases for n = 2, 3, 4 and above. For n=2, since there
    should always be a solution, the function will rescue itself by resetting the initial sites until
    a solution is found.
    """
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
        
        #4. compute the ray objects of each iPoints
        ray_obj_all = findMultiplePositiveRay(x_col, y_col, Sites, Weights, iPoints, simplices_prune, canvas_obj);
        
        #5. compute the coordinate of each polygon
        Polygons_all = DivideNPolygons(Sites, ray_obj_all, canvas_obj, simplices_prune, iPoints);     
        Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);

        #6. find the site with largest Area_target/Area_current ratio, which will displace its neighbor
        Voronoi = findSitesCausingDisplacement(Voronoi);

    return Voronoi

def reComputePowerDiagram(Voronoi):
    """
    Basically the same function as ComputePowerDiagram but accept an exisiting Voronoi object as input.
    """
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

        #4. compute the ray objects of each iPoints
        ray_obj_all = findMultiplePositiveRay(x_col, y_col, Sites, Weights, iPoints, simplices_prune, canvas_obj);

        #5. compute the coordinate of each polygon
        Polygons_all = DivideNPolygons(Sites, ray_obj_all, canvas_obj, simplices_prune, iPoints);
        Voronoi = VORONOI(Polygons_all, Sites, Weights, Values, canvas_obj);

        #6. find the site with largest Area_target/Area_current ratio, which will displace its neighbor
        Voronoi = findSitesCausingDisplacement(Voronoi);

    return Voronoi

def AdaptPositionsWeights(Voronoi):
    """
    Adapt the positions & weights of each site within an existing Voronoi diagram.
    This algorithm also considers the point with largest Area_{target}/Area_{current} ratio
    and displace its neighbor in order to speed up the optimization process.
    """
    K = Voronoi.siteCausingDisp;
    Sites = Voronoi.Sites;
    Weights = Voronoi.Weights;
    #1. move the location of each site to the centroid of each polygon. If there is no closed
    #polygon for that site, randomly reassign a location for that site in order to rescue it.
    for i in range(Voronoi.n_sites):
        if not Voronoi.polygons[i]: #rescue the sites
            Sites[i,:] = RandomPointsInCanvas(1, Voronoi.canvas_obj.xy);
            Weights[i] = Voronoi.canvas_obj.area/Voronoi.n_sites;

        else:
            Sites[i,:] = findCOM(Voronoi.polygons[i].xy); #adapt the position to the centroid

    #2. speed up heuristic: if a site is expected to expand its area greatly, its neighbor should
    #be displaced in advance in order to speed up the optimization process. 
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

    #3. if the square root of weights exceed the minimum distance from the site to the current border
    #of the polygon, the weight should be decreased. 
            distanceBorder = findMinDistToBorder(Sites[i,:], Voronoi.polygons[i]);
            Weights[i] = (np.nanmin( [math.sqrt(Weights[i]) , distanceBorder] ))**2;

    Voronoi.Sites = Sites;
    Voronoi.Weights = Weights;
    return Voronoi

def AdaptWeights(Voronoi):
    """
    Adapt the weight of each site. increase the weight if current area is different from the targeted
    area of each site. However, if the weight is too big such that the site is encroaching its neighbor,
    The weight should be decreased.
    """
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
    """
    compute the error of the Voronoi plot. The error is defined as sum(abs(Area_current - Area_target))/(2*Area_canvas)
    """
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
    """
    If a site need to expand its area greatly in the next round, the other sites surrounding it should
    be displaced a little bit in order to facilitate the process of optimizing the Voronoi diagram.
    This is part of the speed-up heuristic.
    """
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
    """
    find the nearest neighbor of each site and return the DISTANCE (rather than the index) between
    each site and its nearest neighbor
    """
    n_sites = len(Sites);
    NN_dist = np.nanmin(distance_matrix(Sites, Sites)+np.diag([np.nan]*n_sites), axis = 0);
    return NN_dist

def Divide2Polygons(Sites, edge, canvas_obj):
    """
    This function only works for n_site = 2. If the edge formed by the bisector between 2 sites and
    the canvas is computed already, This function returns the coordinates of the 2 polygons formed by
    the canvas & the edge.
    """
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);

    #1. convert edge to Ax+By = C
    [[x1, y1],[x2, y2]] = edge;
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);

    #2. calculate for each site & each corner of the canvas, Ax+By-C >0 or <0
    #to determine if they are on the same side to the edge.
    AboveBelow_sites = (np.dot(Sites, np.array([A,B])) - C) >= 0; #each site is above or below the edge
    AboveBelow_corners = (np.dot(canvas_obj.xy, np.array([A,B])) - C) >= 0; #each corner of the canvas is above or below the edge
    Polygons_all = [];

    #3. assign the corners which is on the same side as each site, and determine the final coordinates
    #of the 2 polygons.
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
    """
    This function only works for n_sites = 3. If the intersection point of the 3 bisectors and the 3 rays
    are computed already, this function find where these 3 rays intersect with canvas, and divide the canvas
    into 3 polygons.
    """
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);

    #0. convert each ray to a_r1*x + b_r1*y + c_r1 = 0
    a_r1, a_r2, a_r3 = ray_obj_all[0].a_r, ray_obj_all[1].a_r, ray_obj_all[2].a_r;
    b_r1, b_r2, b_r3 = ray_obj_all[0].b_r, ray_obj_all[1].b_r, ray_obj_all[2].b_r;
    c_r1, c_r2, c_r3 = ray_obj_all[0].c_r, ray_obj_all[1].c_r, ray_obj_all[2].c_r;
    Coefficient_matrix = np.array([[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]]);
    C_r_matrix = np.array([c_r1, c_r2, c_r3]);

    #1. calculate the edge created by each ray by finding the intersection point of each ray with canvas
    for i in range(3):
        ray_obj_all[i] = FindRayIntersectWithCanvas(ray_obj_all[i], canvas_obj);
    
    #2. determine if each site & each corner is above or below each ray by calculating Ax+By+C >0 or <0 
    #(AB table = above/below table)
    AboveBelow_corners = np.dot(canvas_obj.xy,Coefficient_matrix) + np.tile(C_r_matrix, (n_corners,1)) >= 0; #each corner of the canvas is above or below the edge
    AboveBelow_sites = np.dot(Sites,Coefficient_matrix) + np.tile(C_r_matrix, (n_sites,1)) >= 0; #each site is above or below the edge

    #3. group the corners, edges, intersection points that are on the same side as each site
    #and form the coordinates of each polygon.
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
    """
    If all the intersection points and their corresponding 3 rays are computed already, this function
    find the endpoint of each ray, and groups them into n polygons. The endpoint of each ray could be
    another intersection point or some point on canvas.
    """
    n_corners = canvas_obj.n_corners;
    n_sites = len(Sites);
    n_iPoints = len(simplices_prune);

    #1. for all the intersection points, we first determine whether each site & each corners are above
    #or below each ray.
    AB_corners_table = [[] for _ in range(n_iPoints)];
    AB_sites_table = [[] for _ in range(n_iPoints)];
    for i in range(n_iPoints):
        if not len(ray_obj_all[i]) == 0:
            #1-1. convert each ray to a_r1*x + b_r1*y + c_r1 = 0
            a_r1, a_r2, a_r3 = ray_obj_all[i][0].a_r, ray_obj_all[i][1].a_r, ray_obj_all[i][2].a_r;
            b_r1, b_r2, b_r3 = ray_obj_all[i][0].b_r, ray_obj_all[i][1].b_r, ray_obj_all[i][2].b_r;
            c_r1, c_r2, c_r3 = ray_obj_all[i][0].c_r, ray_obj_all[i][1].c_r, ray_obj_all[i][2].c_r;
            Coefficient_matrix = np.array([[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]]);
            C_r_matrix = np.array([c_r1, c_r2, c_r3]);

            #1-2. calculate the edge created by each ray
            for j in range(3):
                ray_obj_all[i][j] = FindRayIntersectWithCanvasANDiPoints(ray_obj_all[i][j], canvas_obj, iPoints, simplices_prune);
            
            #1-3. calculate Ax+By+C >0 or <0 (AB table = above/below table)
            AB_corners_table[i] = np.dot(canvas_obj.xy,Coefficient_matrix) + np.tile(C_r_matrix, (n_corners,1)) >= 0; #each corner of the canvas is above or below the edge
            AB_sites_table[i] = np.dot(Sites,Coefficient_matrix) + np.tile(C_r_matrix, (n_sites,1)) >= 0; #each site is above or below the edge

    #2. by referring to the AB_corners_table & AB_sites_table, find the coordinates of each polygon.
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
                #2-1. find the potential corners if they are on the same side as each site:
                AB_corners_temp = AB_corners_table[m][:,indices!=loc_012];
                AB_site_temp = AB_sites_table[m][k, indices!=loc_012];
                target_corner_test.append(set(np.where(np.equal(AB_corners_temp, AB_site_temp).sum(axis = 1) == 2)[0]));

                #2-2. append the edges formed by nearby rays.
                n_list = np.arange(3); n_list = n_list[n_list!=loc_012];
                for n in n_list:
                    if not len(ray_obj_all[m][n].edge) == 0:
                        if len(corner_polygon_all[k]) == 0:
                            corner_polygon_all[k] = ray_obj_all[m][n].edge;
                        else:
                            corner_polygon_all[k] = np.vstack([corner_polygon_all[k], ray_obj_all[m][n].edge]);

        #2-3. reconcile the potential corners across different iPoints, find target_corner
        #(some corners may be on the same side as the site with respect to one iPoint, but not other
        # nearby iPoints. We have to make sure that the final corners included are on the same side
        # as each site with respect to all iPoints nearby.)
        if not len(target_corner_test) == 0:
            for p in range(1,len(target_corner_test)):
                target_corner_test[0] = target_corner_test[0].intersection(target_corner_test[p]) 

        #2-4. append the target corners
        if not len(target_corner_test) == 0:
            if len(corner_polygon_all[k]) == 0:
                corner_polygon_all[k] = canvas_obj.xy[list(target_corner_test[0]),:];
            else:
                corner_polygon_all[k] = np.vstack([corner_polygon_all[k], canvas_obj.xy[list(target_corner_test[0]),:]]);

        if not len(corner_polygon_all[k]) == 0:
            #2-5. remove redundant corners of each polygon
            corner_polygon_all[k] = (corner_polygon_all[k]*(10**9)).astype(int)/(10.**9);
            corner_polygon_all[k] = np.unique(corner_polygon_all[k], axis=0);

            #2-6. create and store the polygon object for each site
            if len(corner_polygon_all[k]) >= 3:
                corner_polygon_all[k] = reorderPoints(corner_polygon_all[k]);
                polygon = POLYGON(corner_polygon_all[k]);
                Polygons_all[k]=polygon;

    return Polygons_all

"""
Auxiliary function
"""
def findMinDistToBorder(site, polygon):
    """
    find the minimum distance of the site to the borders of its polygon.
    """
    minDist_list = [];
    for edge in polygon.edges:
        minDist_list.append(findMinDistToEdge(site, edge));
    distanceBorder = min(minDist_list);
    return distanceBorder

def findMinDistToEdge(site, edge):
    """
    find the minimum distance of the site to one specific edge of a polygon.
    """
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
    """
    This function is designed for n=3. This function compute the 3 bisectors among 3 sites,
    find the intersect point (iPoint) of these 3 bisectors, and define the positive direction
    of each ray.
    """
    #1. Find the 3 bisectors and the iPoint
    x_col = Sites[:,0].reshape((-1,1)); y_col = Sites[:,1].reshape((-1,1));
    [[x0], [x1], [x2]] = x_col; [[y0], [y1], [y2]] = y_col; [w0, w1, w2] = Weights;
    A1 = 2*x0; A2 = 2*x1; A3 = 2*x2;
    B1 = 2*y0; B2 = 2*y1; B3 = 2*y2; #C1 = -1; C2 = -1; C3 = -1;
    D1 = -(x0**2+y0**2-w0); D2 = -(x1**2+y1**2-w1); D3 = -(x2**2+y2**2-w2);
    Det = -A1*B2 - A2*B3 - A3*B1 + A3*B2 + A1*B3 + A2*B1;

    # If no solution is possible, reset the location of sites until there is a solution
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

    # If the intersection point is not within the canvas, reset the location of sites until it is within canvas
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

    #2. find the positive direction of the 3 rays.
    ray_obj_all = [];
    if PointIsWithinTriangle(iPoint, Sites):
        #if the iPoint is within the triangle formed by 3 sites, the positive direction of the ray is the one
        #that moves away from the opposite site.
        for k in range(3):
            site_remain = np.arange(3); site_remain = site_remain[site_remain!=k];
            tag1 = site_remain[0]; tag2 = site_remain[1];
            t_bisector = np.array([2*(y_col[tag2][0]-y_col[tag1][0]),2*(x_col[tag1][0]-x_col[tag2][0])]);
            if np.dot(t_bisector,(Sites[k,:] - iPoint)) > 0:
                t_bisector = - t_bisector;
            ray_obj_all.append(RAY(iPoint, t_bisector,[site_label[tag1],site_label[tag2]],site_label[k]));

    else:
        #if the iPoint is not within the triangle formed by 3 sites, the positive direction is defined in a way
        #that can divide the canvas into 3 polygons and each polygon contains 1 site.
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
    """
    This function is designed for n>=4. This function define the positive direction of each ray
    in a way that is the same as findBisectorPositiveRay.
    """    
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
    """
    This function finds the terminal of each ray object.
    """
    #0. get the label of each ray object
    index_iP = ray_obj.index_iP;
    index_ray = ray_obj.index_ray;     
    intersect_coordinates = [];
    simplex = simplices_prune[index_iP];
    a = (simplices_prune == simplex[0]) + (simplices_prune == simplex[1]) + (simplices_prune == simplex[2]);
    arr = np.nonzero((a.sum(axis = 1) == 2))[0]

    #1. check if the ray is intersected with any nearby iPoints
    for i in arr:
        ipoint = iPoints[i, 0:2];
        if isOnRay(ipoint, ray_obj):
            if len(intersect_coordinates) == 0:
                intersect_coordinates = ipoint;
            else:
                intersect_coordinates = np.vstack([intersect_coordinates, ipoint]); 

    #2. if the ray is not intersected with any iPoints, find its intersection with canvas
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

    #3. only keep the nearest point to the origin of that ray
    if len(intersect_coordinates) > 0:
        intersect_coordinates = intersect_coordinates.reshape((-1,2));
    if len(intersect_coordinates) > 1:
        intersect_coordinates = keepNearestPointsOnRay(intersect_coordinates, ray_obj);

    #4. append origin and create an edge
    if len(intersect_coordinates) == 0:
        intersect_coordinates = ray_obj.origin;
    else:
        intersect_coordinates = np.vstack([intersect_coordinates, ray_obj.origin]);

    #5. round and leave only the unique point
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
    """
    If a ray is intersect with multiple iPoints, only keep the one that is closest to the origin of the ray.
    """
    n_candidates = len(intersect_coordinates);
    # find candidates with shortest distance from the origin;
    dist = np.sum((intersect_coordinates - np.tile(ray_obj.origin,(n_candidates,1)))**2, axis = 1);
    dist = (dist*(10**9)).astype(int)/(10.**9);
    dist[dist == 0] = np.inf; # avoid the intersection point that is the origin itself
    intersect_coordinates = intersect_coordinates[np.argmin(dist), :];
    return intersect_coordinates

def lineIntersectWithEdge(bisector, edge_canvas):
    """
    This function finds the intersection points of a line with an edge.
    """
    [P1, P2] = edge_canvas;
    [[x1, y1],[x2, y2]] = edge_canvas;
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    [D, E, F] = bisector;
    Det = A*E - B*D; 
    DetX = C*E - F*B; 
    DetY = A*F - C*D; 

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
    """
    This function finds the intersection points of a bisector with a canvas.
    """
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
    """
    This function determines which iPoints should be kept. 
    Only the iPoints that (1) belong to the lower convex hull and (2) locate within canvas will be kept.
    """
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
        c = delta/gamma; #z = ax+by+c => nf = (a, b, -1)

        # (1) determine if the face belongs to the lower convex Hull
        COM_face = findCOM(Dual_Sites[simplex])
        if np.dot([a, b, -1], (COM_face - COM_DualSites))>0:

            # (2) determine if the iPoint is within canvas
            if PointIsWithinCanvas(np.array([a/2, b/2]), canvas_obj):
                iPoints = np.append(iPoints, np.array([a/2, b/2, -c]), axis=0)
            else:
                prune_list.append(i);
        else:
            prune_list.append(i);

    iPoints = iPoints.reshape(-1,3);
    simplices_prune = np.delete(simplices_prune, prune_list, axis = 0);
    return iPoints, simplices_prune

def isIntersect(edge1, edge2):
    """
    This function determines if 2 edges are intersected with each other. 
    """
    [P1, P2] = edge1; [P3, P4] = edge2;
    [[x1, y1],[x2, y2]] = edge1
    [[x3, y3],[x4, y4]] = edge2
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    D = (y4-y3); E = -(x4-x3); F = x3*(y4-y3)-y3*(x4-x3);
    Det = A*E - B*D;
    DetX = C*E - B*F;
    DetY = A*F - C*D;
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
    """
    This function reorders the corners of a polygon in a counterclockwise manner.
    The input should be a numpy array of nx2.
    """
    orderedPoints = corners;
    COM = orderedPoints.mean(axis = 0);#find center of mass
    orderedPoints=orderedPoints[np.argsort(np.angle((orderedPoints - COM)[:,0] + (orderedPoints - COM)[:,1]*1j))]
    return orderedPoints

def RandomPointsInCanvas(n_Points, canvas):
    """
    This function assigns n random points within the canvas. The algorithm is as followed:
    (1) divide the canvas into multiple triangels
    (2) weight each triangle according to the area
    (3) randomly placed points within the triangles
    Please make sure that the points have been reordered properly!
    """
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
    """
    Calculate polygon area using shoelace formula.
    Please make sure that the corners are reordered before calling PolygonArea function!
    """
    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def findCOM(corners):
    """
    find the center of mass of the polygon corners.
    """
    COM = corners.mean(axis = 0);
    return COM

def isOnSegment(P, edge):
    """
    Determine whether a point is on a segment by checking if Ax+By-C == 0 & falls between the two
    corners which define the edge.
    """
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
    """
    Determine whether a point is on a ray by checking if Ax+By-C == 0 & its relative position with
    respect to the origin is on the same direction as the tangent vector of the ray.
    """
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
    """
    Check if a point is one of the corner of the canvas.
    """
    a = (canvas_obj.xy == P)
    if len(np.nonzero(a[:,0] & a[:,1])[0]) == 1:
        return True
    else:
        return False

def PointIsWithinCanvas(P, canvas_obj):
    """
    Check if a point is within canvas by computing the sum of the angle formed by adjacent corners
    and the point. If a point is within canvas, the sum of the angle should be 2*pi.
    """
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
    """
    Check if a point is within a triangle. The algorithm used is the same as PointIsWithinCanvas.
    However, this function does not require the input to be a POLYGON object.
    """
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

def FindRayIntersectWithCanvas(ray_obj, canvas_obj):
    """
    This function finc the point where a ray object is intersected with canvas by checking the intersection
    with every edge of the canvas.
    """
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
    """
    This function determines if a ray is intersected with an edge.
    """
    [P1, P2] = edge_canvas;
    [[x1, y1],[x2, y2]] = edge_canvas;
    #convert to Ax+By = C, Dx+Ey = F
    A = (y2-y1); B = -(x2-x1); C = x1*(y2-y1)-y1*(x2-x1);
    D, E, F = ray_obj.a_r, ray_obj.b_r, -ray_obj.c_r;
    Det = A*E - B*D; 
    DetX = C*E - B*F; 
    DetY = A*F - C*D;
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