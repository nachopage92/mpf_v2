def NURBS2D_MESHING(test,npts,**kwargs):

    #------------------------------------------------------------------
    #
    # #  NOTA PERSONAL
    # los puntos de control son parametros que generan la geometria
    # los nodos  viene dado por la distribucion de pts knots [0 .. 1]
    #
    # el refinamiento del vector knot viene acompañado adicion de nuevos
    # pts de control. modifica la funcion de forma NURBS.
    #
    # make triangle mesh triangula los PUNTOS DE CONTROL! no los nodos
    # para ello utilizar rutina 'Delaunay' de scipy.spatial
    #
    #
    # #  ALGUNAS VARIABLES DE INTERES PERTENECIENTES A surf
    # # variables relevantes: numero de puntos de control (u y v)
    # size_u = surf.ctrlpts_size_u
    # size_v = surf.ctrlpts_size_v
    #
    # # orden de Bspline/NURBS
    # degree_u = surf.degree_u
    # degree_v = surf.degree_v
    #
    # # vector 'knot'
    # knotvector_u = surf.knotvector_u
    # knotvector_v = surf.knotvector_v
    #
    # # puntos de control
    # w = surf.weights            # weight
    # ctrlpts = surf.ctrlpts      # unweighted ctrlpts
    # ctrlpts2d = surf.ctrlpts2d  # weighted ctrlpts array
    # ctrlptsw = surf.ctrlptsw    # weighted ctrls
    #
    #------------------------------------------------------------------
	
    import os
    
    import numpy as np
    
    from geomdl import NURBS
    from geomdl import exchange
    from geomdl import tessellate
    from geomdl import compatibility
    from geomdl import utilities
    from geomdl import helpers
    from geomdl import operations
    from geomdl.visualization import VisPlotly
    
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    
    from scipy.spatial import Delaunay,ConvexHull
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
    
    #   INPUT DATA
    
    # numero de puntos en cada direccion u,v (x,y)
    cnptx,cnpty = npts   # numero de puntos en direccion x e y
    filename    = test+"_smesh.dat"     # input file 
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
    
    #   SETTING B-SPLINE/NURBS OBJECT
    
    # Create a BSpline surface instance
    surf = NURBS.Surface()
    # Import smesh file (contains ctrlpts,knotvectors,degree)
    imported = exchange.import_smesh("./ctrlpts_data/"+filename)
    surf = imported[0]

    # ------------------------------------------------------------

    ## parametro rutina VisPlotly
    #delta = kwargs.get('delta',[0.05,0.05])
    #delta_u,delta_v = delta

    #density_refinement = kwargs.get('density_refinement',[0,0])
    #density_u,density_v = density_refinement

    ## Refines the knot vector on the v-direction of a surface
    #operations.refine_knotvector(surf,[density_u,density_v])

    ## Set evaluation delta
    #surf.delta_u = delta_u
    #surf.delta_v = delta_v

    ## Evaluate surface
    #surf.evaluate()

    # Informacion de los contornos
    subs = kwargs.get('subs',[[],[],[],[]])
    tol  = kwargs.get('edge_tol',1.0e-12)

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

    # COLLOCATION & SOURCE PTS IN PARAMETRIC SPACE
    
    # collocation pts
    collocation_uknots = np.linspace(0.0,1.0,cnptx)
    collocation_vknots = np.linspace(0.0,1.0,cnpty)
    ## source pts
    #source_uknots = np.linspace(0.0,1.0,snptx-1)
    #source_vknots = np.linspace(0.0,1.0,snpty-1)
    
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
    
    #   GENERACION DE PUNTOS
    
    # collocation pts
    physical_pts=[] ; parametrical_pts=[] ; index=[] ; c=1
    for vnod in collocation_vknots:
      for unod in collocation_uknots:
        uv_knots = [unod,vnod] # aux var
        index.append(str(c))
        parametrical_pts.append(uv_knots)
        physical_pts.append(\
        mapping_from_parametrical_to_physical_coordinates_2d(uv_knots,surf))
        c += 1
    
    #   MALLADOR en espacio parametrico mediante triangulacion Delaunay

    # triangulacion delaunay
    delaunay_tria = Delaunay(parametrical_pts)
    triangulation = [] 
    for tria in delaunay_tria.simplices:
        aux = []
        for tri in tria:
            aux.append(tri)
        triangulation.append(aux)

    # anillo convexo
    convhull = ConvexHull(parametrical_pts)
    # convhull.simplices : elementos simplex (en 2D son lineas)
    # convhull.vertices  : vertices del anillo convexo (antihorario en 2D)
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

    #   CREAR CONTORNOS (Geometry)

    #       a) PUNTOS DE COLOCACION (NODOS) EN EL CONTORNO

    colloc_l1,colloc_l2 = \
    create_matrix_outline(cnptx,cnpty)

    #       b)  CALCULAR VECTOR NORMAL ASOCIADO A CADA PUNTO,
    #               En los contornos del dominio bidimensional
    #               se 'contrae' una dimension,  resultando en
    #               una curva (tambien bidimensional)

    #               b.1) se determinan los puntos de control del contorno
    ctrlpts_l1,ctrlpts_l2 = \
    create_matrix_outline(surf.ctrlpts_size_u,surf.ctrlpts_size_v)

    types=[] ; newsubs=[]
    for pos in range(4):
        # se extrae puntos de la lista encadenada
        ptslist = colloc_l1[colloc_l2[pos]+1:colloc_l2[pos+1]+1]

        if ( pos == 0 ) : # lado inferior
            colloc_knots = collocation_uknots
        elif ( pos == 1 ):   # lado derecho
            colloc_knots = collocation_vknots
        elif ( pos == 2 ) : # lado superior
            colloc_knots = list(collocation_uknots)
            colloc_knots.reverse()
            colloc_knots = np.array(colloc_knots)
        elif ( pos == 3 ):   # lado izquierdo
            colloc_knots = list(collocation_vknots)
            colloc_knots.reverse()
            colloc_knots = np.array(colloc_knots)

        bdrynum = len(subs[pos])+1 # +1 por cada nuevo punto xi
        newbdry1=[] ; newbdry2=[] ; aux=[]
        cbdry = 0   ; new_elem=0
        types.append(pos) ; newsubs.append(subs[pos])

        # se verifica si existen subcontornos
        for idpt,xipt in zip(ptslist,colloc_knots):
            aux.append(idpt) ; cbdry+=1 
            for jpt in subs[pos]:
                if ( jpt == xipt ) :
                    new_elem +=1
                    newbdry1.append(aux)
                    newbdry2.append(cbdry)
                    aux = [idpt] ; cbdry+=1
                    types.append(pos) ; newsubs.append(subs[pos])
        newbdry1.append(aux)
        newbdry2.append(cbdry)
        new_colloc_l2 = [ colloc_l2[pos] ]

        for i in newbdry2:
            new_colloc_l2.append( colloc_l2[pos]+i )
        new_colloc_l1 = [ ]
        for xp in newbdry1:
            new_colloc_l1.extend( xp )

        # actualizacion de las listas encadenadas
        part1 = colloc_l1[:colloc_l2[pos]+1]
        part1.extend(new_colloc_l1)
        part1.extend(colloc_l1[colloc_l2[pos+1]+1:])
        colloc_l1 = part1

        part2 = colloc_l2[:pos]
        part2.extend(new_colloc_l2)
        aux = []
        for i in colloc_l2[pos+2:]:
            aux.append(i+new_elem*2)
        part2.extend(aux)
        colloc_l2 = part2


    # --------- end for ---------------

    #               b.2) loop sobre cada punto del contorno en sentido antihorario
    contador = -1 ; nsides = len(colloc_l2) ; subs = newsubs ; ind = 1
    normal_list1 = [] ; normal_list2 = [ contador ] 
    for pos in range(nsides-1): # recorre cada lado

        ctrlpts_list = ctrlpts_l1[ctrlpts_l2[types[pos]]+1:ctrlpts_l2[types[pos]+1]+1]
        #       TENER EN CONSIDERACION: LOS VERTICES SE REPITEN! 
        #       lado-1 : v1._____.v2 (v1 contabilizado dos veces)  
        #       lado-2 : v2._____.v3 (v2 contabilizado dos veces)  
        #       lado-3 : v3._____.v4 (v3 contabilizado dos veces)  
        #       lado-4 : v4._____.v1 (v4 contabilizado dos veces)  
        
        # transponer surf.ctrlpts2d
        ctrlpts2d = tranpose_array2d(surf.ctrlpts2d)
        # se obtienen los puntos de control asociados a cada
        # curva del contorno
        ctrlptsw=[]
        for ids in ctrlpts_list:
            ctrlptsw.append(ctrlpts2d[ids])
        # _____________________________________________
        # |          -- FIGURA 1 --                    |
        # | Ejemplo:  'contorno inferior'              |
        # |  __________                                | 
        # | |          |       \   xi=0.0       xi=1.0 | 
        # | | eta      |  ======\  ___________________ |
        # | | |__ xi   |  ======/  xi=variable ; eta=0 | 
        # | |__________|       /                       |        
        # |                                            | 
        # |  Espacio 2D   -- contrae-->   Espacio 1D   |
        # |                                            |   
        # | S(xi=0,eta) = C(eta) : contorno izquierdo  |  
        # | S(xi=1,eta) = C(eta) : contorno derecho    |    
        # | S(xi,eta=0) = C(xi)  : contorno inferior   |     
        # | S(xi,eta=1) = C(xi)  : contorno superior   |      
        # |____________________________________________|

        #  --- algunas modificaciones segun  ---
        if ( types[pos] == 0 ) : # lado inferior
            degree = surf.degree_u
            knots = surf.knotvector_u
            colloc_knots = collocation_uknots
        elif ( types[pos] == 1 ):   # lado derecho
            degree = surf.degree_v
            knots = surf.knotvector_v
            colloc_knots = collocation_vknots
        elif ( types[pos] == 2 ) : # lado superior
            degree = surf.degree_u
            knots = surf.knotvector_u
            colloc_knots = list(collocation_uknots)
            colloc_knots.reverse()
            colloc_knots = np.array(colloc_knots)
        elif ( types[pos] == 3 ):   # lado izquierdo
            degree = surf.degree_v
            knots = surf.knotvector_v
            colloc_knots = list(collocation_vknots)
            colloc_knots.reverse()
            colloc_knots = np.array(colloc_knots)
        else:
            from sys import exit
            exit('ERROR: contorno no identificado')

        # crea elemento GEOMDL
        curve = NURBS.Curve()
        curve.degree = degree
        curve.ctrlptsw = ctrlptsw
        curve.knotvector = knots 
        ptslist = colloc_l1[colloc_l2[pos]+1:colloc_l2[pos+1]+1]

        # se recorre cada punto de colocacion del contorno
        for idpt,xipt in zip(ptslist,colloc_knots):
            xipt_ = xipt
            for jpt in subs[pos]:
                if ( jpt == xipt ) :
                    if ( jpt == 0.0 ) : 
                        xipt_ = xipt + tol 
                    elif ( jpt == 1.0 ) : 
                        xipt_ = xipt - tol
                    else :
                        ind += 1
                        xipt_ = xipt - tol*pow(-1.0,ind)
            normal_list1.append( curve.normal(xipt_,normalize=True) )
            # OBS: normal_list2 = colloc_l2
            contador += 1
            # x =  mapping_from_parametrical_to_physical_coordinates_1d(xipt,curve)
        normal_list2.append( contador )

    # D) Creacion de elementos lineas pertenecientes al contorno
    line_elements = []
    for pos in range(4): # recorre cada lado
        pts_list = colloc_l1[colloc_l2[pos]+1:colloc_l2[pos+1]+1]
        for i in range(len(pts_list)-1):
            elem1 = pts_list[i]
            elem2 = pts_list[i+1]
            line_elements.append( [elem1,elem2] )

    # output variable, contiene las listas enlazadas y las coordenadas
    # y normales asociadas a cada punto
    boundary = [normal_list1,colloc_l1,colloc_l2,line_elements]

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
    
    #   OUTPUT : se entregan informacion relacionada a la geometria

    return physical_pts,triangulation,boundary

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

#   END OF NURBS2D_MESHING

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

#   FUNCIONES DE MAPEO ( parametric coords ---> physics coords )

def mapping_from_parametrical_to_physical_coordinates_2d(uv_knots,surf):
    from geomdl import helpers
    from numpy import array,multiply
    # parametric space knots unpack
    uknot,vknot = uv_knots
    # span
    span_u = helpers.find_span_linear(surf.degree_u,surf.knotvector_u,surf.ctrlpts_size_u,uknot) #tambien se puede utilizar find_span_binsearch()
    span_v = helpers.find_span_linear(surf.degree_v,surf.knotvector_v,surf.ctrlpts_size_v,vknot)
    # basis function
    N = helpers.basis_function(surf.degree_u,surf.knotvector_u,span_u,uknot)
    M = helpers.basis_function(surf.degree_v,surf.knotvector_v,span_v,vknot)
    # mapping to phyisical space
    suma = array([0.0,0.0])
    for nj in range(surf.degree_v+1):
        for ni in range(surf.degree_u+1):
            suma += multiply(\
                surf.ctrlpts2d[span_u-surf.degree_u+ni][span_v-surf.degree_v+nj][0:2],\
                N[ni]*M[nj])
            suma = list(suma)
    # S(xi,eta) = sum_{ij} [ CtrlPts_{ij} * R_{ij} ]
    physical_pts = suma
    return physical_pts

def mapping_from_parametrical_to_physical_coordinates_1d(uknot,curve):
    from geomdl import helpers
    from numpy import array,multiply
    # span
    span_u = helpers.find_span_linear(curve.degree,curve.knotvector,curve.ctrlpts_size,uknot) #tambien se puede utilizar find_span_binsearch()
    # basis function
    N = helpers.basis_function(curve.degree,curve.knotvector,span_u,uknot)
    # mapping to phyisical space
    suma = array([0.0,0.0])
    for ni in range(curve.degree+1):
        suma += multiply(curve.ctrlpts[span_u-curve.degree+ni][0:2],N[ni])
        suma = list(suma)
    # S(xi) = sum_{i} [ CtrlPts_{i} * R_{i} ] unidimensional
    physical_pts = suma
    return physical_pts

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

"""
    create_matrix_outline(nx,ny):
    determina los bordes de la matriz de indices, estos
    corresponden a los indices de los puntos pertenecientes
    al contono. 
        Input
    : nx , numero de nodos en la direccion x
    : ny , numero de nodos en la direccion x
        Output
    : id_list1 , lista enlazada 1 (contiene id puntos) 
    : id_list2 , lista enlazada 2 (contiene posiciones de id_list1)
    
"""
def create_matrix_outline(nx,ny):
    from numpy import array
    c = 0 
    idmat=[] #  crear matriz con los indices globales 
                     # asociados a los puntos de control
    for j in range(ny):
        rows=[]
        for i in range(nx):
            rows.append(c)
            c += 1
        idmat.append(rows)
    idmat = array(idmat)
    # obtener puntos de control asociado a cada contorno
    list1 = idmat[0,:]  # lado inferior dominio / superior matrix 
    list2 = idmat[:,-1] # lado derecho dominio / matrix
    list3 = idmat[-1,:] # lado superior dominio / inferior matrix
    list4 = idmat[:,0]  # lado izquierdo dominio / matrix
    
    # crear listas enlazadas (linked list)    _____________________      
    id_list1 = []   # contiene id de ctrlpts |   -- FIGURA 2 --    |
    id_list2 = [-1] # lista de posicion      |                     |
    c = -1        #                          |       _list4_       |
    # van ordenados de la siguiente forma    |      |       |      |                 
    #  1. inferior  : i=0,nx-1 ; j=0         | list1|       |list2 |         
    #  2. derecho   : i=nx-1   ; j=0,ny-1    |      |_______|      |          
    #  3. superior  : i=nx-1,0 ; j=ny-1      |        list3        |    
    #  4. izquierdo : i=0      ; j=ny-1,0    |_____________________|  
    
    #   segun la Figura 2, los el orde de los lados es l3,l2,l4,l1.
    # notar que l4 y l1 estan invertidas (Sentido antihorario)
    aux = list(list3) ; aux.reverse() ; list3 = array(aux)
    aux = list(list4) ; aux.reverse() ; list4 = array(aux)

    for lists in [list1,list2,list3,list4]:
        for lis in lists:
            id_list1.append(lis)
            c+=1         
        id_list2.append(c)

    return id_list1,id_list2

"""
    tranpose_array2d(input_list): 
    obtiene el arreglo transpuesto (o transpuesta de lista)
    del arregle de entrada 'input_list'
"""
def tranpose_array2d(input_list):
    nrow = len(input_list) 
    ncol = len(input_list[0])
    ctrlpts2d_transpose = [ [ 0 for i in range(nrow) ] for i in range(ncol) ] # C-like loop
    for i in range(nrow):
        for j in range(ncol):
            ctrlpts2d_transpose[j][i] = input_list[i][j]
    ctrlpts2d_flatten = []
    for rows in ctrlpts2d_transpose:
        for elem in rows:
            ctrlpts2d_flatten.append(list(elem))
    return ctrlpts2d_flatten

    
#   EXPORTAR DATOS
    
#export = kwargs.get('export',False)
#
##   EXPORTAR MALLA: modificar a gusto
#if ( export ) : 
#    
#    # si no existe la carpeta la crea
#    dir_ = '../DATOS/test2d/mesh/'
#    if ( not os.path.isdir(dir_) ):
#        os.mkdir(dir_)
#    
#    from numpy import prod
#    filename = '../DATOS/test2d/mesh/'+test+'_'+str(prod(npts))+'.dat'
#    f = open(filename,'w')
#    f.write('Coordinates\n')
#    for pts in physical_source_pts:
#      f.write('{0:20.12e} {1:20.12e}\n'.format(pts[0],pts[1]))
#    f.write('End Coordinates\n')
#    f.write('Elements\n')
#    for elems in triangulation.simplices.copy():
#      f.write('{0:5d} {1:5d} {2:5d}\n'.format(elems[0],elems[1],elems[2]))
#    f.write('End Elements\n')
#    f.close()
#
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

## CREAR OBJETOS FIGURAS
#
#plot = kwargs.get('plot',False)
#if ( plot ) : 
#
#    #   GRAFICAR GEOMETRIA, graficar a gusto
#    package = ['Visplotly','matplotlib']
#    package = package[1]
#    
#    if ( package== 'matplotlib' ):
#
#        fig = plt.figure(figsize=(10,6)) #crea elemento 'figure'
#        ax1 = fig.add_subplot(121)       # subfigure1 (left)
#        ax2 = fig.add_subplot(122)       # subfigure2 (right)
#        fig.suptitle('collocation and sources points',fontsize=24)
#        plt.rc('text', usetex=True) # latex things
#        plt.rc('font', family='serif') # more latex things
#    
#        npts = cnptx*cnpty
#        
#        x,y = zip(*physical_pts)
#        ax1.set_xlabel(r'$x$')
#        ax1.set_ylabel(r'$y$')
#        ax1.scatter(x,y,color='blue',label='Physical Pts. ('+str(npts)+')')
#        #for i in range(len(x)):
#        #  ax1.text(x[i],y[i],source_index[i])
#        #ax1.legend()
#        
#        x,y = zip(*parametrical_pts)
#        ax2.set_xlabel(r'$\xi$')
#        ax2.set_ylabel(r'$\eta$')
#        ax2.scatter(x,y,color='blue')
#        #for i in range(len(xsource)):
#        #  ax2.text(xsource[i],ysource[i],source_index[i])
#        #ax2.legend()
#        
#        # graficar contorno
#        plotting = 100
#        bdry1=[];bdry2=[];bdry3=[];bdry4=[]
#        # borde 1 ( u in [0,1] ; v = 0  )
#        bdry1.append(np.linspace(0.0,1.0,plotting))
#        bdry1.append([ 0.0 for i in range(plotting)])
#        bdry1=list(zip(*bdry1))
#        # borde 2 ( u = 1 ; v in [0,1]  )
#        bdry2.append([ 1.0 for i in range(plotting)])
#        bdry2.append(np.linspace(0.0,1.0,plotting))
#        bdry2=list(zip(*bdry2))
#        # borde 3 ( u in [0,1] ; v = 1.0  )
#        bdry3.append(np.linspace(1.0,0.0,plotting))
#        bdry3.append([ 1.0 for i in range(plotting)])
#        bdry3=list(zip(*bdry3))
#        # borde 4 ( u = 0 ; v in [0,1]  )
#        bdry4.append([ 0.0 for i in range(plotting)])
#        bdry4.append(np.linspace(1.0,0.0,plotting))
#        bdry4=list(zip(*bdry4))
#        # joint bdrys
#        bdry_plot=[]
#        for side in [bdry1,bdry2,bdry3,bdry4]:
#          bdry_plot.extend(side)
#        # mapping to physical space
#        physical_bdry_pts = []
#        for bdry_pts in bdry_plot:
#          physical_bdry_pts.append(mapping_to_physical_space_2d(bdry_pts))
#        
#        x,y = zip(*physical_bdry_pts)
#        ax1.plot(x,y,color='black')
#        x,y = zip(*bdry_plot)
#        ax2.plot(x,y,color='black')
#        
#        plt.savefig('./figura.svg')
#        
#        show = kwargs.get('show',False)
#        if ( show ) :
#            plt.show()
#    
#    if ( package== 'VisPlotly' ):
#    
#        # Plot the control point grid and the evaluated surface
#        vis_comp = VisPlotly.VisSurface()
#        surf.vis = vis_comp
#        surf.render()
#        
#        # Good to have something here to put a breakpoint
#        pass

#::::::::::::::::::::::::::::::::

# PLOT CONVHULL

    # --------------------------------------------------------------------

    #   grafica de convex hull
    #points = np.array(parametrical_pts)
    #fig = plt.figure(figsize=(10,6))
    #plt.plot(points[:,0], points[:,1], 'o')
    #for simplex in convhull.simplices:
    #    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    #plt.show()
