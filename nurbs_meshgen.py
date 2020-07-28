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


    ## ---------- una manera, cargar datos externamente: -----------
    ## Import smesh file (contains ctrlpts,knotvectors,degree)
    #imported = exchange.import_smesh("./ctrlpts_data/"+filename)
    #surf = imported[0]
    ## -------------------------------------------------------------

    # ----------- otra manera es cargarlo directamente-------------
    from init_data import nurbs_geom_data
    degree,dim,knotvector,ctrlptsw = nurbs_geom_data(test)

    surf.degree_u = degree[0]
    surf.degree_v = degree[1]
    surf.ctrlpts_size_u = dim[0]
    surf.ctrlpts_size_v = dim[1]
    surf.ctrlptsw = ctrlptsw
    surf.knotvector_u = knotvector[0]
    surf.knotvector_v = knotvector[1]
    #xctrl,yctrl,zctrl,wctrl = zip(*ctrlptsw)
    ## -------------------------------------------------------------

    w         = surf.weights            # weight
    ctrlpts   = surf.ctrlpts      # unweighted ctrlpts
    ctrlpts2d = surf.ctrlpts2d  # weighted ctrlpts array
    ctrlptsw  = surf.ctrlptsw    # weighted ctrls

    # ------------------------------------------------------------

    ## parametro rutina VisPlotly
    #delta = kwargs.get('delta',[0.05,0.05])
    #delta_u,delta_v = delta

    # Refines the knot vector on the v-direction of a surface
    density_refinement = kwargs.get('refinement',[0,0])
    density_u,density_v = density_refinement
    operations.refine_knotvector(surf,[density_u,density_v])

    ## Set evaluation delta
    #surf.delta_u = delta_u
    #surf.delta_v = delta_v

    ## Evaluate surface
    #surf.evaluate()

    # Informacion de los contornos
    subs = kwargs.get('subs',[[],[],[],[]])
    tol  = kwargs.get('edge_tol',1.0e-12)

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.

    # COLLOCATION PTS IN PARAMETRIC SPACE
    collocation_uknots,collocation_vknots = kwargs.get('input_knots',\
            [np.linspace(0.0,1.0,cnptx),np.linspace(0.0,1.0,cnpty)])
    
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

    nnod  = c-1 # numero de nodos totales
    nnorm = 2*(cnptx-1) + 2*(cnpty-1)
    
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

    # puntos del contorno
    colloc_l1,colloc_l2 = \
    create_matrix_outline(cnptx,cnpty)

    # verifica si existen subcontornos
    types=[] ; newsubs=[] ; cstop=[] ; cl1=[] ; cl2=[-1] ; cont=-1 
    for pos in range(4):
        # se extrae puntos de la lista encadenada
        ptslist = colloc_l1[colloc_l2[pos]+1:colloc_l2[pos+1]+1]

        if ( pos == 0 ) : # lado inferior
            colloc_knots = list(collocation_uknots)
        elif ( pos == 1 ):   # lado derecho
            colloc_knots = collocation_vknots
        elif ( pos == 2 ) : # lado superior
            colloc_knots = list(collocation_uknots)
            colloc_knots.reverse()
        elif ( pos == 3 ):   # lado izquierdo
            colloc_knots = list(collocation_vknots)
            colloc_knots.reverse()

        newbdry1=[] ; newbdry2=[] ; aux=[] ; stops=[]
        stop = -1 ; cbdry = -1 ; new_elem = 0
        types.append(pos) ; newsubs.append(subs[pos])
        for idpt,xipt in zip(ptslist,colloc_knots):
            aux.append(idpt) ; cbdry+=1 ; stop+=1
            for jpt in subs[pos]:
                if ( jpt == xipt ) :
                    new_elem +=1
                    newbdry1.append(aux) ; aux= [idpt]
                    newbdry2.append(cont+1) ; cont+=1
                    types.append(pos)
                    newsubs.append(subs[pos])
                    stops.append(stop)
            cont+=1

        cstop.append(stops)

        newbdry1.append(aux)
        newbdry2.append(cont)

        #new_colloc_l2 = [ colloc_l2[pos] ]
        new_colloc_l2 = [ ]
        for xpos in newbdry2:
            new_colloc_l2.append( xpos )
        new_colloc_l1 = [ ]
        for xpts in newbdry1:
            new_colloc_l1.extend( xpts )

        cl1.extend(new_colloc_l1)
        cl2.extend(new_colloc_l2)

    # new linked list
    colloc_l1 = cl1
    colloc_l2 = cl2

    #  crear elementos lineas
    line_elements = []
    for pos in range( len(colloc_l2)-1 ): # recorre cada lado
        pts_list = colloc_l1[colloc_l2[pos]+1:colloc_l2[pos+1]+1]
        for i in range(len(pts_list)-1):
            elem1 = pts_list[i]
            elem2 = pts_list[i+1]
            line_elements.append( [elem1,elem2] )

    normal,inormal = normales(physical_pts,line_elements,nnorm,nnod)

    # output variable, contiene las listas enlazadas y las coordenadas
    # y normales asociadas a cada punto
    boundary = [colloc_l1,colloc_l2,normal,inormal,line_elements]

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

"""
    normales
    calcula el vector normal de los bordes. se calcula 
    utilizando los elementos linea del contorno
"""
def normales(xcoord,linelem,nnorm,nnod):

    B = [ 0.0 for i in range(nnod*2) ]
    x,y = zip(*xcoord)

    for i in range(nnorm): 

        n1,n2 = linelem[i]

        xn = y[n1]-y[n2]
        yn = x[n2]-x[n1]
        mod = pow(pow(xn,2)+pow(yn,2),0.5)
        xn = xn/mod # coseno directo segun el eje x
        yn = yn/mod # coseno directo segun el eje x

        B[n1]=B[n1]+xn
        B[n1+nnod]=B[n1+nnod]+yn
        B[n2]=B[n2]+xn
        B[n2+nnod]=B[n2+nnod]+yn

    ivect = 0 ; cosdir=[]; inorm=[]
    for i in range(nnod):
        rnormal = pow( pow(B[i],2)+pow(B[i+nnod],2) , 0.5 )
        if (rnormal > .5):
            ivect=ivect+1
            cl=B[i]/rnormal
            cm=B[i+nnod]/rnormal
            inorm.append(i)
            cosdir.append([-cl,-cm]) # outward 
            ivect+=1

    if (ivect!=2*nnorm): # una verificacion
        print('warning!');print('ivect != nnorm')
        print('ivect',ivect);print('nnorm',nnorm)

    return cosdir,inorm


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
