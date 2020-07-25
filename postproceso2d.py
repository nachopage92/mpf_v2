#!/usr/bin/env python3
"""
    Meshfree2d_Posproceso :
Clase de Python3 utilizada para procesar los datos
generados por Test2D.

    Autor: Ignacio Apablaza
    Email: iapablaza.b@gmail.com
"""

class Meshfree2d_Postproceso_Data():

    """
        Inicio de clase Meshfree2d_Posproceso
    """
    def __init__(self, test, npts, E, v, plot, show, formato, **kwargs):
        
        # entradas requeridas
        self.test = test

        # folder names, data
        preproceso_folder       = './DATOS/'
        nubes_folder            = './DATOS/'
        sol_aprox_data_folder   = './DATOS/'
        self.sol_aprox_figure_folder = './GRAFICOS/SOL_APROX/'
        self.sol_exact_figure_folder = './GRAFICOS/SOL_EXACT/'
        self.nubes_figure_folder     = './GRAFICOS/CLD/'
        self.shape_figure_folder     = './GRAFICOS/SHP/'

        from os import path,makedirs
        if not path.exists('./GRAFICOS/'):
            makedirs('./GRAFICOS/')
        
        # relevant filenames, geometry collocations pts
        if ( npts != None ):
            filename = test+'_'+npts
        else:
            filename = test

        self.geometry_filename  = preproceso_folder+filename+'.dat'
        self.elementos_filename = preproceso_folder+filename+'.elem'
        self.cloud_filename     = nubes_folder+filename+'.CLD'
        self.shpfcn_filename    = nubes_folder+filename+'.SF'
        self.solution_filename  = sol_aprox_data_folder+filename+'2d.res'

        # formato de figura exportado
        if ( formato != None ):
            fmt = formato
        else:
            fmt = '.svg' # predeterminado

        # results (approx solution)
        aux1 = self.sol_aprox_figure_folder+filename+'_1'+fmt
        aux2 = self.sol_aprox_figure_folder+filename+'_2'+fmt
        aux3 = self.sol_aprox_figure_folder+filename+'_3'+fmt
        self.sol_aprox_export_filename = [aux1,aux2,aux3]

        # directorio de graficos solucion exacta
        aux1 = self.sol_exact_figure_folder+filename+'_1'+fmt
        aux2 = self.sol_exact_figure_folder+filename+'_2'+fmt
        aux3 = self.sol_exact_figure_folder+filename+'_3'+fmt
        self.sol_exact_export_filename = [aux1,aux2,aux3]

        # cloud figure name
        self.cloud_figure = self.nubes_figure_folder+filename+fmt

        # shape function figure name
        self.shpfcn_figure = self.shape_figure_folder+filename+fmt

        # numero de puntos por lado
        self.npts = npts

        # ------ transforma str en boolean
        self.show = False
        if ( show != None ):
            if ( show.lower() == 'true' ):
                self.show = True
        
        # clasificacion del solver segun el test (customizable)
        if ( E != None ):
            self.young_modulus = float(E)
        if ( v != None ):
            self.poisson_coeff = float(v)
        
        return

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    """
        DRIVER
    permite ejecutar todas las rutinas de postproceso mediante una sola
    funci'on. Las funciones disponibles se encuentran anotadas en el
    diccionario FUNCTION_MAP
    """
    def DRIVER(self,function):
        func = self.Meshfree_Postprocess_Function()
        func = func[function]
        func()
        return

    # SELECTOR DE FUNCIONES (ver rutina DRIVER)
    def Meshfree_Postprocess_Function(self):
        mapping = {'PLOT_SOL_APROX':self.GRAFICAR_SOLUCION_NUMERICA_TEST2D,\
                   'PLOT_SOL_EXACT':self.GRAFICAR_SOLUCION_EXACTA_TEST2D,\
                   'ERROR_GLOBAL':self.ERROR_GLOBAL,\
                   'PLOT_CLOUD':self.GRAFICAR_NUBES}
                   #'ERROR_LOCAL':self.ERROR_LOCAL} # pendiente
                    
        return mapping


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    """
        Shape_Function_Data 
    es un metodo que permite almacenar datos relevantes
    de la funcion de forma, entre ellos:
        phi : funcion de forma
        dphix    : derivada parcial de phi respecto a x
        dphix    : derivada parcial de phi respecto a y
        ddphixx  : segunda derivada de phi respecto a x
        ddphixx  : segunda derivada mixta de phi
        ddphiyy  : segunda derivada de phi respecto a y
        pts_list : lista de puntos de la nube asociado al
            i-esimo nodo estrella
        pts_dist : lista de distacia entre los nodos de la
            nube y el i-esimo nodo estrella
    """
    class Shape_Function_Data():
        def __init__(self,phi,dphix,dphiy,ddphixx,ddphixy,ddphiyy,pts_list,pts_dist):
            self.phi = phi 
            self.dphix = dphix
            self.dphiy = dphiy
            self.ddphixx = ddphixx
            self.ddphixy = ddphixy
            self.ddphiyy = ddphiyy
            self.pts_list = pts_list
            self.pts_dist = pts_dist
            return

    """
        sol_type : formato en que se almacena la solucion
                    displ_i = (u,v)_i
                    tracc_i = (sx,sy,tau)_i            
    """
    class sol_type:
        def __init__(self,displ_sol,tracc_sol):
            displ = []
            for u,v in displ_sol:
                displ.append((u,v))
            tracc = []
            for sx,sy,tau in tracc_sol:
                tracc.append((sx,sy,tau))
            self.displ = displ
            self.tracc = tracc

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    """
        Lectura de elementos
    """
    def LECTURA_ELEMENTOS(self):

        # apertura archivo -> lectura de data
        f = open(self.elementos_filename,'r')
        contents = f.read()
        f.close()
        
        # crea elemento iterable / lista 'file_as_lines'
        file_as_list = contents.splitlines()
        
        # elementos trignauglaes
        ntria = int(file_as_list[1])
        self.tria = []
        for i in range(ntria):
            aux = file_as_list[2+i].split()
            self.tria.append((int(aux[0]),int(aux[1]),int(aux[2])))

        # elementosneas 
        nline = int(file_as_list[3+ntria])
        self.line = []
        for i in range(nline):
            aux = file_as_list[4+ntria+i].split()
            self.line.append((int(aux[0]),int(aux[1])))

        return


    """
        Lectura de geometria
    """
    def LECTURA_DATA(self):

        # apertura archivo -> lectura de data
        f = open(self.geometry_filename,'r')
        contents = f.read()
        f.close()
        
        # crea elemento iterable / lista 'file_as_lines'
        file_as_list = contents.splitlines()
        
        # lectura de cada linea
        geo = {} 
        eof = len(file_as_list)
        
        aux=file_as_list[1].split()
        npoin = int(aux[0])
        dcond = int(aux[1])
        fcond = int(aux[2])
        
        # LECTURA DE COORDENADAS
        coord = []
        for i in range(npoin):
            idpt,x,y = file_as_list[3+i].split()
            idpt = int(idpt) ; x=float(x) ; y = float(y)
            coord.append((x,y))
        coord = tuple(coord)

        # LECTURA DE DESPLAZAMIENTOS PREESCRITOS
        dlist=[] ; cond=[] ; displ=[]
        for i in range(dcond):
            idpt,ifx,ify,u,v = file_as_list[4+npoin+i].split()
            dlist.append(int(idpt)) ; cond.append((int(ifx),int(ify))) ; displ.append((float(u),float(v)))
        dlist=tuple(dlist) ; cond=tuple(cond) ; displ=tuple(displ)

        # LECTURA DE ESFUERZOS PREESCRITOS
        flist=[] ; tracc=[]
        for i in range(fcond):
            idpt,Tx,Ty = file_as_list[5+npoin+dcond+i].split()
            flist.append(int(idpt)) ; tracc.append(tuple([float(Tx),float(Ty)]))
        flist=tuple(flist) ; tracc=tuple(tracc)

        # LECTURA DE VECTORES NORMALES 
        normal = [] ; idnorm=[]
        nnorm = int(file_as_list[6+npoin+dcond+fcond])
        for i in range(nnorm):
            idpt,nx,ny = file_as_list[7+npoin+dcond+fcond+i].split()
            idnorm.append(int(idpt)) ; normal.append((float(nx),float(ny)))
        normal = tuple(normal) ; idnorm = tuple(idnorm)
        
        # se define un objeto que contiene la informacion importada
        class geom_type():
            def __init__(self,npoin,coord,normal,idnorm):
                self.npoin  = npoin
                self.coord  = coord
                self.normal = normal
                self.idnorm = idnorm
        class phys_type():
            def __init__(self,dlist,cond,displ,flist,tracc):
                self.dlist  = dlist
                self.cond   = cond
                self.displ  = displ
                self.flist  = flist
                self.tracc  = tracc

        self.geometry = geom_type(npoin,coord,normal,idnorm)
        self.physical = phys_type(dlist,cond,displ,flist,tracc)
        return


    """
        IMPORTAR_NUBES
    """
    def IMPORTAR_NUBES(self):
    
        f = open(self.cloud_filename,'r')
        contents = f.read()
        f.close()
        file_as_list = contents.splitlines()

        nube = []
        for line in file_as_list:
            aux1 = line.split()
            aux2 = []
            np   = int(aux1[0])
            for i in range(np):
                aux2.append(int(aux1[i+1]))
            nube.append(aux2)

        self.nubes = nube
        return 

    """
        IMPORTAR_FUNCION_DE_FORMA()
        lectura de la funcion de forma (archivo con extension '.SF')
    """
    def IMPORTAR_FUNCION_DE_FORMA(self):
    
        f = open(self.shpfcn_filename,'r')
        contents = f.read()
        f.close()
        file_as_list = contents.splitlines()

        npoin = len(self.nubes)

        shpfcn = [] ; icont = 0
        for i in range(npoin):
            phi = [[] for i in range(6)]
            for j in range(6):
                splitline = file_as_list[icont].split()
                icont +=1
                for line in splitline:
                    phi[j].append(float(line))
            icont +=1 
            shpfcn.append(phi)

        self.shpfcn = shpfcn
        return 


    """
        Lectura de geometria
    """
    def LECTURA_GEOMETRIA(self,filename):
        
        # iniciar variable
        coord =[]

        # apertura archivo -> lectura de data
        f = open(filename,'r')
        contents = f.read()
        f.close()
        
        # crea elemento iterable / lista 'file_as_lines'
        file_as_list = contents.splitlines()
        
        # file_as_list[0] = 'NODOS    D.PRESCRITOS    F.PRESCRITAS'
        aux = file_as_list[1].split()
        nnod=int(aux[0]) ; nfix=int(aux[1]) ; nload=int(aux[2])
        #file_as_list[2] = 'PUNTOS     COORDENADAS'
        for i in range(nnod):
            aux = file_as_list[i+3].split()
            coord.append( [ float(aux[1]) , float(aux[2]) ] )

        return coord

    
    


    """
        GRAFICAR_NUBES:
        Rutina principal que lee las nubes y las grafica
    """
    def GRAFICAR_NUBES(self):
        
        # importar datos de la nube de puntos
        self.IMPORTAR_NUBES()
        self.IMPORTAR_FUNCION_DE_FORMA()
        self.LECTURA_DATA()
        self.LECTURA_ELEMENTOS()

        # mostrar 'n_of_figures' funciones de forma 
        # (son muchos datos para graficar!)
        n_of_figures = 10
        npoin=len(self.geometry.coord)

        #steps=int(npoin/n_of_figures)#        <---------------------
        steps = 1 # grafica todas las func. de forma

        from os import path,makedirs
        folders = [self.nubes_figure_folder,self.shape_figure_folder]
        for folder in folders:
            if not path.exists(folder):
                makedirs(folder)

        for j in range(0,npoin,steps):
            cld = self.nubes[j]
            phi = self.shpfcn[j]
            for i in range(len(cld)): # enumeracion de python
                cld[i] = cld[i] - 1   # parte de cero
            title = r'Nube de proximidad, nodo estrella '+str(j+1)+', test '+self.test
            PLOT_CLOUDS(title,j,self.geometry.coord,cld,self.cloud_figure,self.show)
            PLOT_SHAPES(title,j,self.geometry.coord,cld,phi,self.shpfcn_figure,self.show)
        
        return
    

    
    """
         GRAFICAR_SOLUCION_NUMERICATEST2D
    """
    def GRAFICAR_SOLUCION_NUMERICA_TEST2D(self,**kwargs):

        # lectura de data: 
        self.LECTURA_DATA()
        self.LECTURA_SOLUCION()
        self.LECTURA_ELEMENTOS()

        # titulo
        nx,ny = self.npts.split('-')
        titulo = r'Solucion num\'erica '+self.test+', '+nx+' x '+ny+' nodos'

        # integra elementos latex al grafico (e.g. \omega)
        is_latex=kwargs.get('is_latex',True)

        # triangulacion
        new_tria=[]
        for pts in self.tria:
            new_tria.append(tuple(map(sum, zip(pts, (-1,-1,-1)))))
        new_line=[]
        for pts in self.line:
            new_line.append(tuple(map(sum, zip(pts, (-1,-1)))))

        # archivo a exportar
        from os import path,makedirs
        if not path.exists(self.sol_aprox_figure_folder):
            makedirs(self.sol_aprox_figure_folder)

        ## plotear ALGUNOS puntos
        #       ------              EN DESARROLLO!!!
        #from math import ceil # redondeo hacia arriba
        #nxlim = 31 ; nylim = 31
        #nxint,nyint = [ int(nx) , int(ny) ] ; short_list=[]
        #if ( nxint > nxlim or nyint < nylim ):
        #    stepx,stepy = [ int(nxint/nxlim) , int(nyint/nylim) ]
        #    restx,resty = [ nxint%nxlim , nyint%nylim ]
        #    jumpx,jumpy = [ ceil(nxlim/restx) , ceil(nylim/resty) ]

        PLOT_DISPLACEMENT(titulo,self.geometry.coord,self.solution.displ,\
            new_tria,self.sol_aprox_export_filename[0],is_latex,self.show)

        #PLOT_DISPLACEMENT_VECTOR(titulo,self.geometry.coord,self.solution.displ,\
        #    new_line,new_tria,self.sol_aprox_export_filename[2],is_latex,\
        #    self.show)

        PLOT_TRACTION(titulo,self.geometry.coord,self.solution.tracc,\
            new_tria,self.sol_aprox_export_filename[1],is_latex,self.show)

        return


    """
         GRAFICAR_SOLUCION_TEST2D_EXACTA
    """
    def GRAFICAR_SOLUCION_EXACTA_TEST2D(self,**kwargs):

        # numero de puntos por lado
        npts=[20,20]

        # lectura de data
        from nurbs_meshgen import NURBS2D_MESHING
        from exactsol import SOLUCION_EXACTA
        physical_pts,triangulation,boundary = NURBS2D_MESHING(self.test,npts)
        solucion = []
        for x,y in physical_pts:
            solucion.append( SOLUCION_EXACTA(self.test,x,y,\
                E=self.young_modulus,v=self.poisson_coeff) )

        u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = zip(*solucion)

        displacement = []
        for aux in zip(u,v):
            displacement.append(aux)
        traction = []
        for aux in zip(sigmax,sigmay,tauxy):
            traction.append(aux)

        bdry_lines = boundary[-1]

        # titulo 
        titulo = r'Solucion exacta '+self.test

        # integra elementos latex al grafico (e.g. \omega)
        is_latex=True

        # archivo a exportar
        from os import path,makedirs
        if not path.exists(self.sol_exact_figure_folder):
            makedirs(self.sol_exact_figure_folder)


        PLOT_DISPLACEMENT(titulo,physical_pts,displacement,\
            triangulation,self.sol_exact_export_filename[0],is_latex,self.show)

        #PLOT_DISPLACEMENT_VECTOR(titulo,physical_pts,displacement,\
        #    bdry_lines,triangulation,self.sol_exact_export_filename[2],\
        #        is_latex,self.show)

        PLOT_TRACTION(titulo,physical_pts,traction,\
            triangulation,self.sol_exact_export_filename[1],is_latex,self.show)

        return
    
    
    
    """
         LECTURA_SOLUCION 
    """
    def LECTURA_SOLUCION(self):
      
        # apertura archivo -> lectura de data
        f = open(self.solution_filename,'r')
        contents = f.read()
        f.close()
        
        # crea elemento iterable / lista 'file_as_lines'
        file_as_lines = contents.splitlines()

        npoin = self.geometry.npoin

        sol_list1=[] ; displ_sol=[]
        for i in range(npoin):
            idpt,u,v = file_as_lines[i+1].split()
            idpt=int(idpt) ; u=float(u) ; v=float(v)
            displ_sol.append( (u,v) )
            sol_list1.append( idpt )

        sol_list2=[] ; tracc_sol=[]
        for i in range(npoin):
            aux = file_as_lines[i+2+npoin].split()
            idpt,sx,sy,tau = file_as_lines[i+2+npoin].split()
            idpt=int(idpt) ; sx=float(sx) ; sy=float(sy) ; tau=float(tau)
            tracc_sol.append( (sx,sy,tau) )
            sol_list2.append( idpt )

        
        if (sol_list1==sol_list2):
            self.solution = self.sol_type(displ_sol,tracc_sol)
        else:
            print('LECTURA_SOLUCION: Error')
            print(' sol_list1 =! sol_list2 ')
            from sys import exit ; exit('Program stopped')
            # NOTA: Esta situacion no deberia pasar!!

        return

    """
        ERROR_GLOBAL
    """
    def ERROR_GLOBAL(self):
     
        def calculo_solucion_exacta():
            displ = [] ; tracc = []
            from exactsol import SOLUCION_EXACTA
            for x,y in zip(*self.geometry.coord):
                for u,v,ex,ey,exy,sx,sy,tau in SOLUCION_EXACTA(self.test,x,y,E=E,v=v):
                    displ.append((u,v))
                    tracc.append((sx,sy,tau))

            self.solexact = sol_type(displ,tracc)
            return

        def distancia_internodal():
            # la distancia internodal esta dada por la ecuacion:
            #    h_car = sqrt( Area ) / ( sqrt(npoints) - 1 )
            # (An Introduction to Meshfree Methods and Their Programming, por Liu)
            x,y = zip(*self.geometry.coord)
            Lx = max(x)-min(x) ; Ly = max(y)-min(y)
            # OBS: supone area rectangular: A = Lx * Ly
            dc = pow(Lx*Ly,0.5)/(pow(npoin,0.5)-1)
            return dc
     
     #=======================================================
        # idealmente verificar que los datos sol_aprox y sol_exact 
        # sean coherentes, verificar que se cumpla: 
        # len(sol_exact) = len(sol_aprox)

        # ¿es mejor calcularlo o importarlo?
        self.LECTURA_DATA()
        E = self.young_modulus
        v = self.poisson_coeff
        npoin = len(self.geometry.coord)

        # soluciones a comparar
        calculo_solucion_exacta()    # almacenado en self.solexact
        self.LECTURA_SOLUCION()      # almacenado en self.solucion
        
        # calculo del error
        error = {}
        size = len(self.solexact.displ) # o len(self.solexact.tracc), da igual

        if (self.measure_type=='promedio'):
            for key in sol_exact.keys():
                diff = 0.0
                for i in range(size):
                    diff += abs( sol_exact[key][i] - sol_app[key][i] )
                error[key] = diff/size

        elif (self.measure_type=='promedio-cuadratico'):
            for key in sol_exact.keys():
                diff = 0.0
                for i in range(size):
                    diff += pow( sol_exact[key][i] - sol_app[key][i] , 2 )
                error[key] = pow(diff,0.5)/size

        elif (self.measure_type=='inf-norm'):
            for key in sol_exact.keys():
                diff = 0.0
                for i in range(len(sol_exact[key])):
                    diff += max( abs(sol_exact[key][i]-sol_app[key][i]) , diff )
                error[key] = diff

        error['dc'] = distancia_internodal()
        # pendiente: ingresar norma L2,H2,L2_adaptada,H2_adaptada
        self.error = error
        return


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


def PLOT_SHAPES(title,idpt,coords,connect,phi,filename,show):

    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3

    # coordenadas pts nubes
    cld_coords = []
    for pt in connect:
        cld_coords.append(coords[pt])
        
    # extrae data en vectores
    x,y = zip(*coords)
    xcld,ycld=zip(*cld_coords)
    p,dpx,dpy,ddpxx,ddpxy,ddpyy = phi

    # permite utilizar latex
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # crear figura NUBE
    fig=plt.figure(num=None,figsize=(8,6),dpi=80,facecolor='w',edgecolor='k')
    title += r", idpt="+str(idpt+1)
    fig.suptitle(title, fontsize=12)

    #crear subfigura NUBE
    ax1 = fig.add_subplot(231,projection='3d')
    ax2 = fig.add_subplot(232,projection='3d')
    ax3 = fig.add_subplot(233,projection='3d')
    ax4 = fig.add_subplot(234,projection='3d')
    ax5 = fig.add_subplot(235,projection='3d')
    ax6 = fig.add_subplot(236,projection='3d')

    # asignar ejes
    ax1.set_xlabel(r'$x$');ax1.set_ylabel(r'$y$');ax1.set_zlabel(r'$\phi$')
    ax2.set_xlabel(r'$x$');ax2.set_ylabel(r'$y$');ax2.set_zlabel(r'$\frac{\partial \phi}{\partial x}$')
    ax3.set_xlabel(r'$x$');ax3.set_ylabel(r'$y$');ax3.set_zlabel(r'$\frac{\partial \phi}{\partial y}$')
    ax4.set_xlabel(r'$x$');ax4.set_ylabel(r'$y$');ax4.set_zlabel(r'$\frac{\partial^2 \phi}{\partial x^2}$')
    ax5.set_xlabel(r'$x$');ax5.set_ylabel(r'$y$');ax5.set_zlabel(r'$\frac{\partial^2 \phi}{\partial x \partial y}$')
    ax6.set_xlabel(r'$x$');ax6.set_ylabel(r'$y$');ax6.set_zlabel(r'$\frac{\partial^2 \phi}{\partial y^2}$')

    # plot
    ax1.plot_trisurf(xcld,ycld,p)
    ax2.plot_trisurf(xcld,ycld,dpx)
    ax3.plot_trisurf(xcld,ycld,dpy)
    ax4.plot_trisurf(xcld,ycld,ddpxx)
    ax5.plot_trisurf(xcld,ycld,ddpxy)
    ax6.plot_trisurf(xcld,ycld,ddpyy)
    
    # distancia entre graficos
    fig.tight_layout()
    
    # archivo a exportar
    fmt = filename.split('.')[-1] 
    filename=filename.replace('.'+fmt,'_'+str(idpt+1)+'.'+fmt)

    fig.savefig(filename)
    
    if (show):
        plt.show()
    
    # limpiar variables matplotlib
    plt.cla()    # clear axis
    plt.clf()    # clear figure
    plt.close()  # close a figure windows

    return



def PLOT_CLOUDS(title,idpt,coords,connect,filename,show):
        
    import matplotlib.pyplot as plt

    # coordenadas pts nubes
    cld_coords = []
    for pt in connect:
        cld_coords.append(coords[pt])
        
    # extrae data en vectores
    x,y = zip(*coords)
    xcld,ycld=zip(*cld_coords)
    dx = max(x)-min(x) #set aspect
    dy = max(y)-min(y) # dy/dx

    # permite utilizar latex
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # crear figura NUBE
    fig=plt.figure(num=None,figsize=(8,6),dpi=80,facecolor='w',edgecolor='k')
    title += r", idpt="+str(idpt+1)
    fig.suptitle(title, fontsize=12)

    #crear subfigura NUBE
    ax = fig.add_subplot(111)

    # asignar ejes
    ax.set_xlabel( r'$x$' ) ; ax.set_ylabel( r'$y$' ) 
    ax.set_aspect(dy/dx)

    # plot
    ax.scatter(x,y,color='black')              # puntos del dominio
    ax.scatter(xcld,ycld,color='green',s=80)     # nube de puntos
    ax.scatter(xcld[0],ycld[0],color='red',s=80)     # punto estrella

    # obtener anillo convexo de la nube
    from scipy.spatial import ConvexHull
    cnvx = ConvexHull(cld_coords)
    cnvx_coords=[]
    for pt in cnvx.vertices:
        cnvx_coords.append(cld_coords[pt])
    cnvx_coords.append(cld_coords[cnvx.vertices[0]])
    xcnvx,ycnvx=zip(*cnvx_coords)
    ax.plot(xcnvx,ycnvx)

    # distancia entre graficos
    fig.tight_layout()
    
    # archivo a exportar
    fmt = filename.split('.')[-1] 
    filename=filename.replace('.'+fmt,'_'+str(idpt+1)+'.'+fmt)

    fig.savefig(filename)
    
    if (show):
        plt.show()
    
    # limpiar variables matplotlib
    plt.cla()    # clear axis
    plt.clf()    # clear figure
    plt.close()  # close a figure windows

    return




def PLOT_DISPLACEMENT(titulo,colloc_pts,displ_sol,triang,filename,is_latex,show):

    # -----------  preambulo graficos -------------
    import matplotlib.pyplot as plt # paquete MatPlotLib
    import matplotlib.tri as mtri   # rutinas de triangulacion
    from mpl_toolkits.mplot3d import Axes3D  # ejes 3d
    from matplotlib import cm       # import color map
    from numpy import array

    # formato fuente (permite utilizar latex, e.g. r'\omega')
    if (is_latex):
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    # colormap
    cmap = cm.jet # jet es un mapa de color especifico, se puede cambiar
    
    # desempaquetar data
    x,y=list(zip(*colloc_pts)) 
    u,v = zip(*displ_sol)

    # da formato a la triangulacion para ser utilizada por tricontourf
    triangulacion =  mtri.Triangulation(x, y, triang)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # crea figura Matplotlib
    fig = plt.figure(figsize=(10,6))
    ax1 = fig.add_subplot(221) # u      
    ax2 = fig.add_subplot(222) # v      
    ax3 = fig.add_subplot(223,projection='3d') # u      
    ax4 = fig.add_subplot(224,projection='3d') # v      

    # genera titulo
    fig.suptitle(titulo,fontsize=20)
    
    ## ------------ set limit  ---------------------
    #tol = 1e-6

    #zmin = min(u) ; zmax = max(u) ; ucoef=0.15
    #if ( zmax-zmin < tol ):
    #    if ( zmax-zmin == 0.0 ):
    #        print('no se puede generar grafico de desplazamientos')
    #        print('Solucion identicamente cero ; u(i)=0 para todo i')
    #        return
    #    ucoef = pow(zmax-zmin,-1.0)
    #offset = (zmax-zmin)*ucoef
    #ax3.set_zlim([zmin-offset,zmax+offset])

    #zmin = min(v) ; zmax = max(v) ; vcoef=0.15
    #if ( zmax-zmin < tol ):
    #    if ( zmax-zmin == 0.0 ):
    #        print('no se puede generar grafico de desplazamientos')
    #        print('Solucion identicamente cero ; v(i)=0 para todo i')
    #        return
    #    vcoef = pow(zmax-zmin,-1.0)
    #offset = (zmax-zmin)*vcoef
    #ax4.set_zlim([zmin-offset,zmax+offset])
    ## ----------------------------------------------

    # 2D colourmap plot

    # grafico 1
    ax1.set_title(r'Desplazamiento $u$')
    ax1.set_ylabel(r'$y$')
    ax1.set_xlabel(r'$x$')
    element1 = ax1.tricontourf(triangulacion, u, cmap=cmap)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.colorbar(element1,ax=ax1)
    # grafico 2
    ax2.set_title(r'Desplazamiento $v$')
    ax2.set_ylabel(r'$y$')
    ax2.set_xlabel(r'$x$')
    element2 = ax2.tricontourf(triangulacion, v, cmap=cmap)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.colorbar(element2,ax=ax2)

    # 3D plots

    # grafico 3
    ax3.set_ylabel(r'$y$')
    ax3.set_xlabel(r'$x$')
    ax3.plot_trisurf(triangulacion,u)
    ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
    # grafico 4
    ax4.set_ylabel(r'$y$')
    ax4.set_xlabel(r'$x$')
    ax4.plot_trisurf(triangulacion,v)
    ax4.xaxis.set_major_locator(plt.MaxNLocator(3))

    # exportar grafico como archivo
    fig.savefig(filename)
    
    # muestra graficos del terminal
    if (show):
        plt.show()
    # limpiar variables matplotlib
    plt.cla()    # clear axis
    plt.clf()    # clear figure
    plt.close()  # close a figure windows
    
    return 


def PLOT_DISPLACEMENT_VECTOR(titulo,colloc_pts,displ_sol,bdry_elem,tria,filename,is_latex,show):

    # -----------  preambulo graficos -------------
    import matplotlib.pyplot as plt # paquete MatPlotLib
    from mpl_toolkits.mplot3d import Axes3D  # ejes 3d
    from matplotlib import cm       # import color map
    from numpy import array

    # crea figura Matplotlib
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111) # u      

    # formato fuente (permite utilizar latex, e.g. r'\omega')
    if (is_latex):
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    # colormap
    cmap = cm.jet # jet es un mapa de color especifico, se puede cambiar
    
    # desempaquetar data
    x,y=list(zip(*colloc_pts)) 
    u,v = zip(*displ_sol)
    scale = 1.0

    # genera titulo
    fig.suptitle(titulo+', scale = '+str(scale),fontsize=20)

    x_displ = [] ; y_displ = []
    for xpt,ypt,usol,vsol in zip(x,y,u,v) :
        x_displ.append( (xpt+scale*usol,ypt+scale*vsol) )
    x_displ=tuple(x_displ)

    xd,yd = zip(*x_displ)

    # graficar poligonos
    from matplotlib.patches import Circle, Wedge, Polygon
    from matplotlib.collections import PatchCollection

    bdry_coord=[]
    dspl_bdry_coord=[]
    pair1 = [] ;  pair2 = []
    for pts in bdry_elem:
        p1 = pts[0] ; p2 = pts[1]
        pair1.append( (colloc_pts[p1],colloc_pts[p2]) )
        pair2.append( tuple([tuple([sum(x) for x in zip(colloc_pts[p1], [scale*u[p1],scale*v[p1]])]),\
                              tuple([sum(x) for x in zip(colloc_pts[p2], [scale*u[p2],scale*v[p2]])])] ) )
        bdry_coord.append(pair1[-1][0])
        dspl_bdry_coord.append(pair2[-1][0])

    # CONVEX HULL (contorno del dominio)
    patches = []
    patches.append( Polygon(bdry_coord,True) ) 
    patches.append( Polygon(dspl_bdry_coord,True) ) 
    p = PatchCollection(patches, alpha=0.4, color=['Gray','Green'])
    ax.add_collection(p)

    # puntos de colocacion
    for xpt1,xpt2 in zip(colloc_pts,x_displ):
        x1,y1 = xpt1 ; x2,y2 = xpt2
        ax.scatter(x1,y1,color='Gray')
        ax.scatter(x2,y2,color='Green')

    # triangulacion
    from matplotlib.tri import Triangulation 
    triang_1 = Triangulation(x,y,triangles=tria)
    triang_2 = Triangulation(xd,yd,triangles=tria)
    ax.triplot(triang_1, marker="o",color='Gray')
    ax.triplot(triang_2, marker="o",color='Green')
    
    ## vector desplazamiento
    #ax.quiver(x,y,u,v,scale=1,scale_units='xy',color='Gray')

    # exportar grafico como archivo
    fig.savefig(filename)
    
    # muestra graficos del terminal
    if (show):
        plt.show()
    # limpiar variables matplotlib
    plt.cla()    # clear axis
    plt.clf()    # clear figure
    plt.close()  # close a figure windows
    
    return 



def PLOT_TRACTION(titulo,colloc_pts,tracc_sol,triang,filename,is_latex,show):

    # -----------  preambulo graficos -------------
    import matplotlib.tri as mtri   # rutinas de triangulacion
    import matplotlib.pyplot as plt # paquete MatPlotLib
    from mpl_toolkits.mplot3d import Axes3D  # ejes 3d
    from matplotlib import cm       # import color map
    from numpy import array
    
    # formato fuente (permite utilizar latex, e.g. r'\omega')
    if (is_latex):
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    
    # colormap
    cmap = cm.jet # jet es un mapa de color especifico, se puede cambiar
    
    # desempaquetar data
    x,y=list(zip(*colloc_pts)) 
    sigmax,sigmay,tauxy = zip(*tracc_sol)

    # da formato a la triangulacion para ser utilizada por tricontourf
    triangulacion =  mtri.Triangulation(x, y, triang)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # crea figura Matplotlib
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(231) # sigma x
    ax2 = fig.add_subplot(232) # sigma y
    ax3 = fig.add_subplot(233) # tau xy
    ax4 = fig.add_subplot(234,projection='3d') # sigma x
    ax5 = fig.add_subplot(235,projection='3d') # sigma y
    ax6 = fig.add_subplot(236,projection='3d') # tau xy

    # genera titulo
    fig.suptitle(titulo,fontsize=20)

    ## ------------ set limit  ---------------------
    #threshold = 1e-2

    #zmin = min(sigmax) ; zmax = max(sigmax) ; dif=zmax-zmin ; sxcoef=0.15
    #if ( dif < threshold ):
    #    if ( dif == 0.0 ):
    #        print('no se puede generar grafico de desplazamientos')
    #        print('Solucion identicamente cero ; sigma_x(i)=0 para todo i')
    #        return
    #    sxcoef = pow(dif,-1.0)
    #offset = dif*sxcoef
    #ax4.set_zlim([zmin-offset,zmax+offset])

    #zmin = min(sigmay) ; zmax = max(sigmay) ; dif=zmax-zmin ; sycoef=0.15
    #if ( dif < threshold ):
    #    if ( dif == 0.0 ):
    #        print('no se puede generar grafico de desplazamientos')
    #        print('Solucion identicamente cero ; sigma_y(i)=0 para todo i')
    #        return
    #    sycoef = pow(dif,-1.0)
    #offset = dif*sycoef
    #ax5.set_zlim([zmin-offset,zmax+offset])

    #zmin = min(tauxy) ; zmax = max(tauxy) ; dif=zmax-zmin ; taucoef=0.15
    #if ( dif < threshold ):
    #    if ( dif == 0.0 ):
    #        print('no se puede generar grafico de desplazamientos')
    #        print('Solucion identicamente cero ; tau_xy(i)=0 para todo i')
    #        return
    #    taucoef = pow(dif,-1.0)
    #offset = dif*taucoef
    #ax6.set_zlim([zmin-offset,zmax+offset])
    ## ----------------------------------------------

    ## genera titulo
    #fig.suptitle(titulo,fontsize=20)
    
    # grafico 1
    ax1.set_title(r'Esfuerzos $\sigma_{xx}$')
    ax1.set_xlabel(r'$x$');ax1.set_ylabel(r'$y$')
    element1 = ax1.tricontourf(triangulacion,sigmax,cmap=cmap)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.colorbar(element1,ax=ax1)
    # grafico 2
    ax2.set_title(r'Esfuerzos $\sigma_{yy}$')
    ax2.set_xlabel(r'$x$');ax2.set_ylabel(r'$y$')
    element2 = ax2.tricontourf(triangulacion,sigmay,cmap=cmap)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.colorbar(element2,ax=ax2)
    # grafico 3
    ax3.set_title(r'Esfuerzos $\tau_{xy}$')
    ax3.set_xlabel(r'$x$');ax3.set_ylabel(r'$y$')
    element3 = ax3.tricontourf(triangulacion,tauxy,cmap=cmap)
    ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.colorbar(element3,ax=ax3)       
    # 3D plots

    # grafico 4
    ax4.set_ylabel(r'$y$')
    ax4.plot_trisurf(triangulacion,sigmax)
    ax4.xaxis.set_major_locator(plt.MaxNLocator(3))
    # grafico 5
    ax5.set_ylabel(r'$y$')
    ax5.plot_trisurf(triangulacion,sigmay)
    ax5.xaxis.set_major_locator(plt.MaxNLocator(3))
    # grafico 6
    ax6.set_ylabel(r'$y$')
    ax6.plot_trisurf(triangulacion,tauxy)
    ax6.xaxis.set_major_locator(plt.MaxNLocator(3))

    # exportar grafico como archivo
    fig.savefig(filename)
    
    # muestra graficos del terminal
    if (show):
        plt.show()
    # limpiar variables matplotlib
    plt.cla()    # clear axis
    plt.clf()    # clear figure
    plt.close()  # close a figure windows
    
    return 

#::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::

# ...........  PROGRAMA PRINCIPAL ............. 

#::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::

import argparse

# parse arguments from bash script
parser = argparse.ArgumentParser()
parser.add_argument('-t','--test', required=True,help="(srt) nombre del test. Ejemplo : 'patchtest-0'")
parser.add_argument('-f','--func', required=True,help="(str) funcion a ejecutar: 'PLOT_SOL_APROX'/'PLOT_SOL_EXACT'/'PLOT_CLOUD'")
parser.add_argument('-n','--npts',help="(str) numero de puntos por lado (nx-ny). Ejemplo: '5-5'")
parser.add_argument('-E','--young',help='(str) modulo de elasticidad')
parser.add_argument('-v','--poisson',help='(str) coeficiente de Poisson')
parser.add_argument('-p','--plot',help="(str) exporta graficos en ./GRAFICOS ('True'/'False')")
parser.add_argument('-s','--show',help="(str) muestra graficos en terminal ('True'/'False')")
parser.add_argument('-fmt','--format',help="(str) formato de salida figuras generadas. Ejemplo : '.pdf'")

# lectura argumentos
args = parser.parse_args()

#  Utilizar clase Mesfree2d_Postproceso_Data que contienen
# todas las rutinas (almacenadas como metodos en python3).
# Las funciones habilitadas se encuentran contenidas en el
#  diccionario FUNCTION_MAP 
DATA = Meshfree2d_Postproceso_Data(args.test,args.npts,\
    args.young,args.poisson,args.plot,args.show,args.format)

# call function
DATA.DRIVER(args.func)
