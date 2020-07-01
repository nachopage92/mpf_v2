#!/usr/bin/env python3
"""                                                                
    ... PREPROCESADOR PARA SOLVER MPF 2D ...

     # test disponibles:                                            
     tests = ['patchtest-0','patchest-1','patchtest-2','infinite-plate']

     Requiere los siguientes paquetes:  
         - Numpy        :   useful math package
         - Scipy        :   more useful math package
         - Geomdl       :   B-Spline/NURBS package
         - Matplotlib   :   for plotting
                                                                    
    Para ejecutar:

        > run.bat preproceso (Windows)
        > run.sh preproceso (Linux)

    Para obtener informacion acerca de como utilizar este programa
    escribir en terminal:

        > python3 preproceso.py --help
    
"""                                                              

class Preprocess_Data_Type():
    """
    Clase de python3  que contiene  las rutinas de  preproceso
    necesarias para ejecutar los test bidimensionales. De este
    se obtiene la informacion de la geometria y  del sistema a
    resolver (rhs)
    
    Creacion de elemento de clase 'Preprocess_Data_Type':
    Preprocess_Object = \
        Preprocess_Data_Type(test,nptside,meshgen,show,filename)
    donde:
        - test (str)    : nombre del test
        - nptside(str)  : numero de puntos por dirección,
                          ('ncollocx_ncollocy_nsourcex_nsourcey')
        - show(str)     : indica si se mostrara figuras en terminal
        - filename(str) : prefijo nombre archivo a exportar
    """

    def __init__(self,test,nptside,show,plot,E,v,**kwargs):

        # input data
        self.test = test
        self.string_npts = nptside
        aux = nptside.split('-') 
        self.nptside = [ int(aux[0]) , int(aux[1]) ]  # numero de puntos por lado
        self.show    = show     # mostrar graficos en terminal interactiva
        self.plot    = plot     # exporta figuras en ../GRAFICOS/test2d/
        self.formato = kwargs.get('formato','.svg')
        
        self.young_modulus = float(E)
        self.poisson_coeff = float(v)

        # directorios de salida
        self.data_dir = "./DATOS/"
        self.plot_dir = "./GRAFICOS/"

    """
        DISCRETIZACION DEL DOMINIO
    """
    def GENERAR_GEOMETRIA_DOMINIO(self):
        #       discretizacion regular utilizando algoritmo simple aplicado
        #       al dominio parametrico NURBS. 
        from nurbs_meshgen import NURBS2D_MESHING
        self.coord,self.tri,self.bdry = NURBS2D_MESHING(self.test,self.nptside,plot=False)
        # rescato los elementos del contorno en una variable aparte
        # OBS: self.bdry se actualiza mas adelante!
        self.bdry_elem = self.bdry[-1]

        return

    """
    PREESCRIPCION_CONDICIONES_CONTORNO
        establece las condiciones de contorno asociadas a cada borde
        cada test posee condiciones de contorno propias 
    """
    def PREESCRIPCION_CONDICIONES_CONTORNO(self):
        import numpy as np
        
        # revisar 'order_boundary_list' para mas informacion sobre
        # esta funcion y sobre el significado de 'order_list'
        condition_order_list=self.order_boundary_list()

        # 1. Desempaquetar informacion de variable 'self.bdry'
        bdry_coord,bdry_normal = zip(*self.bdry[0])
        colloc_l1 = self.bdry[1]
        colloc_l2 = self.bdry[2]

        # 2. Crear areglo de permutacion 
        #       2.1 eliminar elementos repetidos, 'bdry_pts'
        #           almacena estos puntos (output var)
        bdry_pts = [] ; vrtx_pts=[]
        for elem in colloc_l1:
            if ( not elem in bdry_pts ):
                bdry_pts.append(elem)
            else:
                vrtx_pts.append(elem)
        #       2.2 crear arreglo 'bdry_normal' y 'condition'
        normal_outward = [ (0.0,0.0,0.0) for i in range(len(bdry_pts)) ]
        bdry_condition = [ 0 for i in range(len(bdry_pts)) ]
        #       2.3 a partir de bdry_pts se crea perm (bdry_pts, NO colloc_l1)
        #           perm[ id_global ] = id_bdry
        perm = {}
        c = 0 
        for elem in bdry_pts:
            perm[elem] = c
            c+=1
        
        aux= []
        for elem in bdry_pts:
            aux.append((self.coord[elem][0],self.coord[elem][1],0.0))

        # replace bdry_coord with the new ordering
        bdry_coord = aux

        # 3.   Loop sobre los puntos de contorno,  en este paso
        #    se asignan las condiciones de contorno y el vector
        #    normal unitario. las condiciones de contorno estan
        #    clasificadas de la siguiente manera.
        #
        #      'id number'  |  'physical condition'   | 'output format'
        #     --------------|-------------------------|-----------------
        #           1       |   only_u_displacement   |    ( 1 , 0 )
        #           2       |   only_v_displacement   |    ( 0 , 1 )
        #           3       |   both_uv_displacement  |    ( 1 , 1 )
        #           4       |   tau_traction          |       ---
        ##
        #   output varible: 
        #       - 'self.bdry_pts'     
        #       - 'self.bdry_normal'  
        #       - 'self.bdry_cond'    


        # Tener en consideracion: la asignacion de las condiciones (a)
        # y  vector normal (b) de cada  punto se  realiza por separado

        #   a). Asignacion de vector normal
        #       b.1) crea diccionario 'vrtx_normal'
        vrtx_normal={}
        for i in vrtx_pts:
            vrtx_normal[i] = []
        for i in range(len(colloc_l1)):
            pos = colloc_l1[i]
            if ( pos in vrtx_pts ):
                vrtx_normal[pos].append(bdry_normal[i])

        #   b) la normal en los vertices se calcula como la suma de las
        #       normales asociadas a las caras de ambos vertices. vector
        #       obtenido es normalizado 
        for vrtx in vrtx_normal:
            normales = vrtx_normal[vrtx]
            suma = [ 0.0 , 0.0 , 0.0 ]
            for nrml in vrtx_normal[vrtx]:
                suma = [ sum(i) for i in zip(suma,nrml) ]
            aux = np.array( suma ) # transforma list in np.array
            aux = aux/np.sqrt(np.dot(aux,aux)) #normalizacion vector
            vrtx_normal[vrtx] = tuple(aux)
        
        #   c) finalmente se asigna la normal a la variable de salida
        for i in colloc_l1:
            if (i in vrtx_normal.keys()):
                normal_outward[perm[i]] = vrtx_normal[i]
            else:
                normal_outward[perm[i]] = bdry_normal[perm[i]]

        #   d). Asignacion de condicion
        for pos in condition_order_list:  
            pts_list = colloc_l1[ colloc_l2[pos]+1 : colloc_l2[pos+1]+1 ]
            size = len(pts_list)
            for i in range(colloc_l2[pos]+1,size+colloc_l2[pos]+1):
                # asignacion de las condiciones de contorno
                bdry_condition[perm[colloc_l1[i]]] = \
                    self.boundary_condition_assignment(pos,colloc_l1[i])
                # asignacion de las normales
                if ( colloc_l1[i] in vrtx_normal.keys() ):
                    normal_outward[perm[colloc_l1[i]]] = vrtx_normal[colloc_l1[i]]
                else:
                    normal_outward[perm[colloc_l1[i]]] = bdry_normal[i]

        # global variables
        bdry = {}
        for a,b,c,d in zip(bdry_pts,bdry_condition,bdry_coord,normal_outward):
            bdry[a+1]=(b,c,d)

        # se utiliza actualiza la variable self.bdry
        self.bdry = bdry

        return
            
    """
    FUNCION AUXLIAR 1: 'set_bdry_condition'
      outputs: 
          - cond,  condicion de contorno asociado a cada
                   segun tabla de arriba (cond=1,2,3 o 4)
    """
    def boundary_condition_assignment(self,lado,pt):

        if (self.test=='patchtest-0'):
            cond = 3       # only uv-displacement condition
        elif (self.test=='patchtest-1'):
            ## ------ segun lado ------ ##
            if ( lado == 0 ):   # borde inferior 
                cond = 2       #   -> only_v_displacement
            elif ( lado == 1 ): # borde derecho
                cond = 4       #   -> tau_traction
            elif ( lado == 2 ): # borde superior 
                cond = 4       #   -> tau_traction (zero)
            elif ( lado == 3 ): # borde izquierdo
                cond = 1       #   -> only_u_displacement
            ## ----- segun un punto en particular ----- ##
            if ( pt==0 ): # vertice de simetria (0,0)
                cond = 3  #   -> both_uv_displacement
        elif (self.test=='patchtest-2'):
            ## ------ segun lado ------ ##
            if ( lado == 0 ):   # borde inferior 
                cond = 2       #   -> only_v_displacement
            elif ( lado == 1 ): # borde derecho
                cond = 4       #   -> tau_traction
            elif ( lado == 2 ): # borde superior 
                cond = 4       #   -> tau_traction (zero)
            elif ( lado == 3 ): # borde izquierdo
                cond = 1       #   -> only_u_displacement
            ## ----- segun un punto en particular ----- ##
            if ( pt==0 ): # vertice de simetria (0,0)
                cond = 3  #   -> both_uv_displacement
        else:
            print('ERROR: Test no registrado')
            from sys import exit; exit('Program stopped')

        return cond

    """
    FUNCION AUXLIAR 2: 'set_bdry_condition'
      outputs: 
          - cond,  condicion de contorno asociado a cada
                   segun tabla de arriba (cond=1,2,3 o 4)
    """
    def order_boundary_list(self):
        # order list : contiene las id del contorno ordenadas de tal manera
        #              las  primeras  entradas  son  reemplazadas  por  las 
        #              subsiguientes.( ultimas entradas poseen prioridad )
        if ( self.test == 'patchtest-0' ):
            #     _III_       
            #    |     |     tanto en I, II, III y IV se preescriben 
            # IV |     | II  desplazamientos en direccion x e y (u y v)
            #    |_____|      
            #       I
            priority_order_list = [ 0 , 1 , 2 , 3 ]  # orden irrelevante
        elif ( self.test == 'patchtest-1' ):
            #     _III_       en I   : preescribe v=0 y dot(sigma,n)_x = 0
            #    |     |      en II  : preescribe dot(sigma,n) = traction
            # IV |     | II   en III : preescribe dot(sigma,n) = traction = 0
            #    |_____|      en IV  : preescribe u=0 y dot(sigma,n)_y = 0
            #       I
            priority_order_list = [ 1 , 2 , 0 , 3 ]
        elif ( self.test == 'patchtest-2' ):
            #     _III_       en I   : preescribe v=0 y dot(sigma,n)_x = 0
            #    |     |      en II  : preescribe dot(sigma,n) = traction
            # IV |     | II   en III : preescribe dot(sigma,n) = traction = 0
            #    |_____|      en IV  : preescribe u=0 y dot(sigma,n)_y = 0
            #       I
            priority_order_list = [ 1 , 2 , 0 , 3 ]
        else:
            print('ERROR, test no ingresado')
            from sys import exit ; exit('Program stopped')

        return priority_order_list
    
    """
    GENERAR_RHS
        genera el lado derecho de la ecuacion (rhs), de la forma K u = rhs    
    """
    def GENERAR_RHS(self):
        
        # iniciar variables
        cont=0 ; self.rhs={}
        
        # carga libreria numpy
        from numpy import array,dot

        # recorre cada punto del dominio
        for pts in self.bdry.keys():

            # extract bdry data
            cond,(x,y,z),(nx,ny,nz) =self.bdry[pts] 

            # variables elasticidad lineal (linear elastic deformation)
            from exactsol import SOLUCION_EXACTA
            u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy=SOLUCION_EXACTA(\
                self.test,x,y,E=self.young_modulus,v=self.poisson_coeff)

            # compute Traction from Displacement
            mat = array([[sigmax,tauxy],[tauxy,sigmay]])
            normal = array([nx,ny])
            Tx,Ty = list( dot( mat , normal ) )
            self.rhs[pts] = ((u,v),(Tx,Ty))

        return
        


    """
        EXPORTAR_DATA

    ejemplo de archivo generado: "test-name.dat"

   NODOS  D.PRESCRITOS F.PRESCRITAS 
     9         8         0
 PUNTOS     COORDENADAS
     1  0.000000E+00  0.000000E+00
     2  0.100000E+01  0.000000E+00
     3  0.000000E+00  0.100000E+01
     4  0.100000E+01  0.100000E+01
     5  0.200000E+01  0.000000E+00
     6  0.000000E+00  0.200000E+01
     7  0.100000E+01  0.200000E+01
     8  0.200000E+01  0.100000E+01
     9  0.200000E+01  0.200000E+01
 DESPLAZAMIENTO PRESCRITO
     1     1     1  0.000000E+00  0.000000E+00
     2     1     1  0.100000E+01  0.100000E+01
     3     1     1  0.100000E+01  0.100000E+01
     5     1     1  0.200000E+01  0.200000E+01
     6     1     1  0.200000E+01  0.200000E+01
     7     1     1  0.300000E+01  0.300000E+01
     8     1     1  0.300000E+01  0.300000E+01
     9     1     1  0.400000E+01  0.400000E+01
 FUERZAS PRESCRITAS
 NORMALES
     8
     1 -0.707107E+00 -0.707107E+00
     2  0.000000E+00 -0.100000E+01
     3 -0.100000E+01  0.000000E+00
     5  0.707107E+00 -0.707107E+00
     6 -0.707107E+00  0.707107E+00
     7  0.000000E+00  0.100000E+01
     8  0.100000E+01  0.000000E+00
     9  0.707107E+00  0.707107E+00
END
    """
    def EXPORTAR_DATA(self):
        # --- crea carpeta si es que no existe ---
        import os
        if ( not os.path.isdir(self.data_dir) ):
                os.mkdir(self.data_dir)
        
        # --- exporta los datos ---
        f = open( self.data_dir + self.test+'_'+self.string_npts+'.dat' , 'w' )

        n = len(self.coord)
        
        cond_list = []
        for i in self.bdry.keys():
            cond_list.append(self.bdry[i][0])
        n_dprescrt = cond_list.count(1)+cond_list.count(2)+cond_list.count(3)
        n_fprescrt = cond_list.count(4)

        # -------------------------
        f.write('{0:>10s} {1:>15s} {2:>15s}\n'.format('NODOS','D.PRESCRITOS','F.PRESCRITAS'))
        f.write('{0:>10d} {1:>15d} {2:>15d}\n'.format(n,n_dprescrt,n_fprescrt))

        # -------------------------
        f.write('{0:>10s} {1:>15s}\n'.format('PUNTOS','COORDENADAS'))
        x,y = zip(*self.coord)
        for i in range(n):
            f.write('{0:10d} {1:20.13f} {2:20.13f}\n'.format(i+1,x[i],y[i]))

        # -------------------------
        f.write('{0:<30s}\n'.format('DESPLAZAMIENTO PRESCRITO'))
        for pts in self.bdry.keys():
            # extraer data
            cond,xcoord,xnorm = self.bdry[pts]
            (u,v),(Tx,Ty) = self.rhs[pts]
            if ( cond == 1 ):
                f.write('{0:5d} {1:5d} {2:5d} {3:20.13f} {4:20.13f}\n'.format(pts,1,0,u,v))
            elif ( cond == 2 ):
                f.write('{0:5d} {1:5d} {2:5d} {3:20.13f} {4:20.13f}\n'.format(pts,0,1,u,v))
            elif ( cond == 3 ):
                f.write('{0:5d} {1:5d} {2:5d} {3:20.13f} {4:20.13f}\n'.format(pts,1,1,u,v))

        # -------------------------
        f.write('{0:<30s}\n'.format('FUERZAS PRESCRITAS'))
        for pts in self.bdry.keys():
            # extraer data
            cond,xcoord,xnorm = self.bdry[pts]
            (u,v),(Tx,Ty) = self.rhs[pts]
            if ( cond == 4 ):
                f.write('{0:5d} {1:25.13e} {2:25.13e}\n'.format(pts,Tx,Ty))

        # -------------------------
        f.write('{0:<30s}\n'.format('NORMALES'))
        f.write('{0:5d}\n'.format(len(self.bdry)))
        for pts in self.bdry.keys():
            # extraer data
            cond,xcoord,(nx,ny,nz) = self.bdry[pts]
            f.write('{0:5d} {1:20.13f} {2:20.13f}\n'.format(pts,nx,ny))

    
        # -------------------------
        f.write('END')
        f.close()
    

    def EXPORTAR_GENDATA(self):
        f = open( self.data_dir + self.test+'_'+self.string_npts+'.GEN' , 'w' )
        f.write('IWEIGT\n')
        f.write('0\n')
        f.write('BET\n')
        f.write('0.0e+00\n')
        f.write('PROPIEDADES DEL MATERIAL\n')
        f.write('MODULO DE YOUNG\n')
        f.write('1.0e+03\n')
        f.write('MODULO DE POISSON\n')
        f.write('0.3e+00\n')
        return


    def EXPORTAR_ELEMENTOS(self):
        # --- crea carpeta si es que no existe ---
        import os
        if ( not os.path.isdir(self.data_dir) ):
                os.mkdir(self.data_dir)
        
        # --- exporta los datos ---
        f = open( self.data_dir + self.test+'_'+self.string_npts+'.elem' , 'w' )

        n = len(self.coord)
        
        ntria = len(self.tri)
        nline = len(self.bdry_elem)

        # -------------------------

        f.write('{0:>10s} \n'.format('N. Triangulos'))
        f.write('{0:>10d}\n'.format(ntria))
        for pts in self.tri:
            f.write('{0:10d} {1:10d} {2:10d}\n'.format(pts[0]+1,pts[1]+1,pts[2]+1))

        # -------------------------

        f.write('{0:>10s} \n'.format('N. Lineas'))
        f.write('{0:>10d}\n'.format(nline))
        for pts in self.bdry_elem:
            f.write('{0:10d} {1:10d}\n'.format(pts[0]+1,pts[1]+1))
    
        # -------------------------
        f.write('END')
        f.close()
        return



    def EXPORTAR_el2pt(self):
        # --- crea carpeta si es que no existe ---
        import os
        if ( not os.path.isdir(self.data_dir) ):
                os.mkdir(self.data_dir)
        
        # ----- preambulo --------------------
        npts = len(self.coord)
        ntria = len(self.tri)
        nline = len(self.bdry_elem)

        nfix = 0 ; nload=0 ; nbound=0
        displ=[] ; tracc=[]
        for i in self.bdry:
            a,b,c = self.bdry[i]
            if ( a==1 or a==2 or a==3 ):
                nfix += 1 
            elif ( a==4 ): 
                nbound += 1 
        nbound = nfix+nload

        #:::::::::::::::::::::::::::::::::::

        # --- exporta los datos ---
        f = open(self.data_dir+self.test+'_'+self.string_npts+'_el2pt.dat','w')

        f.write('\n')
        f.write('\n')
        f.write('boundary conditions\n')
        f.write('{0:>10s} {1:>10s} {2:>10s}\n'.format('nfix','nload','nbound'))
        f.write('{0:>10d} {1:>10d} {2:>10d}\n'.format(nfix,nload,nbound))
        f.write('{0:>10s} {1:>10s}\n'.format('nele','npoints'))
        f.write('{0:>10d} {1:>10d}\n'.format(ntria,npts))

        f.write('NODES\n')
        idpt = 0 ; formato = '{0:>10d} {1:>20.13e} {2:>20.13e}\n'
        for xpt in self.coord:
            idpt+=1 
            f.write(formato.format(idpt,xpt[0],xpt[1]))

        f.write('ELEMENTS\n')
        idpt = 0 ; formato='{0:10d} {1:10d} {2:10d} {3:10d}\n'
        for elem in self.tri:
            idpt +=1 
            f.write(formato.format(idpt,elem[0]+1,elem[1]+1,elem[2]+1))


        f.write('DESPLAZAMIENTOS\n')
        formato = '{0:5d} {1:5d} {2:5d} {3:20.13f} {4:20.13f}\n'
        for pts in self.bdry.keys():
            cond,xcoord,xnorm = self.bdry[pts]
            (u,v),(Tx,Ty) = self.rhs[pts]
            if ( cond == 1 ):
                f.write(formato.format(pts,1,0,u,v))
            elif ( cond == 2 ):
                f.write(formato.format(pts,0,1,u,v))
            elif ( cond == 3 ):
                f.write(formato.format(pts,1,1,u,v))

        f.write('FUERZAS\n')
        formato = '{0:5d} {1:20.13f} {2:20.13f}\n'
        for pts in self.bdry.keys():
            cond,xcoord,xnorm = self.bdry[pts]
            if ( cond == 4 ):
                f.write(formato.format(pts,Tx,Ty))
        f.write('\n')

        f.close()

        #:::::::::::::::::::::::::::::::::::

        f = open(self.data_dir+self.test+'_'+self.string_npts+'_el2pt-1.dat','w')

        idpt = 0 ; formato='{0:10d} {1:10d}\n'
        for elem in self.bdry_elem:
            idpt +=1 
            f.write(formato.format(elem[0]+1,elem[1]+1))

        f.close()

        return



    def PLOT_PREPROCESS_RESULTS(self):

        # libreria matplotlib
        import matplotlib.pyplot as plt

        # create 'figure' and 'axis' matplotlib objects
        fig,ax = plt.subplots(figsize=(10,8))
        
        # formato fuente
        plt.rc('text',usetex=True)
        plt.rc('font',family='serif')

        # coordenadas
        x,y = zip(*self.coord)

        # axis settings
        xmax=max(x) ; xmin=min(x) ; ymax=max(y) ; ymin=min(y)
        percent = 0.2 # 20%
        dx=(xmax-xmin)*percent ; dy=(ymax-ymin)*percent 
        ax.set_xlim(xmin-dx,xmax+dx)
        ax.set_ylim(ymin-dy,ymax+dy)
        ax.set_aspect('equal')

        # plot triangulation
        plt.triplot(x,y,triangles=self.tri,color='gray',linestyle='-',linewidth=0.2)

        # plot domain pts
        ax.scatter(x,y,color='black')

        # plot contorno
        for elem in self.bdry_elem:
            aux = []
            for pts in elem:
                aux.append( self.coord[pts] )
            xline,yline=zip(*aux)
            ax.plot(xline,yline,color='black',linestyle='--',linewidth='0.5') 

        # plot id pts
        percent = 0.01 
        dx=(xmax-xmin)*percent ; dy=(ymax-ymin)*percent 
        for i in range(len(self.coord)):
           ax.text(x[i]+dx,y[i]+dy,s=i+1,color='black') 

        # extraer data ( diccionario a lista )
        aux = [] ; bdry_pts =[]
        for i in self.bdry:
            bdry_pts.append(i)
            aux.append(self.bdry[i])
        cond,bdry_coord,bdry_normal = zip(*aux)
        xbdry,ybdry,zbdry = zip(*bdry_coord)
        xnorm,ynorm,zbdry = zip(*bdry_normal)

        # crear listas segun condicion
        dict_cond = { 1:[] , 2:[] , 3:[] , 4:[] }
        dict_color = {1:'blue',2:'green',3:'brown',4:'red'}
        dict_label = {1:r'Desplazamiento $(1,0)$',2:r'Desplazamiento $(0,1)$'\
                        ,3:r'Desplazamiento $(1,1)$',4:r'Esfuerzos $(T_x,T_y)$'}

        for i in range(len(bdry_pts)):
            if ( not cond[i] in [1,2,3,4]  ) :
                from sys import exit
                exit('ERROR: id contorno no registrada.')
            dict_cond[cond[i]].append(bdry_pts[i])

            if ( cond[i] == 4 ) :
                ax.quiver(xbdry[i],ybdry[i],xnorm[i],ynorm[i])
                
        for idcon in dict_cond:
            xcon=[] ; ycon=[] ; cont=0
            for pt in dict_cond[idcon]:
                xcon.append(x[pt-1]) ; ycon.append(y[pt-1]) ; cont+=1
            if ( cont!=0 ):
                ax.scatter(xcon,ycon,color=dict_color[idcon],label=dict_label[idcon])

        # muestra labels
        ax.legend()

        # --- crea carpeta si es que no existe ---
        import os
        if ( not os.path.isdir(self.plot_dir) ):
                os.mkdir(self.plot_dir)

        # aqui hay un error
        export_file = self.plot_dir+self.test+self.formato
        fig.savefig(export_file)
        print('Figura exportada: '+export_file)

        if (self.show):
            plt.show()

        plt.cla()    # clear axis
        plt.clf()    # clear figure
        plt.close()  # close a figure windows

        return
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#   Lectura de parametros via argumentos (Comando), utiliza argsparse
import argparse
import os

# parse arguments from bash/batch script
parser = argparse.ArgumentParser(description='Preproceso utilizado por \
         programa de resolucion numerica de ecuaciones diferenciales \
         utilizando metodos sin malla de colocacion puntual')
parser.add_argument('-f','--func', required=True, help="(str) Funcion a ejecutar: 'Preproceso'/'Graficar_Solucion_Exacta'")
parser.add_argument('-t','--test', required=True, help="(str) Ingresar nombre del test")
parser.add_argument('-n','--npts', help="(str) Numero de puntos de colocacion en la direccion x e y ('xcollocnpt-ycollocnpt')")
parser.add_argument('-s','--show', help="(str) Muestra graficos interactivos en terminal")
parser.add_argument('-p','--plot', help="(str) Exporta figura en ./GRAFICOS/")
parser.add_argument('-E','--young', help="(str) Modulo de Young")
parser.add_argument('-v','--poisson', help="(str) Coeficiente de Poisson")

# lectura argumentos, se almacena en 'args'
args = parser.parse_args()

# verificacion de carpeta contenedora. las crea si es que no existen
directorios=["./DATOS/","./GRAFICOS/"]
for dirs in directorios:
    if ( not os.path.isdir(dirs) ):
        os.mkdir(dirs)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#   1.    ingreso datos de entrada
meshfree_test2d_collocation_data = Preprocess_Data_Type(\
        args.test,args.npts,args.show,args.plot,args.young,args.poisson)
#   2.    generar informacion del dominio y su contorno
#           -   coordenadas puntos de discretizacion
#           -   triangulacion de elementos 
#           -   identificacion de puntos de contorno y 
#           -   su vector normal asociado
meshfree_test2d_collocation_data.GENERAR_GEOMETRIA_DOMINIO()
#   3.    establecer condiciones de contorno
meshfree_test2d_collocation_data.PREESCRIPCION_CONDICIONES_CONTORNO()
#   4.    generar datos rhs, usados en solver
#         resolucion de un sistema de ecuaciones lineales
#         A + x = b
#         utilizando la nomenclatura meshfree
#         [ K_matrix ] * u_vector = right_hand_side_vector (rhs)
meshfree_test2d_collocation_data.GENERAR_RHS()
#   5.    exportar data
meshfree_test2d_collocation_data.EXPORTAR_DATA()
#   6.    exportar informacion necesario para el solver
meshfree_test2d_collocation_data.EXPORTAR_GENDATA()
#   7.    exportar elementos del dominio
meshfree_test2d_collocation_data.EXPORTAR_ELEMENTOS()
#   8.    exportar datos de entrada para el2pt.f90
meshfree_test2d_collocation_data.EXPORTAR_el2pt()
#   8.  Plotear data (opcional)
if (args.plot!=None):
    if (args.plot.lower()=="true"):
        meshfree_test2d_collocation_data.PLOT_PREPROCESS_RESULTS()

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

