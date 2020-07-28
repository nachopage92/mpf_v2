"""
FUNCION AUXLIAR 1: 'set_bdry_condition'
  outputs: 
      - cond,  condicion de contorno asociado a cada
               segun tabla de arriba (cond=1,2,3 o 4)
"""
def boundary_condition_assignment(test,lado,pt):

    if ( test=='patchtest-0' ):
        cond = 3       # only uv-displacement condition
    elif ( test=='patchtest-1' ):
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
    elif ( test=='patchtest-2' ):
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
    elif ( test=='cantilever' ):
        ## ------ segun lado ------ ##
        if ( lado == 0 ):   # borde inferior 
            cond = 4       #   -> tau_traction (zero)
        elif ( lado == 1 ): # borde derecho
            cond = 4       #   -> tau_traction
        elif ( lado == 2 ): # borde superior 
            cond = 4       #   -> tau_traction (zero)
        elif ( lado == 3 ): # borde izquierdo
            cond = 3       #   -> uv_displacement
    elif ( test=='infinite-plate' ):
        if ( lado == 0):
            cond = 2
        if ( lado == 1 ):
            cond = 4
        if ( lado == 2 ):
            cond = 4
        if ( lado == 3 ):
            cond = 1
        if ( lado == 4 ):
            cond = 4
    else:
        print('boundary_condition_assignment():')
        print('ERROR: Test no registrado')
        from sys import exit; exit('Program stopped')

    return cond

"""
FUNCION AUXLIAR 2: 'set_bdry_condition'
  outputs: 
      - cond,  condicion de contorno asociado a cada
               segun tabla de arriba (cond=1,2,3 o 4)
"""
def order_boundary_list(test):
    # order list : contiene las id del contorno ordenadas de tal manera
    #              las  primeras  entradas  son  reemplazadas  por  las 
    #              subsiguientes.( ultimas entradas poseen prioridad )
    if ( test == 'patchtest-0' ):
        #     _III_       
        #    |     |     tanto en I, II, III y IV se preescriben 
        # IV |     | II  desplazamientos en direccion x e y (u y v)
        #    |_____|      
        #       I
        priority_order_list = [ 0 , 1 , 2 , 3 ]  # orden irrelevante
    elif ( test == 'patchtest-1' ):
        #     _III_       en I   : preesc. v=0 y dot(sigma,n)_x = tx
        #    |     |      en II  : preesc. dot(sigma,n) = traction
        # IV |     | II   en III : preesc. dot(sigma,n) = traction
        #    |_____|      en IV  : preesc. u=0 y dot(sigma,n)_y = ty
        #       I
        priority_order_list = [ 1 , 2 , 0 , 3 ]
    elif ( test == 'patchtest-2' ):
        #     _III_       en I   : preesc. v=0 y dot(sigma,n)_x = tx
        #    |     |      en II  : preesc. dot(sigma,n) = traction
        # IV |     | II   en III : preesc. dot(sigma,n) = traction
        #    |_____|      en IV  : preesc. u=0 y dot(sigma,n)_y = ty
        #       I
        priority_order_list = [ 1 , 2 , 0 , 3 ]
    elif ( test=='cantilever'):
        #    _______III______       en I   : preesc. dot(sigma,n) = traction
        #    |               |      en II  : preesc. dot(sigma,n) = traction
        # IV |               | II   en III : preesc. dot(sigma,n) = traction
        #    1_______________|      en IV  : preesc. u=0 y v=0
        #            I
        priority_order_list = [ 0 , 1 , 2 , 3 ]
    elif ( test=='infinite-plate'):
        #    ___III___       en I   : preesc. v=0 y dot(sigma,n)_x = tx
        # IV |       |      en II  : preesc. dot(sigma,n) = traction
        #    --\     |II    en III : preesc. dot(sigma,n) = traction
        #     V |____|      en IV  : preesc. u=0 y dot(sigma,n)_y = ty
        #         I         en V   : preesc. dot(sigma,n) = traction 
        priority_order_list = [ 4 , 2 , 1 , 0 , 3 ]
    else:
        print(' order_boundary_list(): ')
        print('ERROR, test no ingresado')
        from sys import exit ; exit('Program stopped')

    return priority_order_list


def additional_info( test ):
    # ingresa subcontornos en caso de haberlos
    edge_tol   = 1e-12
    refinement = [0,0]
    subs     = [[],[],[],[]]
    if ( test == 'infinite-plate' ):
        subs = [[],[0.5],[],[]]
        refinement = [0,0]
    return subs, edge_tol, refinement

def input_knots( test , npts ):
    from numpy import linspace,sin,pi
    nx,ny = npts
    uknots,vknots = [linspace(0.0,1.0,nx),linspace(0.0,1.0,ny)]
    if (test=='infinite-plate'):
        vknots = sin(linspace(-0.5*pi,0.5*pi,ny))
        vknots = [ 0.5*(1.0+val) for val in vknots ]

    return [uknots,vknots]

    
    
