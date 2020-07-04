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
        #     _III_       en I   : preescribe v=0 y dot(sigma,n)_x = 0
        #    |     |      en II  : preescribe dot(sigma,n) = traction
        # IV |     | II   en III : preescribe dot(sigma,n) = traction = 0
        #    |_____|      en IV  : preescribe u=0 y dot(sigma,n)_y = 0
        #       I
        priority_order_list = [ 1 , 2 , 0 , 3 ]
    elif ( test == 'patchtest-2' ):
        #     _III_       en I   : preescribe v=0 y dot(sigma,n)_x = 0
        #    |     |      en II  : preescribe dot(sigma,n) = traction
        # IV |     | II   en III : preescribe dot(sigma,n) = traction = 0
        #    |_____|      en IV  : preescribe u=0 y dot(sigma,n)_y = 0
        #       I
        priority_order_list = [ 1 , 2 , 0 , 3 ]
    elif ( test=='cantilever'):
        #    _______III______       en I   : preescribe dot(sigma,n) = traction = 0
        #    |               |      en II  : preescribe dot(sigma,n) = traction
        # IV |               | II   en III : preescribe dot(sigma,n) = traction = 0
        #    1_______________|      en IV  : preescribe u=0 y v=0
        #            I
        priority_order_list = [ 0 , 1 , 2 , 3 ]
    else:
        print('ERROR, test no ingresado')
        from sys import exit ; exit('Program stopped')

    return priority_order_list

