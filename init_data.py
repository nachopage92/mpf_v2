def nurbs_geom_data(test):

    if (test=='patchtest-0'):
        degree_u = 1
        degree_v = 1
        dim_u = 2
        dim_v = 2
        knotvector_u = [0.0,0.0,1.0,1.0]
        knotvector_v = [0.0,0.0,1.0,1.0]
        ctrlpts  = [[0.0,0.0],\
                    [0.0,2.0],\
                    [2.0,0.0],\
                    [2.0,2.0]]
        weight = [1.0,1.0,1.0,1.0]
    
    elif (test=='patchtest-1'):
        degree_u = 1
        degree_v = 1
        dim_u = 2
        dim_v = 2
        knotvector_u = [0.0,0.0,1.0,1.0]
        knotvector_v = [0.0,0.0,1.0,1.0]
        ctrlpts  = [[0.0,0.0],\
                    [0.0,1.0],\
                    [1.0,0.0],\
                    [1.0,1.0]]
        weight = [1.0,1.0,1.0,1.0]
    
    elif (test=='patchtest-2'):
        degree_u = 1
        degree_v = 1
        dim_u = 2
        dim_v = 2
        knotvector_u = [0.0,0.0,1.0,1.0]
        knotvector_v = [0.0,0.0,1.0,1.0]
        ctrlpts  = [[0.0,0.0],\
                    [0.0,1.0],\
                    [1.0,0.0],\
                    [1.0,1.0]]
        weight = [1.0,1.0,1.0,1.0]
    
    elif (test=='cantilever'):
        degree_u = 1
        degree_v = 1
        dim_u = 2
        dim_v = 2
        knotvector_u = [0.0,0.0,1.0,1.0]
        knotvector_v = [0.0,0.0,1.0,1.0]
        ctrlpts  = [[0.0, -0.5],\
                    [0.0,0.5],\
                    [4.0,-0.5],\
                    [4.0,0.5]]
        weight = [1.0,1.0,1.0,1.0]
    
    elif (test=='infinite-plate'):
        from numpy import sqrt
        degree_u = 2
        degree_v = 2
        dim_u = 3 
        dim_v = 4
        knotvector_u = [0.0,0.0,0.0,1.0,1.0,1.0]
        knotvector_v = [0.0,0.0,0.0,0.5,1.0,1.0,1.0]
        ctrlpts =  [[1.0,0.0],\
                    [1.0,sqrt(2.0)-1.0],\
                    [sqrt(2.0)-1.0,1.0],\
                    [0.0,1.0],\
                    [2.5,0.0],\
                    [2.5,0.75],\
                    [0.75,2.5],\
                    [0.0,2.5],\
                    [4.0,0.0],\
                    [4.0,4.0],\
                    [4.0,4.0],\
                    [0.0,4.0]]
        weight =  [1.0,0.5*(1.0+1.0/sqrt(2.0)),0.5*(1.0+1.0/sqrt(2.0)),1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]


    else:
        print('test no ingresado:')
        from sys import exit ; exit('programa detenido!')

    ctrlptsw = []
    for pts,w in zip(ctrlpts,weight):
        aux = [] 
        for pt in pts:
            aux.append(pt*w)
        aux.append(w)
        ctrlptsw.append(aux)
    
    return [degree_u,degree_v],[dim_u,dim_v],[knotvector_u,knotvector_v],ctrlptsw
