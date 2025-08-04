from lrgsglib.core import *
#
move_to_rootf(print_tf=True)
#
navrg=50
shapeflippin = 'single'
p1=np.linspace(0.005,0.06,8)
p2=np.linspace(0.06,0.2,20)
p3=np.linspace(0.2,0.5,7)
p_comb = np.concatenate([p1, p2, p3])
#
for N in [32, 64, 128]:
    #
    for p in p_comb:
        #
        Pinf=np.zeros(navrg)
        Pinf2=np.zeros(navrg)
        Fluct=np.zeros(navrg)
        Fluct2=np.zeros(navrg)
        #
        lattice = Lattice2D(N, pflip=p, geometry='triangular')
        file_out = f'{lattice.path_phtra.name}{shapeflippin}_p_{lattice.pflip:.3g}_{navrg}'
        if os.path.exists(file_out):
            continue
        #
        for cont, avrg in tqdm(enumerate(range(navrg))):
            lattice = Lattice2D(N, pflip=p, geometry='triangular')
            file1 = f'{lattice.path_phtra.name}{shapeflippin}_p_{lattice.pflip:.3g}_{avrg+1}'
            try:
                file2 = f'{lattice.path_phtra.name}{shapeflippin}_p_{lattice.pflip:.3g}_{avrg}'
                os.remove(file2)
            except FileNotFoundError:
                pass
            if shapeflippin == 'triangle':
                lattice.flip_sel_edges(lattice.neg_weights_dict.NEG_WEIGHTS_DICT_H_PTRIA)
            elif shapeflippin == 'single':
                lattice.flip_random_fract_edges()
            lattice.compute_k_eigvV()
            lattice.calc_fluct_Pinf()

            Fluct[cont]=lattice.eigV_fluct
            Fluct2[cont]=lattice.eigV_fluct**2
            
            Pinf[cont]=lattice.Pinf
            Pinf2[cont]=lattice.Pinf**2
            
            x=[np.sum(Pinf)/(1+avrg),np.sum(Fluct)/(1+avrg), np.sum(Fluct2)/(1+avrg), np.var(Fluct[Fluct!=0]), np.sum(Pinf2)/(1+avrg),lattice.pflip,int(avrg+1)]
            np.savetxt(file1, np.atleast_2d(x), fmt='%.5g', delimiter=',')
