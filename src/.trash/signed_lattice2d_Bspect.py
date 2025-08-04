from lrgsglib.core import *

plist = np.concatenate([[0.001, 0.01, 0.05, 0.075, 0.080, 0.085, 0.090, 0.091], 
                np.linspace(0.090, 0.106, num=17),
                np.linspace(0.12, 1, num=9)])
lattice2d = nx.grid_2d_graph(48, 48, periodic=True)
os.makedirs("data/mat_spectra/", exist_ok=True)
for p in plist:
    saveFile = open(f"data/mat_spectra/p={p:.3g}.bin", "wb")
    for na in tqdm(range(1000)):
        adj = nx.adjacency_matrix(lattice2d).toarray()
        indx1 = np.where(adj == 1)
        #
        total_indices = len(indx1[0])
        num_indices_to_select = int(p * total_indices)
        random_indices = np.random.choice(total_indices, num_indices_to_select, replace=False)
        # Get the randomly selected indices from the original indices array
        selected_indices = (indx1[0][random_indices], indx1[1][random_indices])
        adj[selected_indices] = 0
        graph = nx.from_numpy_array(adj)
        # Compute the Laplacian matrix
        laplacian_matrix = nx.laplacian_matrix(graph).toarray()
        saveFile.write(bytes(eigvalsh(laplacian_matrix)))
    saveFile.close()