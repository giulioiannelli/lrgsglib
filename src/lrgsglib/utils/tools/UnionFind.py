class UnionFind:
    """
    A Union-Find or Disjoint Set Union (DSU) data structure with path compression and union by rank.
    
    It provides an efficient way to manage a partition of a set into disjoint subsets and is useful 
    for dealing with connectivity queries, particularly in graph algorithms.

    Attributes:
    -----------
    parent (List[int]): Stores the parent of each element. Initially, each element is its own parent.
    rank (List[int]): Represents the rank of each element, used to keep the tree flat by attaching 
                       the root of the smaller tree under the root of the larger tree.
    """

    def __init__(self, n):
        """
        Initializes the UnionFind structure with n elements.

        Parameters:
        -----------
        n (int): The number of elements in the Union-Find structure.
        
        Returns:
        --------
        None
        """
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, p):
        """
        Finds the representative of the set containing 'p' with path compression.
        Path compression flattens the structure of the tree whenever `find` is used,
        so that future operations also benefit from the flatter tree.

        Parameters:
        -----------
        p (int): The element whose set representative is to be found.
        
        Returns:
        --------
        int: The representative of the set containing 'p'.
        """
        if self.parent[p] != p:
            self.parent[p] = self.find(self.parent[p])
        return self.parent[p]

    def union(self, p, q):
        """
        Merges the sets containing 'p' and 'q'. It uses the union by rank strategy,
        which attaches the tree with less depth (smaller rank) under the root of the deeper tree
        (larger rank) to keep the tree flat.

        Parameters:
        -----------
        p (int): An element of the first set.
        q (int): An element of the second set.
        
        Returns:
        --------
        None
        """
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP != rootQ:
            if self.rank[rootP] > self.rank[rootQ]:
                self.parent[rootQ] = rootP
            elif self.rank[rootP] < self.rank[rootQ]:
                self.parent[rootP] = rootQ
            else:
                self.parent[rootQ] = rootP
                self.rank[rootP] += 1