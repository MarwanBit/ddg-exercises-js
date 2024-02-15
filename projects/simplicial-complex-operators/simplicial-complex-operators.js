"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // iterate through the vertices, edges, and faces and assign the indices
                for (let i = 0; i < mesh.vertices.length; i++) {
                        mesh.vertices[i].index = i;
                }
                for (let i = 0; i < mesh.edges.length; i++) {
                        mesh.edges[i].index = i;
                }
                for (let i = 0; i < mesh.faces.length; i++) {
                        mesh.faces[i].index = i;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                //We create an adjacency matrix of size |E| x |V|
                // we first create a triplet 
                let T = new Triplet(mesh.edges.length, mesh.vertices.length);
                // now we iterate through the edges
                for (let e of mesh.edges) {
                        // retrieve the indices of the vertices belonging to the edge
                        let v1Index = e.halfedge.vertex.index;
                        let v2Index = e.halfedge.next.vertex.index;
                        // add 1's to the sparse matrix 
                        T.addEntry(1, e.index, v1Index);
                        T.addEntry(1, e.index, v2Index);
                }
                let VertexEdgeAdjacencyMatrix = SparseMatrix.fromTriplet(T);
                return VertexEdgeAdjacencyMatrix;
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                // We create an adjacency matrix of size |F| x |V|
                let T = new Triplet(mesh.faces.length, mesh.edges.length);
                // now we iterate through each face
                for (let face of mesh.faces) {
                        // iterate through each edge of each face
                        for (let edge of face.adjacentEdges()) {
                                let edgeIndex = edge.index;
                                T.addEntry(1, face.index, edgeIndex);
                        }
                }
                let EdgeFaceAdjacencyMatrix = SparseMatrix.fromTriplet(T);
                return EdgeFaceAdjacencyMatrix;
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                // first we must create our vertex vector
                let VertexVector = DenseMatrix.zeros(this.mesh.vertices.length, 1);
                // now iterate through each vertex in the subset and update VertexVector
                for (let vertexIndex of subset.vertices) {
                        VertexVector.set(1, vertexIndex, 0);
                }
                return VertexVector;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                // first we must create our edge vector
                let EdgeVector = DenseMatrix.zeros(this.mesh.edges.length, 1);
                // now iterate through each edge in the subset and update VertexVector
                for (let edgeIndex of subset.edges) {
                        EdgeVector.set(1, edgeIndex, 0);
                }
                return EdgeVector;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                // first we must create our face vector
                let FaceVector = DenseMatrix.zeros(this.mesh.faces.length, 1);
                // now iterate through each face in the subset and update VertexVector
                for (let faceIndex of subset.faces) {
                        FaceVector.set(1, faceIndex, 0);
                }
                return FaceVector;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // create a star subset which starts as the current subset
                // now we create our star object
                let star = MeshSubset.deepCopy(subset);;
                
                //now all we must do is iterate through each vertex and edge in our subset
                // and append the adjacent edges and face
                for (let i in Array.from(subset.vertices)) {
                        let vIndex = Array.from(subset.vertices)[i];
                        let v = this.mesh.vertices[vIndex];
                        // now we iterate through the edges adjacent to this vertex;
                        // now we create an array for the edges
                        for (let adjacentEdgeIndex in Array.from(v.adjacentEdges())) {
                                let adjacentEdge = Array.from(v.adjacentEdges())[adjacentEdgeIndex];
                                star.addEdge(adjacentEdge.halfedge.edge.index);
                        }   
                        // add all adjacent faces 
                        for (let adjacentFaceIndex in Array.from(v.adjacentFaces())) {
                                let adjacentFace = Array.from(v.adjacentFaces())[adjacentFaceIndex];
                                star.addFace(adjacentFace.halfedge.face.index);
                        }
                }

                // now for edges we must add all adjacent vertices and faces
                for (let i in Array.from(subset.edges)) {
                        let eIndex = Array.from(subset.edges)[i];
                        let e = this.mesh.edges[eIndex];
                        // now we iterate through the edges adjacent to this vertex;
                        let f1Index = e.halfedge.face.index;
                        star.addFace(f1Index);

                        // Make sure we are not on the bounday
                        if (!e.onBoundary()){
                                let f2Index = e.halfedge.twin.face.index;
                                star.addFace(f2Index);
                                 e.halfedge.twin;
                        }
                
                }

                return star; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // What we must do is close each of the face, and edges in our subset to make it 
                // a subcomplex, which involed adding all vertices and edges in a boundary loop of F
                // and then for each edge adding it's two vertices.
                
                // first we create a closure subset starting with our original mesh
                let closure = MeshSubset.deepCopy(subset);
 
                // now let us loop through each edge and add the corresponding vertices
                for (let i in Array.from(subset.edges)) {
                        let eIndex = Array.from(subset.edges)[i];
                        let e = this.mesh.edges[eIndex];
                        // now we iterate through the edges adjacent to this vertex;
                        let v1Index = e.halfedge.vertex.index;
                        let v2Index = e.halfedge.next.vertex.index;
                        closure.addVertex(v1Index);
                        closure.addVertex(v2Index);
                }

                // now let us loop through each face's edges
                for (let i in Array.from(subset.faces)) {
                        let fIndex = Array.from(subset.faces)[i];
                        let f = this.mesh.faces[fIndex];
                        // now let us iterate through each edge of the face
                        for (let adjacentEdgeIndex in Array.from(f.adjacentEdges())) {
                                let adjacentEdge = Array.from(f.adjacentEdges())[adjacentEdgeIndex];
                                closure.addEdge(adjacentEdge.halfedge.edge.index);
                                //additionally we must add the vertices of this edge
                                let v1Index = adjacentEdge.halfedge.vertex.index;
                                let v2Index = adjacentEdge.halfedge.next.vertex.index;
                                closure.addVertex(v1Index);
                                closure.addVertex(v2Index);
                        }
                }

                return closure; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // Remembering that the link is just Cl(st(s)) - st(cl(s)) gives us our result
                let a1 = this.closure(this.star(subset));
                let a2 = this.star(this.closure(subset));
                a1.deleteSubset(a2);
                return a1; 
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // if a subset is a complex than cl(S) = S
                return this.closure(subset).equals(subset);
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // first check what the maximal degree of the pure complex can be
                let possibleComplex = this.isComplex(subset);
                let possibleTwoComplex = possibleComplex && (Array.from(subset.faces).length != 0);
                let possibleOneComplex = possibleComplex && (Array.from(subset.edges).length != 0);
                let possibleZeroComplex = possibleComplex && (Array.from(subset.vertices).length != 0);

                // lets make a copy of our subset to keep track of if we have a pure complex or not
                let trackingComplex = MeshSubset.deepCopy(subset);

                // now let's make a copy of our subset and delete from to keep track if we have 
                // a pure complex or not, if so deleting edges and vertices from the faces should result
                // in an empty subsetafterwards.
                if (possibleTwoComplex) {
                        // iterate through each of the faces
                        for (let i in Array.from(subset.faces)){
                                let fIndex = Array.from(subset.faces)[i];
                                // now let's get the closure of f.
                                let vSet = new Set();
                                let eSet = new Set();
                                let fSet = new Set();
                                fSet.add(fIndex);
                                let fClosure = new MeshSubset(vSet, eSet, fSet);
                                fClosure = this.closure(fClosure);
                                // now we remove fStar from trackingComplex
                                trackingComplex.deleteSubset(fClosure);
                        }
                        // now we check if our trackingComplex is empty, if so than we know it is a 2 complex
                        // make an empty complex for cmparison
                        let emptyComplex = new MeshSubset();
                        if (emptyComplex.equals(trackingComplex)) {
                                return 2;
                        } else {
                                return -1;
                        }
                } else if (!(possibleTwoComplex) && possibleOneComplex) {
                        // iterate through each of the edges
                        for (let i in Array.from(subset.edges)){
                                let eIndex = Array.from(subset.edges)[i];
                                // now let's get the closure of f.
                                let vSet = new Set();
                                let eSet = new Set();
                                let fSet = new Set();
                                eSet.add(eIndex);
                                let eClosure = new MeshSubset(vSet, eSet, fSet);
                                eClosure = this.closure(eClosure);
                                // now we remove fStar from trackingComplex
                                trackingComplex.deleteSubset(eClosure);
                        }
                        // now we check if our trackingComplex is empty, if so than we know it is a 2 complex
                        // make an empty complex for cmparison
                        let emptyComplex = new MeshSubset();
                        if (emptyComplex.equals(trackingComplex)) {
                                return 1;
                        } else {
                                return -1;
                        }
                } else if (!(possibleTwoComplex || possibleOneComplex) && possibleZeroComplex) {
                        // in this case we know that we must possess only vertices thus we can just return 0
                        return 0;
                }

                return -1;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // Here we are going to use the definition that bd(K') is the closure of the set of all simplices that are
                // proper faces of exactly one simplex of K'
                
                // lets first make the empty subset that will be our boundary
                let boundary = new MeshSubset();

                // depends on the degree of our simplex (if it is a 2-simplex we check all edges and vertices)
                if (this.isPureComplex(subset) == 2) {
                        // now we iterate through each face
                        for (let i in Array.from(subset.faces)){
                                let fIndex = Array.from(subset.faces)[i];
                                let f = this.mesh.faces[fIndex];
                                console.log(f.adjacentEdges());
                                for (let i in Array.from(f.adjacentEdges())) {
                                        let adjacentEdgeIndex = Array.from(f.adjacentEdges())[i];
                                        boundary.addEdge(adjacentEdgeIndex);
                                }
                                
                        }
                } else if (this.isPureComplex(subset) == 1) {

                } else if (this.isPureComplex(subset) == 0) {

                }

                return boundary; 
        }
}
