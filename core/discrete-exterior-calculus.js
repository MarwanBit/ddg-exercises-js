"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// because each vertex if assumed to have area 1 this is merely a diagonal matrix each each entry
		// *_i,i being the area of the barycentric dual of that particular vertex

		// what we do here is merely create a |V|*|V| size matrix and make the entry 
		// first let's get the size of the dictionary
		let vSize = Object.keys(vertexIndex).length;
		// now we make our Hodge Star operator!
		let hodgeStar = new Triplet(vSize, vSize);
		// now we are going to iterate through each vertex in our vertexIndex dictionary
		for (let v in vertexIndex){
			let vIndex = vertexIndex[v];
			let vDualArea = geometry.barycentricDualArea(geometry.mesh.vertices[vIndex]);
			// add the dual area of v_i to the *_i,i entry of the hodge star
			hodgeStar.addEntry(vDualArea, vIndex, vIndex);
		}

		return SparseMatrix.fromTriplet(hodgeStar); // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// because each vertex if assumed to have area 1 this is merely a diagonal matrix each each entry
		// *_i,i being is the ratio of the dual edge length to the primal edge length l_dual/l_primal

		// what we do here is merely create a |e|*|e| size matrix and make the entry 
		// first let's get the size of the dictionary
		let eSize = Object.keys(edgeIndex).length;
		// now we make our Hodge Star operator!
		let hodgeStar = new Triplet(eSize, eSize);
		// now we are going to iterate through each edge in our edgeIndex dictionary
		for (let e in edgeIndex){
			let eIndex = edgeIndex[e];
			// now we find the cotangent of each of the angles and use the cotangent formula
			// to get l_dual / l_primal
			let edge1 = geometry.mesh.edges[eIndex].halfedge;
			let edge2 = edge1.twin;

			//now that we have both sides of the halfedge we can get both angles
			let cotan1 = geometry.cotan(edge1);
			let cotan2 = geometry.cotan(edge2);

			//now we know that the value should be 1/2 (cotan1 + cotan2)
			let lengthRatio = 0.5*(cotan1 + cotan2);
			// add the dual area of v_i to the *_i,i entry of the hodge star
			hodgeStar.addEntry(lengthRatio, eIndex, eIndex);
		}

		// now we return the matrix
		return SparseMatrix.fromTriplet(hodgeStar);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		// because each vertex if assumed to have area 1 this is merely a diagonal matrix with each entry
		// *_i,i being the 1/area(face) where face is the ith face
		// what we do here is merely create a |f|*|f| size matrix and make the entry 
		// first let's get the size of the dictionary
		let fSize = Object.keys(faceIndex).length;
		// now we make our Hodge Star operator!
		let hodgeStar = new Triplet(fSize, fSize);
		// now we are going to iterate through each face in our faceIndex dictionary
		for (let f in faceIndex){
			let fIndex = faceIndex[f];
			let face = geometry.mesh.faces[fIndex];
			let area = (1/geometry.area(face));
			// add the area 1/area(face) to the i,i'th entry of *
			hodgeStar.addEntry(area, fIndex, fIndex);
		}

		return SparseMatrix.fromTriplet(hodgeStar); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		// To build the d_{0} operator we must first create a signed incidence matrix for the exterior derivative
		// we know the size will be |E| * |V| so let's make a sparse matrix
		let eSize = Object.keys(edgeIndex).length;
		let vSize = Object.keys(vertexIndex).length;
		// now we create our exterior derivative for 0 forms which we will call d_0
		let d_0 = new Triplet(eSize, vSize);

		//now what we must do is iterate through each edge and edit the oriented adjacency matrix;
		for (let eIndex in Object.values(edgeIndex)) {
			// we must convert eIndex to an integer
			eIndex = parseInt(eIndex);
			let edge = geometry.mesh.edges[eIndex];
			// now we must get the vertex the current edge is "leaving"
			let v1 = edge.halfedge.vertex;
			// now we must get the vertex the current edge is "moving towards"
			let v2 = edge.halfedge.next.vertex;
			// now for v1 and v2 we must get their corresponding indices
			let v1Index = vertexIndex[v1];
			let v2Index = vertexIndex[v2];
			// now finally we just add these entries into our d_0 matrix
			d_0.addEntry(-1, eIndex, v1Index);
			d_0.addEntry(1, eIndex, v2Index);
			//This completes the matrix!
		}

		return SparseMatrix.fromTriplet(d_0); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		// To build the d_{1} operator we must first create a signed incidence matrix for the exterior derivative
		// we know the size will be |F| * |E| so let's make a sparse matrix
		let fSize = Object.keys(faceIndex).length;
		let eSize = Object.keys(edgeIndex).length;
		// now we create our exterior derivative for 1 forms which we will call d_1
		let d_1 = new Triplet(fSize, eSize);

		//now what we must do is iterate through each face and edit the oriented adjacency matrix;
		for (let fIndex in Object.values(faceIndex)) {
			// we must convert fIndex to an integer
			fIndex = parseInt(fIndex);
			let face = geometry.mesh.faces[fIndex];

			//now we must iterate through the edges adjacent to this face
			// and check their relative orientation to the face
			for (let adjacentEdge in Array.from(face.adjacentEdges())){
				// we check the orientation by checking it is in the faces adjacent halfedges
				adjacentEdge = Array.from(face.adjacentEdges())[adjacentEdge];

				// now we get the index of this edge
				let adjacentEdgeIndex = edgeIndex[adjacentEdge];
				let relativeOrientation = false;

				//check if the adjacentEdge is in face.adjacentHalfedges to determine relative orientation
				for (let adjacentHalfEdgeIndex in Array.from(face.adjacentHalfedges())) {
					let adjacentHalfEdge = Array.from(face.adjacentHalfedges())[adjacentHalfEdgeIndex];
					if (adjacentEdge.halfedge == adjacentHalfEdge) {
						relativeOrientation = true;
					} else if (adjacentEdge.halfedge.twin == adjacentHalfEdge) {
						relativeOrientation = false;
					}
				}
				

				// check if adjacentEdge.edge is equal to any of these
				// it seems like the error is here, we need to check orientation, one way to do this is to see
				// that for two halfedges h1, h2 that h1.next.vertex = h2.next.vertex (indicating they have the same)
				// "direction"
				if (relativeOrientation) {
					d_1.addEntry(1, fIndex, adjacentEdgeIndex);
				} else {
					d_1.addEntry(-1, fIndex, adjacentEdgeIndex);
				}
			}
		}
		return SparseMatrix.fromTriplet(d_1); // placeholder
	}
}
