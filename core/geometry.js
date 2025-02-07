"use strict";

class Geometry {
	/**
	 * This class represents the geometry of a {@link module:Core.Mesh Mesh}. This includes information such
	 * as the position of vertices as well as methods to compute edge lengths, corner
	 * angles, face area, normals, discrete curvatures etc.
	 * @constructor module:Core.Geometry
	 * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
	 * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is true.
	 * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @property {Object} positions A dictionary mapping each vertex to a normalized position.
	 */
	constructor(mesh, positions, normalizePositions = true) {
		this.mesh = mesh;
		this.positions = {};
		for (let i = 0; i < positions.length; i++) {
			let v = this.mesh.vertices[i];
			let p = positions[i];

			this.positions[v] = p;
		}

		if (normalizePositions) {
			normalize(this.positions, mesh.vertices);
		}
	}

	/**
	 * Computes the vector along a halfedge.
	 * @method module:Core.Geometry#vector
	 * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vector(h) {
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];

		return b.minus(a);
	}

	/**
	 * Computes the length of an edge.
	 * @method module:Core.Geometry#length
	 * @param {module:Core.Edge} e The edge whose length needs to be computed.
	 * @returns {number}
	 */
	length(e) {
		return this.vector(e.halfedge).norm();
	}

	/**
	 * Computes the midpoint of an edge.
	 * @method module:Core.Geometry#midpoint
	 * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
	 * @returns {number}
	 */
	midpoint(e) {
		let h = e.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.twin.vertex];

		return (a.plus(b)).over(2);
	}

	/**
	 * Computes the mean edge length of all the edges in a mesh.
	 * @method module:Core.Geometry#meanEdgeLength
	 * @returns {number}
	 */
	meanEdgeLength() {
		let sum = 0;
		let edges = this.mesh.edges;
		for (let e of edges) {
			sum += this.length(e);
		}

		return sum / edges.length;
	}

	/**
	 * Computes the area of a face.
	 * @method module:Core.Geometry#area
	 * @param {module:Core.Face} f The face whose area needs to be computed.
	 * @returns {number}
	 */
	area(f) {
		if (f.isBoundaryLoop()) return 0.0;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return 0.5 * u.cross(v).norm();
	}

	/**
	 * Computes the total surface area of a mesh.
	 * @method module:Core.Geometry#totalArea
	 * @returns {number}
	 */
	totalArea() {
		let sum = 0.0;
		for (let f of this.mesh.faces) {
			sum += this.area(f);
		}

		return sum;
	}

	/**
	 * Computes the normal of a face.
	 * @method module:Core.Geometry#faceNormal
	 * @param {module:Core.Face} f The face whose normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	faceNormal(f) {
		if (f.isBoundaryLoop()) return undefined;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return u.cross(v).unit();
	}

	/**
	 * Computes the centroid of a face.
	 * @method module:Core.Geometry#centroid
	 * @param {module:Core.Face} f The face whose centroid needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	centroid(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		return a.plus(b).plus(c).over(3);
	}

	/**
	 * Computes the circumcenter of a face.
	 * @method module:Core.Geometry#circumcenter
	 * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	circumcenter(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		let ac = c.minus(a);
		let ab = b.minus(a);
		let w = ab.cross(ac);

		let u = (w.cross(ab)).times(ac.norm2());
		let v = (ac.cross(w)).times(ab.norm2());
		let x = (u.plus(v)).over(2 * w.norm2());

		return x.plus(a);
	}

	/**
	 * Computes an orthonormal bases for a face.
	 * @method module:Core.Geometry#orthonormalBases
	 * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
	 * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
	 */
	orthonormalBases(f) {
		let e1 = this.vector(f.halfedge).unit();

		let normal = this.faceNormal(f);
		let e2 = normal.cross(e1);

		return [e1, e2];
	}

	/**
	 * Computes the angle (in radians) at a corner.
	 * @method module:Core.Geometry#angle
	 * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
	 * @returns {number} The angle clamped between 0 and π.
	 */
	angle(c) {
		// first let us get the half edge associated with this angle 
		let adjacent = this.vector(c.halfedge.next).unit();
		let hypotenuse = this.vector(c.halfedge.next.next).unit().negated();

		// now we use the dot product formula to calculate the angle between them
		return Math.acos(adjacent.dot(hypotenuse));
	}

	/**
	 * Computes the cotangent of the angle opposite to a halfedge.
	 * @method module:Core.Geometry#cotan
	 * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
	 * @returns {number}
	 */
	cotan(h) {
		// first we get the length of this opposite half edge
		let v = this.vector(h.next.twin);
		// then we must get the length of adjacent halfedge
		let u = this.vector(h.next.next);
		// get this length using the cross product formula
		let cotan = u.dot(v)/(u.cross(v).norm());


		return cotan; // placeholder
	}

	/**
	 * Computes the signed angle (in radians) between two adjacent faces.
	 * @method module:Core.Geometry#dihedralAngle
	 * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
	 * the dihedral angle is computed.
	 * @returns {number} The dihedral angle.
	 */
	dihedralAngle(h) {
		// get the normal of the face associated with this halfedge
		let face_normal_1 = this.faceNormal(h.face);
		// now get the edge and length corresponding to the halfedge
		let edge = this.vector(h).unit();
		//get the normal of the adjacent face to this edge
		let face_normal_2 = this.faceNormal(h.twin.face);

		// now we use the angle formula given the following
		let theta = Math.atan2(edge.dot(face_normal_1.cross(face_normal_2)), face_normal_1.dot(face_normal_2));

		return theta; // placeholder
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	barycentricDualArea(v) {
		// to compute the barycentric dual we iterate through all faces adjacent to this vertex
		// and sum up their areas, and afterwards divide by 3
		let area = 0;
		for (let f of v.adjacentFaces()){
			area += this.area(f);
		}

		// divide by 3 to get the barycentric dual area
		return (1/3)*area; // placeholder
	}

	/**
	 * Computes the circumcentric dual area of a vertex.
	 * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
	 * @method module:Core.Geometry#circumcentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
	 * @returns {number}
	 */
	circumcentricDualArea(v) {
		let a = 0;
		for (let c of v.adjacentCorners()) {
			let hedge1 = c.halfedge.next;
			let hedge2 = c.halfedge.next;

			let edge1 = this.vector(hedge1);
			let edge2 = this.vector(hedge2);

			let alpha = this.cotan(hedge1);
			let beta = this.cotan(hedge2);

			a += edge1.norm2()*alpha + edge2.norm2()*beta;
		}

		return a/8; 
	}

	/**
	 * Computes the normal at a vertex using the "equally weighted" method.
	 * @method module:Core.Geometry#vertexNormalEquallyWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalEquallyWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);

			n.incrementBy(normal);
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "face area weights" method.
	 * @method module:Core.Geometry#vertexNormalAreaWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAreaWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);
			let face_area = this.area(f);

			// here we weight each normal by the area of its associated triangular face
			// for some reason the scaling operation isn't working so I'm just going to manually scale the normal vectors
			normal.x *= face_area;
			normal.y *= face_area;
			normal.z *= face_area;

			// now we increment our normal vector which is the sum of the face weighted areas
			n.incrementBy(normal);
		}

		// normalize and return
		n.normalize();

		return n; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "tip angle weights" method.
	 * @method module:Core.Geometry#vertexNormalAngleWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAngleWeighted(v) {
		// what we will do first is loop over the faces containing the vertex i 
		// it seems like the problem ehre is with the angle functions output
		let n = new Vector();
		for (let corner of v.adjacentCorners()) {
			let angle = this.angle(corner);
			let normal = this.faceNormal(corner.halfedge.face);

			//now we will scale the normal with the associated angle
			normal.x *= angle;
			normal.y *= angle;
			normal.z *= angle;

			n.incrementBy(normal);
		}

		n.normalize();

		return n; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "gauss curvature" method.
	 * @method module:Core.Geometry#vertexNormalGaussCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalGaussCurvature(v) {
		// it seems like the problem here is with the dihedral angle function (prob )
		let n = new Vector();
		// it seems like the problem here maybe over the edges being iterated (i.e) orentation of the edge vector??
		for (let e of v.adjacentEdges()) {
			let angle = this.dihedralAngle(e.halfedge);
			let edge = this.vector(e.halfedge).unit();

			let scaling_factor = angle;
			edge.x *= scaling_factor;
			edge.y *= scaling_factor;
			edge.z *= scaling_factor;

			n.incrementBy(edge);
		}

		//now scale it by 1/2
		n.x *= (0.5);
		n.y *= (0.5);
		n.z *= (0.5);

		n.normalize();
		return n; 
	}

	/**
	 * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
	 * @method module:Core.Geometry#vertexNormalMeanCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalMeanCurvature(v) {
		let n = new Vector();
		for (let e of v.adjacentEdges()) {
			let edge = this.vector(e.halfedge);
			let cotan_1 = this.cotan(e.halfedge);
			let cotan_2 = this.cotan(e.halfedge.twin);

			let scaling_factor = cotan_1 + cotan_2;
			edge.x *= scaling_factor;
			edge.y *= scaling_factor;
			edge.z *= scaling_factor;

			n.incrementBy(edge);
		}

		//now scale it by 1/2
		n.x *= (0.5);
		n.y *= (0.5);
		n.z *= (0.5);

		n.normalize();

		return n; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "inscribed sphere" method.
	 * @method module:Core.Geometry#vertexNormalSphereInscribed
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalSphereInscribed(v) {
		let n = new Vector();
		for (let c of v.adjacentCorners()) {
			let edge_1 = this.vector(c.halfedge.next);
			let edge_2 = this.vector(c.halfedge.next.next);
			
			// now swap the orientation of one of them so they're the same
			edge_2.negated();

			// now let's get the scaling factor for the new "normal"
			let scaling_factor = 1/(edge_1.norm()*edge_1.norm()*edge_2.norm()*edge_2.norm());
			// now let's get the vector associated with this new "normal"
			let normal = edge_1.cross(edge_2);

			// since scaling isn't working we scale normal by scaling_factor
			normal.x *= scaling_factor;
			normal.y *= scaling_factor;
			normal.z *= scaling_factor;

			// finally we add this to our n vector
			n.incrementBy(normal);
		}

		//normalize n
		n.normalize();

		return n; 
	}

	/**
	 * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
	 * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
	 * @method module:Core.Geometry#angleDefect
	 * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
	 * @returns {number}
	 */
	angleDefect(v) {
		let angle_defect = 2*Math.PI;

		//Now lets iterate through the corners adjacent to v and calculate the angle!
		for (let c of v.adjacentCorners()) {
			let angle = this.angle(c);
			angle_defect -= angle;
		}

		return angle_defect; // placeholder
	}

	/**
	 * Computes the (integrated) scalar gauss curvature at a vertex.
	 * @method module:Core.Geometry#scalarGaussCurvature
	 * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
	 * @returns {number}
	 */
	scalarGaussCurvature(v) {
		return this.angleDefect(v);
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	scalarMeanCurvature(v) {
		let mean_curvature = 0;
		for (let e of v.adjacentEdges()) {
			let edge_length = this.length(e);
			let dihedral_angle = this.dihedralAngle(e.halfedge);
			mean_curvature += edge_length*dihedral_angle;
		}

		// finally scale the mean curvature by 1/2
		mean_curvature = mean_curvature/2;

		return mean_curvature;
	}

	/**
	 * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
	 * @method module:Core.Geometry#totalAngleDefect
	 * @returns {number}
	 */
	totalAngleDefect() {
		// utilize the Gauss-Bonnet Theorem
		return 2*Math.PI*this.mesh.eulerCharacteristic();
	}

	/**
	 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
	 * @method module:Core.Geometry#principalCurvatures
	 * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
	 * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
	 */
	principalCurvatures(v) {
		let K = this.scalarGaussCurvature(v);
		let H = this.scalarMeanCurvature(v);
		let k1 = Math.max(H + Math.sqrt(H*H - K), H - Math.sqrt(H*H - K));
		let k2 = Math.min(H + Math.sqrt(H*H - K), H - Math.sqrt(H*H - K));

		return [k2, k1]; // placeholder
	}

	/**
	 * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#laplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	laplaceMatrix(vertexIndex) {
		// TODO

		return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
	 * of a mesh.
	 * @method module:Core.Geometry#massMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	massMatrix(vertexIndex) {
		// TODO

		return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#complexLaplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	complexLaplaceMatrix(vertexIndex) {
		// TODO

		return ComplexSparseMatrix.identity(1, 1); // placeholder
	}
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
	// compute center of mass
	let N = vertices.length;
	let cm = new Vector();
	for (let v of vertices) {
		let p = positions[v];

		cm.incrementBy(p);
	}
	cm.divideBy(N);

	// translate to origin and determine radius
	let radius = -1;
	for (let v of vertices) {
		let p = positions[v];

		p.decrementBy(cm);
		radius = Math.max(radius, p.norm());
	}

	// rescale to unit radius
	if (rescale) {
		for (let v of vertices) {
			let p = positions[v];

			p.divideBy(radius);
		}
	}
}