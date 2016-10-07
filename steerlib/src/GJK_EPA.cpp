#include "obstacles/GJK_EPA.h"
#include "obstacles/triangulate.h"


SteerLib::GJK_EPA::GJK_EPA()
{
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	//both shapes are convex
	if (isConvex(_shapeA) && isConvex(_shapeB)){
		//std::cout << "both are convex" << std::endl;
		return intersectConvex(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB);
	}

	//at least one shape is concave
	return intersectConcave(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB);
}

//check if the shape is convex
bool SteerLib::GJK_EPA::isConvex(const std::vector<Util::Vector>& shape){

	//std::cout << "new shape -------------------------" << std::endl;

	int size = shape.size();
	if (shape.size() < 4){
		return true;
	}
	bool sign = false;
	for (int i = 0; i < size; i++){
		const Util::Vector &a = shape[i];
		const Util::Vector &b = shape[(i + 1) % size];
		const Util::Vector &c = shape[(i + 2) % size];
		/*
		std::cout << "a.x = " << a.x << std::endl;
		std::cout << "a.y = " << a.y << std::endl;
		std::cout << "a.z = " << a.z << std::endl;
		std::cout << "b.x = " << b.x << std::endl;
		std::cout << "b.y = " << b.y << std::endl;
		std::cout << "b.z = " << b.z << std::endl;
		std::cout << "c.x = " << c.x << std::endl;
		std::cout << "c.y = " << c.y << std::endl;
		std::cout << "c.z = " << c.z << std::endl;
		*/
		const float dx1 = a.x - b.x;
		const float dz1 = a.z - b.z;
		const float dx2 = c.x - b.x;
		const float dz2 = c.z - b.z;

const float zcrossproduct = (dx1*dz2) - (dx2*dz1);
if (i == 0){
	//set the initial sign
	sign = (zcrossproduct < 0);
}
else{
	if (sign != (zcrossproduct < 0)){
		return false;
	}
}
	}
	return true;
}

bool SteerLib::GJK_EPA::intersectConvex(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> simplex;
	if (gjk(_shapeA, _shapeB, simplex)){
		//std::cerr << "gjk is true !!!!!!!!!!!!" << std::endl;
		epa(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplex);
		return true;
	}
	return false;
}

bool SteerLib::GJK_EPA::gjk(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	Util::Vector centerA = getCenter(_shapeA);
	Util::Vector centerB = getCenter(_shapeB);
	Util::Vector direction = centerB - centerA;
	Util::Vector PointA, PointB, PointC, d;
	Util::Vector AB, AO;
	Util::Vector origin(0, 0, 0);

	simplex.push_back(getSimplexPoint(_shapeA, _shapeB, direction));
	direction = -direction;

	PointA = simplex[0];
	PointB = simplex[1];

	AB = PointB - PointA;
	AO = origin - PointA;

	while (true){
		simplex.push_back(getSimplexPoint(_shapeA, _shapeB, direction));
		//make sure that the last point of simplex points passes the origin, else return false, 
		//since it is on the edge of Minkowski Difference and will not contain the origin anymore.
		if (dot(simplex.back(), direction) <= 0){
			return false;
		}
		else{
			if (simplexContainOrigin(simplex, direction)){
				return true;
			}
			else{
				simplex.push_back(getSimplexPoint(_shapeA, _shapeB, direction));
			}
		}
	}

	
	return false;
}

Util::Vector SteerLib::GJK_EPA::getCenter(const std::vector<Util::Vector>& shape){

	Util::Vector center(0, 0, 0);
	for (int i = 0; i < shape.size(); i++){
		center[0] += shape[i][0];
		center[2] += shape[i][2];
	}
	center[0] = center[0] / shape.size();
	center[2] = center[2] / shape.size();

	return center;
}

Util::Vector SteerLib::GJK_EPA::getSimplexPoint(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, Util::Vector& direction){
	
	Util::Vector farthestPointA, farthestPointB;

	//get points on the edge of the shapes in opposite directions
	farthestPointA = getFarthestPoint(_shapeA, direction);
	farthestPointB = getFarthestPoint(_shapeB, -direction);

	return farthestPointA - farthestPointB;
}

Util::Vector SteerLib::GJK_EPA::getFarthestPoint(const std::vector<Util::Vector>& shape, Util::Vector& direction){

	Util::Vector p = shape[0];
	float distance = dot(shape[0], direction);
	for (int i = 1; i < shape.size(); i++){
		if (dot(shape[i], direction) > distance){
			distance = dot(shape[i], direction);
			p[0] = shape[i][0];
			p[2] = shape[i][2];
		}
	}
	return p;
}

bool SteerLib::GJK_EPA::simplexContainOrigin(std::vector<Util::Vector>& simplex, Util::Vector& direction){

	Util::Vector A, B, C, AO, AB, AC, abPerp, acPerp;
	Util::Vector origin(0, 0, 0);
	A = simplex.back();
	B = simplex[0];
	AO = origin - A;
	AB = B - A;

	if (simplex.size() == 2){
		//only 2 points, it is a line
		if (isOriginOntheLine(A, B)){
			return true;
		}
		else{
			//set a new direction to be triple product of (AB,AO,AB), it is the perpendicular to AB in the direction of the origin
			direction = AO*(dot(AB, AB)) - AB*(dot(AB, AO));
		}
	}
	else{
		//3 points
		C = simplex[1];
		AC = C - A;
		//get normals
		abPerp = AB*(dot(AB, AC)) - AC*(dot(AB, AB));
		acPerp = AC*(dot(AC, AB)) - AB*(dot(AC, AC));
		//check if origin is outside of edge AB
		if (dot(abPerp, AO) > 0){
			//then C is useless, remove it
			simplex.erase(simplex.begin()+1);
			//set the new direction to be abPerp
			direction = abPerp;
		}
		else{
			if (dot(acPerp, AO) > 0){
				//check if origin is outside of edge AC
				//B is useless, remove it
				simplex.erase(simplex.begin());
				//set new direction to be acPerp
				direction = acPerp;
			}
			else{
				return true;
			}
		}
	}
	return false;
}

bool SteerLib::GJK_EPA::isOriginOnSimplex(std::vector<Util::Vector>& simplex, float& return_penetration_depth, Util::Vector& return_penetration_vector){

	for (int i = 0; i < simplex.size()-1; i++){
		for (int j = i + 1; j < simplex.size(); j++){
			//for every pair of simplex points, check if the origin is on the edge
			if (isOriginOntheLine(simplex[i], simplex[j])){
				//std::cerr << "it is on the line between point" << simplex[i] << std::endl;
				//std::cerr << "and" << simplex[j] << std::endl;
				updatePenetration(simplex[i], simplex[j], return_penetration_depth, return_penetration_vector);
				return true;
			}
		}
	}
	return false;
}

bool SteerLib::GJK_EPA::isOriginOntheLine(Util::Vector& PointA, Util::Vector& PointB){

	if (PointA.x == 0){
		if (PointB.x != 0){
			return false;
		}
		else{
			if (PointA.z * PointB.z > 0){
				return false;
			}
			else
			{
				return true;
			}
		}
	}

	//it is a vertical line, x must be 0
	if (PointA.x == PointB.x){
		if (PointA.x != 0){
			return false;
		}
	}
	
	if (PointA.z == 0){
		if (PointB.z != 0){
			return false;
		}
		else{
			if (PointA.x * PointB.x > 0){
				return false;
			}
			else{
				return true;
			}
		}
	}
	
	if (PointA.z == PointB.z){
		//horizontal line, z must be zero
		if (PointA.z != 0){
			return false;
		}
	}

	float slope = (PointB.z - PointA.z) / (PointB.x - PointA.x);

	if (((PointB.z / PointB.x) == slope) && ((PointA.z / PointA.x) == slope)){
		if (PointA.x*PointA.z*PointB.x*PointB.z > 0){
			return true;
		}
	}
	return false;
}

void SteerLib::GJK_EPA::updatePenetration(Util::Vector& PointA, Util::Vector& PointB, float& return_penetration_depth, Util::Vector& return_penetration_vector){
	if (PointA.lengthSquared() <= PointB.lengthSquared()){
		return_penetration_depth = PointA.length();
		return_penetration_vector = PointA;
	}
	else{
		return_penetration_depth = PointB.length();
		return_penetration_vector = PointB;
	}
	return_penetration_vector = return_penetration_vector / return_penetration_vector.norm();
}

void SteerLib::GJK_EPA::findClosestEdge(std::vector<Util::Vector>& simplex, Util::Vector& normal, int& index, float& distance){
	//std::cerr << "start of findCloestEdge function ----------- = " << std::endl;
	//std::cerr << "simplex.size = " << simplex.size() << std::endl;
	Util::Vector a, b, e, n, n_norm;
	double d;
	int j = 0;
	for (int i = 0; i < simplex.size(); i++){
		//if i is the last point, then take j as the first point
		if (i == simplex.size() - 1){
			j = 0;
		}
		else{
			j = i + 1;
		}
		//get two points and the edge linked them
		a = simplex[i];
		b = simplex[j];
		e = b - a;

		//std::cerr << "e = " << e << std::endl;

		//get the vector from the edge towards the origin
		//a is originToA, since origin = (0 ,0 ,0)
		//this is a tripleProduct of (e, originToA, e)
		n = a*dot(e, e) - e*dot(e, a);

		//std::cerr << "n = " << n << std::endl;

		//normalize the vector
		n_norm = n / n.norm();
		//std::cerr << "n_norm = " << n_norm << std::endl;
		//get the distance from the origin to the edge
		d = dot(n_norm, a);
		//std::cerr << "d = " << d << std::endl;

		//check to get the least distance and update all the attributes
		if (d < distance){
			distance = d;
			normal = n_norm;
			index = j;
			
		}
	}
}

void SteerLib::GJK_EPA::epa(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	Util::Vector normal;
	int index;
	

	if (isOriginOnSimplex(simplex, return_penetration_depth, return_penetration_vector)){
		//if origin is on the edge of any pair of simplex points, update p_depth as distance from the closer point to origin
		//and update p_vector as this point
		//std::cerr << "origin is on the simplex" << std::endl;
		return;
	}
	else{
		while (true){
			float distance = FLT_MAX;

			//at first we have 3 points in simplex, find the closest edge of it.
			findClosestEdge(simplex, normal, index, distance);
			//std::cerr << "normal = " << normal << std::endl;
			//std::cerr << "distance = " << distance << std::endl;
			//std::cerr << "index = " << index << std::endl;
			//get a new support point in the new direction of normal
			Util::Vector supportP = getSimplexPoint(_shapeA, _shapeB, normal);
			//std::cerr << "p = " << supportP << std::endl;
			double d = dot(supportP, normal);
			//std::cerr << "d = " << d << std::endl;

			//if the distance is close enough, return
			//else insert the support point to the index we calculate in the simplex points and loop again
			if (d - distance < 0.00001){
				return_penetration_depth = d;
				return_penetration_vector = normal;
				return;
			}
			else{
				simplex.insert(simplex.begin() + index, supportP);
			}
		}
	}
	
	
}



bool SteerLib::GJK_EPA::intersectConcave(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> trianglesA, trianglesB;
	triangulatePolygon(_shapeA, trianglesA);
	triangulatePolygon(_shapeB, trianglesB);
	float maxDepth = 0;
	Util::Vector maxVector;

	for (size_t i = 0; i < trianglesA.size(); i += 3){
		std::vector<Util::Vector> triangleA(trianglesA.begin() + i, trianglesA.begin() + i + 3);
		for (size_t j = 0; j < trianglesB.size(); j += 3){
			std::vector<Util::Vector> triangleB(trianglesB.begin() + j, trianglesB.begin() + j + 3);
			if (intersectConvex(return_penetration_depth, return_penetration_vector, triangleA, triangleB)){
				if (return_penetration_depth > maxDepth){
					maxDepth = return_penetration_depth;
					maxVector = return_penetration_vector;
					
				}
			}
		}
	}

	if (maxDepth > 0){
		return_penetration_depth = maxDepth;
		return_penetration_vector = maxVector;
		return true;
	}

	return false;
}

bool SteerLib::GJK_EPA::triangulatePolygon(const std::vector<Util::Vector>& shape, std::vector<Util::Vector>& triangles) {
	Vector2dVector shape2D;
	Vector2dVector triangles2D;
	for (const Util::Vector & vec : shape) {
		Vector2d vec2D(vec.x, vec.z);
		shape2D.push_back(vec2D);
	}
	if (Triangulate::Process(shape2D, triangles2D)) {
		for (const Vector2d & vec : triangles2D) {
			Util::Vector vec2(vec.GetX(), 0, vec.GetY());
			triangles.push_back(vec2);
		}
		return true;
	}
	return false;
}