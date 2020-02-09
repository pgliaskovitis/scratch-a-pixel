/*
 * Copyright (C) 2012  www.scratchapixel.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <atomic>
#include <queue>
#include <cmath>

struct AccelerationStats {
	std::atomic<uint32_t> numPrimaryRays;
	std::atomic<uint32_t> numRayMeshTests;
	std::atomic<uint32_t> numRayBBoxTests;
	std::atomic<uint32_t> numRayBoundingVolumeTests;
	std::atomic<uint32_t> numRayBoundingVolumeIntersections;
	std::atomic<uint32_t> numRayTriangleTests;
	std::atomic<uint32_t> numRayTriangleIntersections;
};

// [comment]
// The most basic acceleration class (the parent class of all the other acceleration structures)
// could have a *pure* virtual intersect() method but instead we decided in this implementation
// to have it supporting the basic ray-mesh intersection routine.
// [/comment]
class AccelerationStructure
{
public:
	// [comment]
	// We transfer owner ship of the mesh list to the acceleration structure. This makes
	// more sense from a functional/structure stand point because the objects/meshes themselves
	// should be destroyed/deleted when the acceleration structure is being deleted
	// Ideally this means the render function() itself should be bounded (in terms of lifespan)
	// to the lifespan of the acceleration structure (aka we should wrap the accel structure instance
	// and the render method() within the same object, so that when this object is deleted,
	// the render function can't be called anymore.
	// [/comment]
	AccelerationStructure(std::vector<std::unique_ptr<const TriangleMesh>>& m) : meshes(std::move(m))
	{
		stats.numPrimaryRays = 0;
		stats.numRayTriangleTests = 0;
		stats.numRayTriangleIntersections = 0;
		stats.numRayBBoxTests = 0;
		stats.numRayBoundingVolumeTests = 0;
		stats.numRayBoundingVolumeIntersections = 0;
	}
	virtual ~AccelerationStructure() {}
	virtual bool intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
	{
		// [comment]
		// Because we don't want to change the content of the mesh itself, just get a pointer to it so
		// it's safer to make it const (which doesn't mean we can't change its assignment just that
		// we can't do something like intersectedMesh->color = blue. You would get something like:
		// "read-only variable is not assignable" error message at compile time)
		// [/comment]
		const TriangleMesh* intersectedMesh = nullptr;
		float t = kInfinity;
		for (const auto& mesh: meshes) {
			stats.numRayMeshTests++;
			if (mesh->intersect(orig, dir, t,
					    stats.numRayTriangleTests,
					    stats.numRayTriangleIntersections) && t < tHit) {
				intersectedMesh = mesh.get();
				tHit = t;
			}
		}

		return (intersectedMesh != nullptr);
	}
	AccelerationStats& getStats() { return stats;}
protected:
		const std::vector<std::unique_ptr<const TriangleMesh>> meshes;
		mutable struct AccelerationStats stats;
};

// [comment]
// Implementation of the ray-bbox method. If the ray intersects the bbox of a mesh then
// we test if the ray intersects the mesh contained by the bbox itself.
// [/comment]
class BBoxAcceleration : public AccelerationStructure
{
public:
	BBoxAcceleration(std::vector<std::unique_ptr<const TriangleMesh>>& m) : AccelerationStructure(m) {}

	// [comment]
	// Implement the ray-bbox acceleration method. The method consist of intersecting the
	// ray against the bbox of the mesh first, and if the ray inteesects the boudning box
	// then test if the ray intersects the mesh itself. It is obvious that the ray can't
	// intersect the mesh if it doesn't intersect its boudning volume (a box in this case)
	// [/comment]
	virtual bool intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
	{
		const TriangleMesh* intersectedMesh = nullptr;
		const Vec3f invDir = 1 / dir;
		const Vec3b sign(dir.x < 0, dir.y < 0, dir.z < 0);
		float t = kInfinity;
		for (const auto& mesh : meshes) {
			// If you intersect the box
			stats.numRayBBoxTests++;
			if (mesh->bbox.intersect(orig, invDir, sign, t)) {
				// Then test if the ray intersects the mesh and if it does then first check
				// if the intersection distance is the nearest and if we pass that test as well
				// then update tNear variable with t and keep a pointer to the intersected mesh
				stats.numRayMeshTests++;
				if (mesh->intersect(orig, dir, t,
						    stats.numRayTriangleTests,
						    stats.numRayTriangleIntersections)) {
					intersectedMesh = mesh.get();
					tHit = t;
				}
			}
		}

		// Return true if the variable intersectedMesh is not null, false otherwise
		return (intersectedMesh != nullptr);
	}
};

// [comment]
// Implementation of the Bounding Volume Hieratchy (BVH) acceleration structure
// [/comment]
class BVH : public AccelerationStructure
{
	static const uint8_t kNumPlaneSetNormals = 7;
	static const Vec3f planeSetNormals[kNumPlaneSetNormals];
	struct Extents
	{
		Extents()
		{
			for (uint8_t i = 0;  i < kNumPlaneSetNormals; ++i)
			d[i][0] = kInfinity, d[i][1] = -kInfinity;
		}
		void extendBy(const Extents& e)
		{
			for (uint8_t i = 0;  i < kNumPlaneSetNormals; ++i) {
				if (e.d[i][0] < d[i][0]) d[i][0] = e.d[i][0];
				if (e.d[i][1] > d[i][1]) d[i][1] = e.d[i][1];
			}
		}
		/* inline */
		Vec3f centroid() const
		{
			return Vec3f(
				d[0][0] + d[0][1] * 0.5,
				d[1][0] + d[1][1] * 0.5,
				d[2][0] + d[2][1] * 0.5);
		}
		bool intersect(const float*, const float*, float&, float&, uint8_t&) const;
		float d[kNumPlaneSetNormals][2];
		AccelerationStats* stats;
		const TriangleMesh* mesh;
	};

	struct Octree
	{
		Octree(const Extents& sceneExtents)
		{
			float xDiff = sceneExtents.d[0][1] - sceneExtents.d[0][0];
			float yDiff = sceneExtents.d[1][1] - sceneExtents.d[1][0];
			float zDiff = sceneExtents.d[2][1] - sceneExtents.d[2][0];
			float maxDiff = std::max(xDiff, std::max(yDiff, zDiff));
			Vec3f minPlusMax(
			sceneExtents.d[0][0] + sceneExtents.d[0][1],
			sceneExtents.d[1][0] + sceneExtents.d[1][1],
			sceneExtents.d[2][0] + sceneExtents.d[2][1]);
			bbox[0] = (minPlusMax - maxDiff) * 0.5;
			bbox[1] = (minPlusMax + maxDiff) * 0.5;
			root = new OctreeNode;
		}

		~Octree() { deleteOctreeNode(root); }

		void insert(const Extents* extents) { insert(root, extents, bbox, 0); }
		void build() { build(root, bbox); };
		void bind_intersection_test_counters(AccelerationStats *stats)
		{
			root->stats = stats;
			root->nodeExtents.stats = stats;
		}

		struct OctreeNode
		{
			OctreeNode* child[8] = { nullptr };
			std::vector<const Extents *> nodeExtentsList; // pointer to the objects extents
			Extents nodeExtents; // extents of the octree node itself
			bool isLeaf = true;
			AccelerationStats* stats;
		};

		struct QueueElement
		{
			const OctreeNode *node; // octree node held by this element in the queue
			float t; // distance from the ray origin to the extents of the node
			QueueElement(const OctreeNode *n, float tn) : node(n), t(tn) {}
			// priority_queue behaves like a min-heap
			friend bool operator < (const QueueElement &a, const QueueElement &b) { return a.t > b.t; }
		};

		OctreeNode* root = nullptr; // make unique so we don't have to manage deallocation
		BBox<> bbox;

	private:

		void deleteOctreeNode(OctreeNode*& node)
		{
			for (uint8_t i = 0; i < 8; i++) {
				if (node->child[i] != nullptr) {
					deleteOctreeNode(node->child[i]);
				}
			}
			delete node;
		}

		void insert(OctreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth)
		{
			if (node->isLeaf) {
				if (node->nodeExtentsList.size() == 0 || depth == 16) {
					node->nodeExtentsList.push_back(extents);
				}
				else {
					node->isLeaf = false;
					// Re-insert extents held by this node
					while (node->nodeExtentsList.size()) {
						insert(node, node->nodeExtentsList.back(), bbox, depth);
						node->nodeExtentsList.pop_back();
					}
					// Insert new extent
					insert(node, extents, bbox, depth);
				}
			}
			else {
				// Need to compute in which child of the current node this extents should
				// be inserted into
				Vec3f extentsCentroid = extents->centroid();
				Vec3f nodeCentroid = (bbox[0] + bbox[1]) * 0.5;
				BBox<> childBBox;
				uint8_t childIndex = 0;
				// x-axis
				if (extentsCentroid.x > nodeCentroid.x) {
					childIndex = 4;
					childBBox[0].x = nodeCentroid.x;
					childBBox[1].x = bbox[1].x;
				}
				else {
					childBBox[0].x = bbox[0].x;
					childBBox[1].x = nodeCentroid.x;
				}
				// y-axis
				if (extentsCentroid.y > nodeCentroid.y) {
					childIndex += 2;
					childBBox[0].y = nodeCentroid.y;
					childBBox[1].y = bbox[1].y;
				}
				else {
					childBBox[0].y = bbox[0].y;
					childBBox[1].y = nodeCentroid.y;
				}
				// z-axis
				if (extentsCentroid.z > nodeCentroid.z) {
					childIndex += 1;
					childBBox[0].z = nodeCentroid.z;
					childBBox[1].z = bbox[1].z;
				}
				else {
					childBBox[0].z = bbox[0].z;
					childBBox[1].z = nodeCentroid.z;
				}

				// Create the child node if it doesn't exsit yet and then insert the extents in it
				if (node->child[childIndex] == nullptr) {
					node->child[childIndex] = new OctreeNode;
					node->child[childIndex]->stats = node->stats;
					node->child[childIndex]->nodeExtents.stats = node->stats;
				}
				insert(node->child[childIndex], extents, childBBox, depth + 1);
			}
		}

		void build(OctreeNode*& node, const BBox<>& bbox)
		{
			if (node->isLeaf) {
				for (const auto& e: node->nodeExtentsList) {
					node->nodeExtents.extendBy(*e);
				}
			}
			else {
				for (uint8_t i = 0; i < 8; ++i) {
					if (node->child[i]) {
						BBox<> childBBox;
						Vec3f centroid = bbox.centroid();
						// x-axis
						childBBox[0].x = (i & 4) ? centroid.x : bbox[0].x;
						childBBox[1].x = (i & 4) ? bbox[1].x : centroid.x;
						// y-axis
						childBBox[0].y = (i & 2) ? centroid.y : bbox[0].y;
						childBBox[1].y = (i & 2) ? bbox[1].y : centroid.y;
						// z-axis
						childBBox[0].z = (i & 1) ? centroid.z : bbox[0].z;
						childBBox[1].z = (i & 1) ? bbox[1].z : centroid.z;

						// Inspect child
						build(node->child[i], childBBox);

						// Expand extents with extents of child
						node->nodeExtents.extendBy(node->child[i]->nodeExtents);
					}
				}
			}
		}
	};

	std::vector<Extents> extentsList;
	Octree* octree = nullptr;
public:
	BVH(std::vector<std::unique_ptr<const TriangleMesh>>& m);
	bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&) const;
	~BVH() { delete octree; }
};

const Vec3f BVH::planeSetNormals[BVH::kNumPlaneSetNormals] = {
	Vec3f(1, 0, 0),
	Vec3f(0, 1, 0),
	Vec3f(0, 0, 1),
	Vec3f( sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	Vec3f(-sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	Vec3f(-sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	Vec3f( sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f)
};

BVH::BVH(std::vector<std::unique_ptr<const TriangleMesh>>& m) : AccelerationStructure(m)
{
	Extents sceneExtents; // that's the extent of the entire scene which we need to compute for the octree
	extentsList.reserve(meshes.size());
	for (uint32_t i = 0; i < meshes.size(); ++i) {
		for (uint8_t j = 0; j < kNumPlaneSetNormals; ++j) {
			for (uint32_t k = 0; k < meshes[i]->numTris * 3; ++k) {
				float d = planeSetNormals[j].dotProduct(meshes[i]->P[k]);
				// set dNear and dFar
				if (d < extentsList[i].d[j][0]) extentsList[i].d[j][0] = d;
				if (d > extentsList[i].d[j][1]) extentsList[i].d[j][1] = d;
			}
		}
		sceneExtents.extendBy(extentsList[i]); // expand the scene extent of this object's extent
		extentsList[i].mesh = meshes[i].get(); // the extent itself needs to keep a pointer to the object its holds
		extentsList[i].stats = &this->stats;
	}

	// Now that we have the extent of the scene we can start building our octree
	// Using C++ make_unique function here but you don't need to, just to learn something...
	octree = new Octree(sceneExtents);
	octree->bind_intersection_test_counters(&this->stats);

	for (uint32_t i = 0; i < meshes.size(); ++i) {
		octree->insert(&extentsList[i]);
	}

	// Build from bottom up
	octree->build();
}

bool BVH::Extents::intersect(
	const float* precomputedNumerator,
	const float* precomputedDenominator,
	float& tNear,   // tn and tf in this method need to be contained
	float& tFar,    // within the range [tNear:tFar]
	uint8_t& planeIndex) const
{
	this->stats->numRayBoundingVolumeTests++;
	for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
		float tNearExtents = (d[i][0] - precomputedNumerator[i]) / precomputedDenominator[i];
		float tFarExtents = (d[i][1] - precomputedNumerator[i]) / precomputedDenominator[i];
		if (precomputedDenominator[i] < 0) std::swap(tNearExtents, tFarExtents);
		if (tNearExtents > tNear) tNear = tNearExtents, planeIndex = i;
		if (tFarExtents < tFar) tFar = tFarExtents;
		if (tNear > tFar) return false;
	}
	this->stats->numRayBoundingVolumeIntersections++;
	return true;
}

bool BVH::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
{
	tHit = kInfinity;
	const TriangleMesh* intersectedMesh = nullptr;
	float precomputedNumerator[BVH::kNumPlaneSetNormals];
	float precomputedDenominator[BVH::kNumPlaneSetNormals];
	for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
		precomputedNumerator[i] = planeSetNormals[i].dotProduct(orig);
		precomputedDenominator[i] = planeSetNormals[i].dotProduct(dir);
	}

	uint8_t planeIndex;
	float tNear = 0, tFar = kInfinity; // tNear, tFar for the intersected extents
	if (!octree->root->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNear, tFar, planeIndex) || tFar < 0)
		return false;
	tHit = tFar;
	std::priority_queue<BVH::Octree::QueueElement> queue;
	queue.push(BVH::Octree::QueueElement(octree->root, 0));
	while (!queue.empty() && queue.top().t < tHit) {
		const Octree::OctreeNode *node = queue.top().node;
		queue.pop();
		if (node->isLeaf) {
			for (const auto& e: node->nodeExtentsList) {
				float t = kInfinity;
				if (e->mesh->intersect(orig, dir, t) && t < tHit) {
					tHit = t;
					intersectedMesh = e->mesh;
				}
			}
		}
		else {
			for (uint8_t i = 0; i < 8; ++i) {
				if (node->child[i] != nullptr) {
					float tNearChild = 0, tFarChild = tFar;
					if (node->child[i]->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNearChild, tFarChild, planeIndex)) {
						float t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
						queue.push(BVH::Octree::QueueElement(node->child[i], t));
					}
				}
			}
		}
	}

	return (intersectedMesh != nullptr);
}

// [comment]
// Implementation of the Grid acceleration structure
// [/comment]
class Grid : public AccelerationStructure
{
	struct Cell
	{
		Cell() {}
		struct TriangleDesc
		{
			TriangleDesc(const TriangleMesh* m, const uint32_t &t) : mesh(m), tri(t) {}
			const TriangleMesh* mesh;
			uint32_t tri;
		};

		void insert(const TriangleMesh* mesh, uint32_t t)
		{ triangles.push_back(Grid::Cell::TriangleDesc(mesh, t)); }

		bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&, const TriangleMesh*&) const;

		std::vector<TriangleDesc> triangles;
		AccelerationStats* stats;
	};
public:
	Grid(std::vector<std::unique_ptr<const TriangleMesh>>& m);
	~Grid()
	{
		for (uint32_t i = 0; i < resolution[0] * resolution[1] * resolution[2]; ++i) {
			if (cells[i] != NULL) {
				delete cells[i];
			}
		}
		delete [] cells;
	}
	bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&) const;
	Cell **cells;
	BBox<> bbox;
	Vec3ui resolution;
	Vec3f cellDimension;
};

Grid::Grid(std::vector<std::unique_ptr<const TriangleMesh>>& m) : AccelerationStructure(m)
{
	uint32_t totalNumTriangles = 0;
	for (const auto& m : meshes) {
		bbox.extendBy(m->bbox[0]);
		bbox.extendBy(m->bbox[1]);
		totalNumTriangles += m->numTris;
	}
	// Create the grid
	Vec3f size = bbox[1] - bbox[0];
	float cubeRoot = powf(totalNumTriangles / (size.x * size.y * size.z), 1. / 3.f);
	for (uint8_t i = 0; i < 3; ++i) {
		resolution[i] = std::floor(size[i] * cubeRoot);
		if (resolution[i] < 1) resolution[i] = 1;
		if (resolution[i] > 128) resolution[i] = 128;
		cellDimension[i] = size[i] / resolution[i];
	}

	// [comment]
	// Allocate memory - note that we don't create the cells yet at this
	// point but just an array of pointers to cell. We will create the cells
	// dynamically later when we are sure to insert something in them
	// [/comment]
	uint32_t numCells = resolution.x * resolution.y * resolution.z;
	cells = new Grid::Cell* [numCells];
	memset(cells, 0x0, sizeof(Grid::Cell*) * numCells);

	for (const auto& m : meshes) {
		for (uint32_t i = 0, off = 0; i < m->numTris; ++i, off += 3) {
			Vec3f min(kInfinity), max(-kInfinity);
			const Vec3f& v0 = m->P[m->trisIndex[off]];
			const Vec3f& v1 = m->P[m->trisIndex[off + 1]];
			const Vec3f& v2 = m->P[m->trisIndex[off + 2]];
			for (uint8_t j = 0; j < 3; ++j) {
				if (v0[j] < min[j]) min[j] = v0[j];
				if (v1[j] < min[j]) min[j] = v1[j];
				if (v2[j] < min[j]) min[j] = v2[j];
				if (v0[j] > max[j]) max[j] = v0[j];
				if (v1[j] > max[j]) max[j] = v1[j];
				if (v2[j] > max[j]) max[j] = v2[j];
			}
			// Convert to cell coordinates
			min = min - bbox[0];
			max = max - bbox[0];
			for (uint8_t j = 0; j < 3; ++j) {
				min[j] = min[j] / cellDimension[j];
				max[j] = max[j] / cellDimension[j];
			}
			uint32_t zmin = scratch::utils::clamp<uint32_t>(std::floor(min[2]), 0, resolution[2] - 1);
			uint32_t zmax = scratch::utils::clamp<uint32_t>(std::floor(max[2]), 0, resolution[2] - 1);
			uint32_t ymin = scratch::utils::clamp<uint32_t>(std::floor(min[1]), 0, resolution[1] - 1);
			uint32_t ymax = scratch::utils::clamp<uint32_t>(std::floor(max[1]), 0, resolution[1] - 1);
			uint32_t xmin = scratch::utils::clamp<uint32_t>(std::floor(min[0]), 0, resolution[0] - 1);
			uint32_t xmax = scratch::utils::clamp<uint32_t>(std::floor(max[0]), 0, resolution[0] - 1);
			// Loop over the cells the triangle overlaps and insert
			for (uint32_t z = zmin; z <= zmax; ++z) {
				for (uint32_t y = ymin; y <= ymax; ++y) {
					for (uint32_t x = xmin; x <= xmax; ++x) {
						uint32_t index = z * resolution[0] * resolution[1] + y * resolution[0] + x;
						if (cells[index] == NULL) cells[index] = new Grid::Cell;
						cells[index]->stats = &this->stats;
						cells[index]->insert(m.get(), i);
					}
				}
			}
		}
	}
}

bool Grid::Cell::intersect(
	const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId,
	float& tHit, const TriangleMesh*& intersectedMesh) const
{
	bool intersected = false;
	uint32_t intersectedTri;
	float uhit, vhit;
	for (uint32_t i = 0; i < triangles.size(); ++i) {
		// [comment]
		// Be sure that rayId is never 0 - because all mailbox values
		// in the array are initialized with 0 too
		// [/comment]
		if (rayId != triangles[i].mesh->mailbox[triangles[i].tri]) {
			const TriangleMesh *mesh = triangles[i].mesh;
			uint32_t j = triangles[i].tri * 3;
			const Vec3f &v0 = mesh->P[mesh->trisIndex[j    ]];
			const Vec3f &v1 = mesh->P[mesh->trisIndex[j + 1]];
			const Vec3f &v2 = mesh->P[mesh->trisIndex[j + 2]];
			float t, u, v;
			stats->numRayTriangleTests++;
			if (scratch::geometry_utils::rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v)) {
				if (t < tHit) {
					tHit = t;
					uhit = u;
					vhit = v;
					intersected = true;
					intersectedTri = triangles[i].tri;
					intersectedMesh = triangles[i].mesh;
					stats->numRayTriangleIntersections++;
				}
			}
		} else {
			intersectedMesh = triangles[i].mesh;
		}
	}

	if (intersected) {
		intersectedMesh->mailbox[intersectedTri] = rayId;
	}
	return (intersectedMesh != nullptr);
}

bool Grid::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
{
	const Vec3f invDir = 1 / dir;
	const Vec3b sign(dir.x < 0, dir.y < 0, dir.z < 0);
	float tHitBox;
	if (!bbox.intersect(orig, invDir, sign, tHitBox)) return false;

	// initialization step
	Vec3i exit, step, cell;
	Vec3f deltaT, nextCrossingT;
	for (uint8_t i = 0; i < 3; ++i) {
		// convert ray starting point to cell coordinates
		float rayOrigCell = ((orig[i] + dir[i] * tHitBox) -  bbox[0][i]);
		cell[i] = scratch::utils::clamp<uint32_t>(std::floor(rayOrigCell / cellDimension[i]), 0, resolution[i] - 1);
		if (dir[i] < 0) {
			deltaT[i] = -cellDimension[i] * invDir[i];
			nextCrossingT[i] = tHitBox + (cell[i] * cellDimension[i] - rayOrigCell) * invDir[i];
			exit[i] = -1;
			step[i] = -1;
		}
		else {
			deltaT[i] = cellDimension[i] * invDir[i];
			nextCrossingT[i] = tHitBox + ((cell[i] + 1)  * cellDimension[i] - rayOrigCell) * invDir[i];
			exit[i] = resolution[i];
			step[i] = 1;
		}
	}

	// Walk through each cell of the grid and test for an intersection if
	// current cell contains geometry
	const TriangleMesh* intersectedMesh = nullptr;
	while (1) {
		uint32_t o = cell[2] * resolution[0] * resolution[1] + cell[1] * resolution[0] + cell[0];
		if (cells[o] != nullptr) {
			cells[o]->intersect(orig, dir, rayId, tHit, intersectedMesh);
			//if (intersectedMesh != nullptr) { ray.color = cells[o]->color; }
		}
		uint8_t k =
			((nextCrossingT[0] < nextCrossingT[1]) << 2) +
			((nextCrossingT[0] < nextCrossingT[2]) << 1) +
			((nextCrossingT[1] < nextCrossingT[2]));
		static const uint8_t map[8] = {2, 1, 2, 1, 2, 2, 0, 0};
		uint8_t axis = map[k];

		if (tHit < nextCrossingT[axis]) break;
		cell[axis] += step[axis];
		if (cell[axis] == exit[axis]) break;
		nextCrossingT[axis] += deltaT[axis];
	}

	return (intersectedMesh != nullptr);
}
