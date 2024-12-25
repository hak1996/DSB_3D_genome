import numpy as np
from tqdm import tqdm
import open3d as o3d
import pandas as pd
from ReadGTrack import GtrackFile


class TADInfo:
    def __init__(self, vertex_id, x, y, z, radius, bp):
        self.id = vertex_id
        self.center = (x, y, z)
        self.radius = radius
        self.bp = bp
        self.neighbors = set()  # 存储相邻顶点的集合

class Sphere:
    def __init__(self, id, center, radius):
        self.center = np.array(center, dtype=np.float32)
        self.radius = radius
        self.id = id

    def get_aabb(self):
        return AABB(self.center - self.radius, self.center + self.radius)


class AABB:
    def __init__(self, min_point, max_point):
        self.min_point = np.array(min_point, dtype=np.float32)
        self.max_point = np.array(max_point, dtype=np.float32)

    def contains_point(self, point):
        return np.all(self.min_point <= point) and np.all(point <= self.max_point)

    def intersects_aabb(self, other):
        return np.all(self.min_point <= other.max_point) and np.all(other.min_point <= self.max_point)

class BVHNode:
    def __init__(self, spheres):
        self.spheres = spheres
        self.left = None
        self.right = None
        self.aabb = self._compute_aabb()

    def _compute_aabb(self):
        min_point = np.min([sphere.center - sphere.radius for sphere in self.spheres], axis=0)
        max_point = np.max([sphere.center + sphere.radius for sphere in self.spheres], axis=0)
        return AABB(min_point, max_point)

    def build(self):
        if len(self.spheres) > 1:
            # Split spheres into two groups based on the center of the AABB
            axis = np.argmax(self.aabb.max_point - self.aabb.min_point)
            sorted_spheres = sorted(self.spheres, key=lambda s: s.center[axis])
            mid = len(sorted_spheres) // 2
            self.left = BVHNode(sorted_spheres[:mid])
            self.right = BVHNode(sorted_spheres[mid:])
            self.left.build()
            self.right.build()
        else:
            # Leaf node
            pass

    def query(self, point):
        if not self.aabb.contains_point(point):
            return None
        if self.left is None and self.right is None:
            # Leaf node
            for sphere in self.spheres:
                if np.linalg.norm(point - sphere.center) <= sphere.radius:
                    return sphere
        else:
            # Internal node
            sphere = self.left.query(point)
            if sphere is not None:
                return sphere
            return self.right.query(point)


class NucleusModel():
    def __init__(self):
        self.TADs = {}
        self.spheres = []

    def AddTAD(self, id, x, y, z, r, bp):
        if id not in self.TADs:
            self.TADs[id] = TADInfo(id, x, y, z, r, bp)

    def AddEdge(self, id1, id2):
        self.TADs[id1].neighbors.add(id2)
        self.TADs[id2].neighbors.add(id1)

    def AddSphere(self, id, x, y, z, r):
        bead = Sphere(id, [x, y, z], r)
        self.spheres.append(bead)

    def ConstructTree(self):
        self.root = BVHNode(self.spheres)
        self.root.build()

    def get_neighbors(self, vertex_id):
        # 获取与指定顶点相连的所有顶点
        if vertex_id in self.TADs:
            return self.TADs[vertex_id].neighbors
        else:
            return None  # 如果顶点不存在，则返回 None

    def IsPointsInTAD(self, points):
        point_cloud = o3d.geometry.PointCloud()
        point_cloud.points = o3d.utility.Vector3dVector(points)
        pcd_tree = o3d.geometry.KDTreeFlann(point_cloud)
        mask = np.zeros(len(points), dtype=bool)
        for TAD in self.TADs.values():
            xyz = TAD.center
            radius = TAD.radius
            [k, idx, _] = pcd_tree.search_radius_vector_3d(xyz, radius)
            mask[idx] = True
        return mask

    def InWhichTAD(self, point):
        sphere = self.root.query(point)
        id = None
        if sphere is not None:
            id = sphere.id
        return id

def ConstructNucleus():
    # read gtrack
    gtrack_file = 'IMR90_noLADs.gtrack'
    gtrack = GtrackFile()
    interaction_ids = gtrack.read_gtrack(gtrack_file)

    #read nucleus model
    file_dir = "IMR90_noLADs_1.parquet"
    df = pd.read_parquet(file_dir, engine='pyarrow')
    beadID = df.beadID.tolist()
    bp = df.bp.tolist()
    x = df.x.tolist()
    y = df.y.tolist()
    z = df.z.tolist()
    r = df.radius.tolist()
    radius = [float(item) for item in r]

    nucleus = NucleusModel()
    for i, id in enumerate(beadID):
        nucleus.AddTAD(id, x[i]*1000.0, y[i]*1000.0, z[i]*1000.0, radius[i]*1000.0, bp[i])
        nucleus.AddSphere(id, x[i]*1000.0, y[i]*1000.0, z[i]*1000.0, radius[i]*1000.0)
    nucleus.ConstructTree()

    for edge in interaction_ids:
        nucleus.AddEdge(edge[0], edge[1])

    return nucleus

if __name__ == "__main__":
    nucleus = ConstructNucleus()
    print(len(nucleus.TADs))
    count_have_edge = 0
    for TAD in nucleus.TADs.values():
        if len(TAD.neighbors) > 0:
            count_have_edge += 1
    print(count_have_edge)
